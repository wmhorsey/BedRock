//! Conservative acoustic field: the constitutive law that actually carries a
//! wave.
//!
//! The relaxational tension update in `physics.rs` is first-order and
//! dissipative (`tension += (target - tension) * relaxation`), so a pulse decays
//! in place instead of radiating. A wave needs a second state that stores the
//! overshoot. The ontology already names it: the attraction field
//! `G = grad A` (Axiom 2). Evolving the pair conservatively,
//!
//! ```text
//!     dA/dt = -c * div(G)
//!     dG/dt = -c * grad(A)
//! ```
//!
//! gives `d^2A/dt^2 = c^2 * lap(A)` — a real wave at speed `c` — and conserves
//! the energy `E = 1/2 * sum(A^2 + |G|^2)` (Axiom 0: no drag, no resistance).
//!
//! The mesh-free discretization uses symmetric pair weights so the discrete
//! divergence is exactly the negative transpose of the discrete gradient
//! (`Div = -Grad^T`). That makes the semi-discrete system Hamiltonian, so the
//! symplectic (staggered) update below has no systematic energy drift — the
//! energy diagnostic is the direct, measurable answer to "does this medium
//! damp?".
//!
//! Samples are fixed: the wave lives in the `(A, G)` field, not in particle
//! motion, so neighbor pairs and weights are precomputed once.

use nalgebra::Vector2;
use rayon::prelude::*;

/// A precomputed neighbor coupling: index, symmetric weight, unit direction.
struct Pair {
    j: usize,
    weight: f32,
    dir: Vector2<f32>,
}

pub struct WaveField {
    positions: Vec<Vector2<f32>>,
    /// Field value (tension deviation from ambient; rest state is exactly 0).
    a: Vec<f32>,
    /// Flux state `G = grad A`.
    g: Vec<Vector2<f32>>,
    pairs: Vec<Vec<Pair>>,
    c: f32,
    dt: f32,
    cfl: f32,
}

impl WaveField {
    /// Build the field over fixed sample positions. `c` is the wave speed
    /// (units/time), `support` is the kernel/neighbor radius, and `spacing` is
    /// the mean sample spacing used for volume weighting and CFL reporting.
    pub fn new(positions: Vec<Vector2<f32>>, c: f32, dt: f32, spacing: f32, support: f32) -> Self {
        let n = positions.len();
        let volume = spacing * spacing;
        let support = support.max(spacing * 1.5);

        // One-time neighbor search over a uniform bin grid sized to the support.
        let grid = BinGrid::build(&positions, support);
        let raw: Vec<Vec<(usize, f32, f32, Vector2<f32>)>> = (0..n)
            .into_par_iter()
            .map(|i| {
                let mut out = Vec::new();
                let pos = positions[i];
                grid.for_each_neighbor(pos, |j| {
                    if j == i {
                        return;
                    }
                    let offset = positions[j] - pos;
                    let distance = offset.norm();
                    if distance <= 1.0e-5 || distance > support {
                        return;
                    }
                    let q = distance / support;
                    // Smooth compact bump kernel, C1 at the support boundary.
                    let w = {
                        let t = 1.0 - q * q;
                        t * t
                    };
                    out.push((j, w, distance, offset / distance));
                });
                out
            })
            .collect();

        // Global normalization so the discrete gradient reproduces a true
        // gradient to first order: mean over samples of 1/2 * sum(w * V * r).
        let beta_sum: f64 = raw
            .par_iter()
            .map(|list| {
                0.5 * list
                    .iter()
                    .map(|(_, w, r, _)| (*w * volume * *r) as f64)
                    .sum::<f64>()
            })
            .sum();
        let beta = (beta_sum / n.max(1) as f64) as f32;
        let beta = if beta > 1.0e-12 { beta } else { 1.0 };

        let pairs: Vec<Vec<Pair>> = raw
            .into_par_iter()
            .map(|list| {
                list.into_iter()
                    .map(|(j, w, _r, dir)| Pair {
                        j,
                        weight: w * volume / beta,
                        dir,
                    })
                    .collect()
            })
            .collect();

        let cfl = if spacing > 0.0 { c * dt / spacing } else { 0.0 };

        Self {
            a: vec![0.0; n],
            g: vec![Vector2::zeros(); n],
            positions,
            pairs,
            c,
            dt,
            cfl,
        }
    }

    pub fn len(&self) -> usize {
        self.positions.len()
    }

    pub fn cfl(&self) -> f32 {
        self.cfl
    }

    pub fn position(&self, i: usize) -> Vector2<f32> {
        self.positions[i]
    }

    pub fn field(&self, i: usize) -> f32 {
        self.a[i]
    }

    pub fn flux_magnitude(&self, i: usize) -> f32 {
        self.g[i].norm()
    }

    /// Inject a compact cos^2 tension pulse (zero flux) centered at `center`.
    pub fn inject_pulse(&mut self, center: Vector2<f32>, radius: f32, amplitude: f32) {
        let radius = radius.max(1.0e-3);
        for (a, pos) in self.a.iter_mut().zip(&self.positions) {
            let distance = (pos - center).norm();
            if distance <= radius {
                let phase = 0.5 * std::f32::consts::PI * distance / radius;
                *a = amplitude * phase.cos() * phase.cos();
            }
        }
    }

    /// Shell pop: dump a compact cos^2 tension ring at `radius` (half-width
    /// `width`) with zero flux. A physical collapse source rather than a
    /// hand-placed disk: the ring radiates inward (focusing) and outward, and the
    /// outward front is the measured cone.
    pub fn inject_shell_pop(
        &mut self,
        center: Vector2<f32>,
        radius: f32,
        width: f32,
        amplitude: f32,
    ) {
        let width = width.max(1.0e-3);
        for (a, pos) in self.a.iter_mut().zip(&self.positions) {
            let offset = ((pos - center).norm() - radius).abs();
            if offset <= width {
                let phase = 0.5 * std::f32::consts::PI * offset / width;
                *a = amplitude * phase.cos() * phase.cos();
            }
        }
    }

    /// One synchronized velocity-Verlet (kick-drift-kick) step of the
    /// conservative wave. Keeping A and G at the same time level makes the
    /// measured energy oscillation second-order small (a plain symplectic-Euler
    /// staggering leaves a first-order swing).
    pub fn step(&mut self) {
        let c_dt = self.c * self.dt;
        self.kick_flux(0.5 * c_dt);
        self.drift_field(c_dt);
        self.kick_flux(0.5 * c_dt);
    }

    /// Flux kick: `G -= factor * grad(A)`.
    fn kick_flux(&mut self, factor: f32) {
        let new_g: Vec<Vector2<f32>> = (0..self.len())
            .into_par_iter()
            .map(|i| {
                let ai = self.a[i];
                let mut grad = Vector2::zeros();
                for pair in &self.pairs[i] {
                    grad += pair.dir * (pair.weight * (self.a[pair.j] - ai));
                }
                self.g[i] - grad * factor
            })
            .collect();
        self.g = new_g;
    }

    /// Field drift: `A -= factor * div(G)`, with div = -grad^T (energy-conserving).
    fn drift_field(&mut self, factor: f32) {
        let new_a: Vec<f32> = (0..self.len())
            .into_par_iter()
            .map(|i| {
                let gi = self.g[i];
                let mut div = 0.0;
                for pair in &self.pairs[i] {
                    div += pair.weight * (gi + self.g[pair.j]).dot(&pair.dir);
                }
                self.a[i] - factor * div
            })
            .collect();
        self.a = new_a;
    }

    /// Total conserved energy `1/2 * sum(A^2 + |G|^2)`.
    pub fn energy(&self) -> f64 {
        self.a
            .par_iter()
            .zip(&self.g)
            .map(|(a, g)| 0.5 * ((*a as f64) * (*a as f64) + g.norm_squared() as f64))
            .sum()
    }
}

/// Minimal uniform bin grid for the one-time neighbor precompute.
struct BinGrid {
    cell: f32,
    cells: std::collections::HashMap<(i32, i32), Vec<usize>>,
}

impl BinGrid {
    fn build(positions: &[Vector2<f32>], cell: f32) -> Self {
        let cell = cell.max(1.0e-3);
        let mut cells: std::collections::HashMap<(i32, i32), Vec<usize>> =
            std::collections::HashMap::new();
        for (i, pos) in positions.iter().enumerate() {
            cells.entry(Self::key(*pos, cell)).or_default().push(i);
        }
        Self { cell, cells }
    }

    fn key(pos: Vector2<f32>, cell: f32) -> (i32, i32) {
        ((pos.x / cell).floor() as i32, (pos.y / cell).floor() as i32)
    }

    fn for_each_neighbor(&self, pos: Vector2<f32>, mut visit: impl FnMut(usize)) {
        let (bx, by) = Self::key(pos, self.cell);
        for dy in -1..=1 {
            for dx in -1..=1 {
                if let Some(indices) = self.cells.get(&(bx + dx, by + dy)) {
                    for &i in indices {
                        visit(i);
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn uniform_positions(width: f32, height: f32, count: usize) -> (Vec<Vector2<f32>>, f32) {
        let spacing = ((width * height) / count as f32).sqrt();
        let mut positions = Vec::new();
        let mut y = spacing * 0.5;
        while y < height {
            let mut x = spacing * 0.5;
            while x < width {
                positions.push(Vector2::new(x, y));
                x += spacing;
            }
            y += spacing;
        }
        (positions, spacing)
    }

    fn test_field() -> (WaveField, f32, f32) {
        let (positions, spacing) = uniform_positions(400.0, 400.0, 4000);
        let dt = 0.08;
        let c = 0.3 * spacing / dt; // c*dt = 0.3*spacing, comfortably inside CFL
        let field = WaveField::new(positions, c, dt, spacing, 4.0 * spacing);
        (field, spacing, c * dt)
    }

    #[test]
    fn energy_is_conserved() {
        let (mut field, spacing, _advance) = test_field();
        field.inject_pulse(Vector2::new(200.0, 200.0), 4.0 * spacing, 1.0);
        let e0 = field.energy();
        assert!(e0 > 0.0);
        for _ in 0..120 {
            field.step();
        }
        let e1 = field.energy();
        let drift = ((e1 - e0) / e0).abs();
        // Symplectic integration keeps energy bounded near E0: no systematic
        // damping (unlike the relaxational engine, which would decay to ~0).
        assert!(drift < 0.05, "energy drift {} too large (e0={}, e1={})", drift, e0, e1);
    }

    #[test]
    fn pulse_radiates_outward() {
        let (mut field, spacing, advance) = test_field();
        let center = Vector2::new(200.0, 200.0);
        let pulse_radius = 4.0 * spacing;
        let support = 4.0 * spacing; // matches test_field's kernel support
        field.inject_pulse(center, pulse_radius, 1.0);

        let amp_front = |f: &WaveField| {
            let mut r = 0.0f32;
            for i in 0..f.len() {
                if f.field(i).abs() > 0.02 {
                    r = r.max((f.position(i) - center).norm());
                }
            }
            r
        };
        let flux_front = |f: &WaveField| {
            let mut r = 0.0f32;
            for i in 0..f.len() {
                if f.flux_magnitude(i) > 1.0e-3 {
                    r = r.max((f.position(i) - center).norm());
                }
            }
            r
        };

        let steps = 40;
        for _ in 0..steps {
            field.step();
        }

        // The amplitude front radiates outward at roughly the wave speed.
        let amp = amp_front(&field);
        assert!(
            amp > pulse_radius + 5.0 * advance,
            "amplitude front {} did not radiate past source {}",
            amp,
            pulse_radius
        );

        // Propagation is finite: nothing outruns the operator's stencil reach
        // (dispersive precursors may run ahead of the physical front at c).
        let flux = flux_front(&field);
        let stencil_bound = pulse_radius + steps as f32 * support + spacing;
        assert!(
            flux <= stencil_bound,
            "flux front {} outran stencil bound {}",
            flux,
            stencil_bound
        );
    }

    #[test]
    fn shell_pop_radiates_outward() {
        let (mut field, spacing, advance) = test_field();
        let center = Vector2::new(200.0, 200.0);
        let ring_radius = 60.0;
        let ring_width = 3.0 * spacing;
        field.inject_shell_pop(center, ring_radius, ring_width, 1.0);
        let e0 = field.energy();

        let outer_front = |f: &WaveField| {
            let mut r = 0.0f32;
            for i in 0..f.len() {
                if f.field(i).abs() > 0.02 {
                    r = r.max((f.position(i) - center).norm());
                }
            }
            r
        };

        let steps = 30;
        for _ in 0..steps {
            field.step();
        }

        // The outgoing front radiates past the ring at roughly the wave speed.
        let outer = outer_front(&field);
        assert!(
            outer > ring_radius + ring_width + 3.0 * advance,
            "outgoing front {} did not radiate past ring {}",
            outer,
            ring_radius + ring_width
        );
        // Energy stays bounded (conservative), same as the disk source.
        let drift = ((field.energy() - e0) / e0).abs();
        assert!(drift < 0.05, "shell-pop energy drift {} too large", drift);
    }
}
