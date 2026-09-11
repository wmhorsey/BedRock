//! Causal-front "cone test": a falsification harness for the finite-speed
//! substrate.
//!
//! The ontology (Axiom 10, SimulationSpec section 3) asserts a finite
//! propagation speed enforced by a per-step causal horizon
//! `ell_causal <= v_eff * dt`. This harness injects a single localized tension
//! pulse into an otherwise quiescent uniform substrate and measures the outward
//! disturbance front over time.
//!
//! To separate the pulse-induced signal from static lattice noise (the jittered
//! sample layout produces a small background force everywhere), the harness runs
//! a **control** simulation in lockstep: identical seed and initial layout, but
//! no pulse. Every measured quantity is the difference between the perturbed and
//! control runs, so only the causal signal remains.
//!
//! A finite-speed substrate must trace a causal cone: the disturbance leading
//! edge grows at most linearly, never exceeds `perturb_radius + step *
//! causal_horizon`, and is isotropic. The cone speed should be governed by the
//! physical `signal_speed`, not by the cosmetic tuning gains — that invariance
//! is what turns the pile of knobs into testable physics.

use std::error::Error;
use std::f32::consts::{PI, TAU};
use std::fs::File;
use std::io::{BufWriter, Write};

use nalgebra::Vector2;

use crate::config::Config;
use crate::physics::{self, ParticleKind, Simulation};

/// One measured time slice of the propagating front (perturbed minus control).
struct Sample {
    step: usize,
    time: f32,
    /// Outer radius where the force difference exceeds the strong-signal
    /// threshold (the detectable amplitude front; diffusion-limited).
    signal_front_max: f32,
    /// Median per-bin strong-signal front radius.
    signal_front_median: f32,
    /// Outer radius where the tension difference exceeds its threshold.
    tension_front: f32,
    /// Outer radius of any force difference above the tiny floor: the causal
    /// leading edge that must respect the cone.
    leading_edge: f32,
    /// Hard causal bound `perturb_radius + step * causal_horizon`.
    causal_radius: f32,
    /// Coefficient of variation of the per-bin leading edge (0 == round cone).
    isotropy_cv: f32,
    disturbed: usize,
}

struct ConeResult {
    samples: Vec<Sample>,
    horizon: f32,
    spacing: f32,
    /// Worst amount the leading edge exceeded `causal_radius + spacing`.
    worst_overshoot: f32,
    /// Last step whose leading edge stayed clear of the boundary buffer.
    last_clean_step: usize,
}

pub fn run(config: Config) -> Result<(), Box<dyn Error>> {
    let result = measure(&config);
    write_csv(&config.cone_csv, &result.samples)?;
    report(&config, &result);
    Ok(())
}

fn measure(config: &Config) -> ConeResult {
    let (w, h) = config.world_size();
    let center = Vector2::new(w * 0.5, h * 0.5);
    let ambient = config.ambient_tension;
    let perturb_radius = config.cone_perturb_radius.max(1.0);
    let amplitude = config.cone_perturb_amplitude;
    let bins = config.cone_bins.max(1);
    let dt = config.dt;
    let horizon = physics::causal_horizon(config);
    let spacing = physics::substrate_spacing(config).max(0.5);

    // Two identically seeded runs. `perturbed` carries the pulse; `control` does
    // not. Both start from the same quiescent uniform substrate so their
    // difference is exactly the causal signal.
    let mut perturbed = Simulation::new(config.clone());
    let mut control = Simulation::new(config.clone());
    make_quiescent(&mut perturbed, ambient);
    make_quiescent(&mut control, ambient);
    inject_pulse(&mut perturbed, center, perturb_radius, amplitude, ambient);
    perturbed.grid.rebuild(&perturbed.particles);
    control.grid.rebuild(&control.particles);

    let nearest_boundary = center.x.min(w - center.x).min(center.y).min(h - center.y);
    let stop_radius = 0.85 * nearest_boundary;

    // Thresholds calibrated once from the first step's peak force difference so
    // the "disturbed" definition scales with the actual response magnitude.
    let mut signal_threshold = 0.0f32;
    let mut edge_floor = 0.0f32;
    let tension_threshold = (0.01 * amplitude).max(1.0e-4);

    let mut samples: Vec<Sample> = Vec::with_capacity(config.cone_steps);
    let mut bin_signal = vec![0.0f32; bins];
    let mut bin_edge = vec![0.0f32; bins];
    let mut worst_overshoot = f32::MIN;
    let mut last_clean_step = 0usize;

    for step in 1..=config.cone_steps {
        perturbed.update();
        control.update();

        if step == 1 {
            let peak = perturbed
                .particles
                .iter()
                .zip(&control.particles)
                .map(|(p, c)| (p.last_force - c.last_force).norm())
                .fold(0.0f32, f32::max);
            signal_threshold = (0.01 * peak).max(1.0e-6);
            edge_floor = (1.0e-4 * peak).max(1.0e-9);
        }

        for slot in bin_signal.iter_mut() {
            *slot = 0.0;
        }
        for slot in bin_edge.iter_mut() {
            *slot = 0.0;
        }
        let mut tension_front = 0.0f32;
        let mut disturbed = 0usize;

        for (p, c) in perturbed.particles.iter().zip(&control.particles) {
            let offset = p.position - center;
            let radius = offset.norm();
            let force_diff = (p.last_force - c.last_force).norm();
            let tension_diff = (p.tension - c.tension).abs();

            let angle = offset.y.atan2(offset.x).rem_euclid(TAU);
            let bin = ((angle / TAU) * bins as f32) as usize % bins;

            if tension_diff > tension_threshold {
                tension_front = tension_front.max(radius);
            }
            if force_diff > edge_floor {
                bin_edge[bin] = bin_edge[bin].max(radius);
            }
            if force_diff > signal_threshold {
                disturbed += 1;
                bin_signal[bin] = bin_signal[bin].max(radius);
            }
        }

        let signal_front_max = bin_signal.iter().copied().fold(0.0f32, f32::max);
        let signal_front_median = median(&bin_signal);
        let leading_edge = bin_edge.iter().copied().fold(0.0f32, f32::max);
        let isotropy_cv = coefficient_of_variation(&bin_edge);
        let causal_radius = perturb_radius + step as f32 * horizon;
        worst_overshoot = worst_overshoot.max(leading_edge - (causal_radius + spacing));

        samples.push(Sample {
            step,
            time: step as f32 * dt,
            signal_front_max,
            signal_front_median,
            tension_front,
            leading_edge,
            causal_radius,
            isotropy_cv,
            disturbed,
        });

        if leading_edge < stop_radius {
            last_clean_step = step;
        } else {
            // The cone is nearing the domain edge; stop before boundary
            // reflections contaminate the measurement.
            break;
        }
    }

    if last_clean_step == 0 {
        last_clean_step = samples.len();
    }

    ConeResult {
        samples,
        horizon,
        spacing,
        worst_overshoot,
        last_clean_step,
    }
}

/// Overwrite a freshly seeded simulation into a clean, quiescent uniform
/// substrate: no cavities, ambient tension everywhere, zero velocity, cleared
/// derived state. Sample positions (including their fixed jitter) are kept so
/// the control and perturbed runs share an identical layout.
fn make_quiescent(sim: &mut Simulation, ambient: f32) {
    sim.cavities.clear();
    for particle in &mut sim.particles {
        particle.kind = ParticleKind::Substrate;
        particle.fixed = false;
        particle.velocity = Vector2::zeros();
        particle.last_force = Vector2::zeros();
        particle.tension = ambient;
        particle.base_tension = ambient;
        particle.ejecta_timer = 0.0;
        particle.local_density = 0.0;
        particle.shell_reactivity = 0.0;
        particle.merge_factor = 0.0;
        particle.influence_scale = 1.0;
        particle.stored_energy = 0.0;
        particle.shell_overlap = 0.0;
        particle.spike_drive = 0.0;
    }
}

/// Add a single compactly-supported cos^2 tension pulse centered at `center`.
/// Compact support keeps the causal bound honest (the initial disturbance is
/// strictly inside `radius`); the soft edge keeps the source isotropic so the
/// isotropy check reflects the dynamics, not a lattice-aligned discontinuity.
fn inject_pulse(
    sim: &mut Simulation,
    center: Vector2<f32>,
    radius: f32,
    amplitude: f32,
    ambient: f32,
) {
    for particle in &mut sim.particles {
        let distance = (particle.position - center).norm();
        if distance <= radius {
            let phase = 0.5 * PI * distance / radius;
            let bump = amplitude * phase.cos() * phase.cos();
            particle.tension = ambient + bump;
            particle.base_tension = particle.tension;
        }
    }
}

fn write_csv(path: &str, samples: &[Sample]) -> Result<(), Box<dyn Error>> {
    if path.is_empty() {
        return Ok(());
    }
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "step,time,leading_edge,causal_radius,signal_front_max,signal_front_median,tension_front,isotropy_cv,disturbed"
    )?;
    for s in samples {
        writeln!(
            writer,
            "{},{:.5},{:.4},{:.4},{:.4},{:.4},{:.4},{:.5},{}",
            s.step,
            s.time,
            s.leading_edge,
            s.causal_radius,
            s.signal_front_max,
            s.signal_front_median,
            s.tension_front,
            s.isotropy_cv,
            s.disturbed,
        )?;
    }
    Ok(())
}

fn report(config: &Config, result: &ConeResult) {
    let dt = config.dt;
    let horizon = result.horizon;
    let spacing = result.spacing;
    let (w, h) = config.world_size();

    let phys = config.signal_speed * dt;
    let inner = phys.min(config.interaction_radius);
    let classification = if spacing >= inner {
        "spacing-limited: floored by substrate spacing (raise particle_count or signal_speed*dt)"
    } else if phys <= config.interaction_radius {
        "physical: set by signal_speed * dt"
    } else {
        "NUMERICS-LIMITED: clamped by interaction_radius (raise interaction_radius or lower signal_speed*dt)"
    };

    // Linear-fit the fronts over the clean window, skipping the source-size
    // startup transient.
    let warmup = 4usize;
    let clean: Vec<&Sample> = result
        .samples
        .iter()
        .filter(|s| s.step <= result.last_clean_step)
        .collect();
    let slope_edge = fit_slope(&clean, warmup, |s| s.leading_edge);
    let slope_signal = fit_slope(&clean, warmup, |s| s.signal_front_max);

    let final_cv = clean.last().map(|s| s.isotropy_cv).unwrap_or(0.0);

    let causality = if result.worst_overshoot <= 0.0 {
        "PASS"
    } else if result.worst_overshoot <= spacing {
        "WARN"
    } else {
        "FAIL"
    };
    let isotropy = if final_cv < 0.15 {
        "PASS"
    } else if final_cv < 0.30 {
        "WARN"
    } else {
        "FAIL"
    };

    println!();
    println!("=== BedRock Cone Test (causal-front falsification) ===");
    println!("domain            : {:.0} x {:.0}", w, h);
    println!("particles         : {}", config.particle_count);
    println!("substrate spacing : {:.3}", spacing);
    println!("dt                : {:.4}", dt);
    println!(
        "signal_speed      : {:.2}  (signal_speed*dt = {:.3})",
        config.signal_speed, phys
    );
    println!("interaction_radius: {:.2}", config.interaction_radius);
    println!(
        "causal horizon    : {:.3} units/step  [{}]",
        horizon, classification
    );
    println!(
        "pulse             : radius {:.2}, amplitude {:.2} over ambient {:.2}",
        config.cone_perturb_radius, config.cone_perturb_amplitude, config.ambient_tension
    );
    println!(
        "steps measured    : {} (clean window ends at step {})",
        result.samples.len(),
        result.last_clean_step
    );
    println!("-----------------------------------------------------");
    println!(
        "cone edge speed    : {:.3} units/step = {:.2} units/time",
        slope_edge,
        slope_edge / dt
    );
    println!(
        "signal front speed : {:.3} units/step = {:.2} units/time",
        slope_signal,
        slope_signal / dt
    );
    println!(
        "causal limit       : {:.3} units/step = {:.2} units/time",
        horizon,
        horizon / dt
    );
    if horizon > 0.0 {
        println!("edge / causal ratio: {:.2}", slope_edge / horizon);
    }
    println!("-----------------------------------------------------");
    println!(
        "causality (no superluminal leakage): {}  (worst overshoot {:.3} units beyond causal+spacing)",
        causality, result.worst_overshoot
    );
    println!(
        "isotropy  (cone CV at step {})      : {:.3}  {}",
        result.last_clean_step, final_cv, isotropy
    );
    if !config.cone_csv.is_empty() {
        println!("per-step CSV written to            : {}", config.cone_csv);
    }
    println!("-----------------------------------------------------");
    println!("Interpretation:");
    println!("  The causal limit is fixed by signal_speed, a physical parameter.");
    println!("  To promote this from a knob to a law, sweep the cosmetic gains");
    println!("  (attraction_strength, shell_overlap_gain, spike_threshold, ...)");
    println!("  and confirm the cone edge speed and the causality verdict stay");
    println!("  invariant.");
}

fn fit_slope(clean: &[&Sample], warmup: usize, value: impl Fn(&Sample) -> f32) -> f32 {
    let mut points: Vec<(f32, f32)> = clean
        .iter()
        .filter(|s| s.step >= warmup)
        .map(|s| (s.step as f32, value(s)))
        .collect();
    if points.len() < 3 {
        points = clean.iter().map(|s| (s.step as f32, value(s))).collect();
    }
    linear_slope(&points)
}

fn linear_slope(points: &[(f32, f32)]) -> f32 {
    let n = points.len() as f32;
    if n < 2.0 {
        return 0.0;
    }
    let sum_x: f32 = points.iter().map(|p| p.0).sum();
    let sum_y: f32 = points.iter().map(|p| p.1).sum();
    let sum_xx: f32 = points.iter().map(|p| p.0 * p.0).sum();
    let sum_xy: f32 = points.iter().map(|p| p.0 * p.1).sum();
    let denom = n * sum_xx - sum_x * sum_x;
    if denom.abs() < 1.0e-9 {
        return 0.0;
    }
    (n * sum_xy - sum_x * sum_y) / denom
}

fn median(values: &[f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        0.5 * (sorted[mid - 1] + sorted[mid])
    } else {
        sorted[mid]
    }
}

fn coefficient_of_variation(values: &[f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let n = values.len() as f32;
    let mean = values.iter().sum::<f32>() / n;
    if mean <= 1.0e-9 {
        return 0.0;
    }
    let variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f32>() / n;
    variance.sqrt() / mean
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    fn cone_config() -> Config {
        let mut config = Config::parse_from(["bedrock"]);
        config.mode = crate::config::SimulationMode::TwoD;
        config.width = 400;
        config.height = 400;
        config.particle_count = 8000;
        // spacing = sqrt(400*400/8000) = 4.47 ; signal_speed*dt = 90*0.08 = 7.2 ;
        // interaction_radius = 16 -> causal horizon = 7.2 (physical).
        config.signal_speed = 90.0;
        config.cone_steps = 30;
        config.cone_perturb_radius = 10.0;
        config.cone_perturb_amplitude = 1.5;
        config.cone_bins = 16;
        config.cone_csv = String::new();
        config
    }

    #[test]
    fn horizon_is_physical_for_cone_config() {
        let config = cone_config();
        let horizon = physics::causal_horizon(&config);
        assert!((horizon - config.signal_speed * config.dt).abs() < 1.0e-3);
    }

    #[test]
    fn front_respects_causal_bound() {
        let result = measure(&cone_config());
        // No disturbance may appear beyond the causal cone (allowing one lattice
        // spacing of discretization slack on top of the built-in tolerance).
        assert!(
            result.worst_overshoot <= result.spacing,
            "superluminal leakage: worst overshoot {} exceeds spacing {}",
            result.worst_overshoot,
            result.spacing
        );
    }

    #[test]
    fn cone_edge_advances_but_stays_subluminal() {
        let result = measure(&cone_config());
        let horizon = result.horizon;
        let clean: Vec<&Sample> = result
            .samples
            .iter()
            .filter(|s| s.step <= result.last_clean_step)
            .collect();
        let slope = fit_slope(&clean, 4, |s| s.leading_edge);
        // The cone edge must actually propagate outward...
        assert!(slope > 0.0, "cone edge did not advance: slope {}", slope);
        // ...but never faster than the causal horizon per step.
        assert!(
            slope <= horizon + result.spacing,
            "cone edge speed {} exceeds causal horizon {}",
            slope,
            horizon
        );
    }
}
