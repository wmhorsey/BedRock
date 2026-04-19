use std::f32::consts::TAU;

use nalgebra::Vector2;

use crate::{
    config::Config,
    neighbor_search::SpatialHashGrid,
    physics::{Cavity, Particle, ParticleKind},
};

#[derive(Debug, Clone)]
pub struct BreachEvent {
    pub cavity_id: usize,
    pub direction: Vector2<f32>,
    pub strength: f32,
}

pub fn detect_breaches(
    particles: &[Particle],
    cavities: &[Cavity],
    grid: &SpatialHashGrid,
    config: &Config,
) -> Vec<BreachEvent> {
    let mut events = Vec::new();

    for cavity in cavities {
        if cavity.cooldown_steps > 0 {
            continue;
        }

        let mut shell_indices = Vec::new();
        grid.query_neighbors(cavity.center, &mut shell_indices);

        let mut internal = 0.0;
        let mut external = 0.0;
        let mut internal_count = 0.0;
        let mut external_count = 0.0;
        let mut angle_bins = [0.0f32; 24];

        for index in shell_indices {
            let particle = &particles[index];
            if particle.kind == ParticleKind::TrueVoid {
                continue;
            }

            let offset = particle.position - cavity.center;
            let distance = offset.norm();
            if distance <= f32::EPSILON {
                continue;
            }

            let inward_force = particle.last_force.dot(&(-offset / distance)).max(0.0);
            let crowding_term = particle.shell_reactivity + (particle.local_density - 1.0).max(0.0) * 0.35;

            if distance >= cavity.radius && distance <= cavity.radius + cavity.shell_thickness {
                internal += inward_force + crowding_term;
                internal_count += 1.0;

                let angle = offset.y.atan2(offset.x).rem_euclid(TAU);
                let bin = ((angle / TAU) * angle_bins.len() as f32) as usize % angle_bins.len();
                angle_bins[bin] += inward_force
                    + (particle.tension - config.ambient_tension).max(0.0)
                    + crowding_term;
            } else if distance <= cavity.radius + cavity.shell_thickness * 3.0 {
                external += inward_force + crowding_term * 0.2;
                external_count += 1.0;
            }
        }

        if internal_count < 6.0 || external_count < 6.0 {
            continue;
        }

        internal /= internal_count;
        external = (external / external_count).max(0.0001);

        let curvature = 1.0 / cavity.radius.max(1.0);
        let core_depth = if cavity.core_tension == 0.0 {
            0.0
        } else {
            (config.ambient_tension - cavity.core_tension).max(0.0)
        };

        let ratio_factor = config.breach_threshold_base
            * (1.0 + 0.3 * config.ambient_tension + 0.2 * core_depth)
            / (1.0 + config.curvature_scale * curvature * 12.0);

        if internal <= external * ratio_factor {
            continue;
        }

        let (dominant_bin, dominant_value) = angle_bins
            .iter()
            .copied()
            .enumerate()
            .max_by(|(_, left), (_, right)| left.partial_cmp(right).unwrap())
            .unwrap();

        if dominant_value <= 0.0 {
            continue;
        }

        let angle = TAU * dominant_bin as f32 / angle_bins.len() as f32;
        let direction = Vector2::new(angle.cos(), angle.sin());
        let strength = ((internal / external) - ratio_factor).max(0.0);

        events.push(BreachEvent {
            cavity_id: cavity.id,
            direction,
            strength,
        });
    }

    events
}