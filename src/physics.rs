use std::collections::BTreeMap;
use std::f32::consts::PI;

use nalgebra::Vector2;
use rand::{rngs::StdRng, Rng, SeedableRng};
use rayon::prelude::*;

use crate::{
    breach_detection::{detect_breaches, BreachEvent},
    config::Config,
    neighbor_search::SpatialHashGrid,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParticleKind {
    Substrate,
    TrueVoid,
    #[allow(dead_code)]
    Depression,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum QuarkClass {
    Up,
    Charm,
    Top,
}

impl QuarkClass {
    pub fn label(self) -> &'static str {
        match self {
            Self::Up => "up",
            Self::Charm => "charm",
            Self::Top => "top",
        }
    }
}

#[derive(Debug, Clone)]
pub struct Cavity {
    pub id: usize,
    pub center: Vector2<f32>,
    pub radius: f32,
    pub shell_thickness: f32,
    pub core_tension: f32,
    pub class: QuarkClass,
    pub pop_count: u32,
    pub cooldown_steps: u32,
}

#[derive(Debug, Clone)]
pub struct Particle {
    pub position: Vector2<f32>,
    pub velocity: Vector2<f32>,
    pub tension: f32,
    pub kind: ParticleKind,
    pub fixed: bool,
    pub last_force: Vector2<f32>,
    pub base_tension: f32,
    pub ejecta_timer: f32,
    pub local_density: f32,
    pub shell_reactivity: f32,
    pub merge_factor: f32,
    pub influence_scale: f32,
    pub stored_energy: f32,
    pub shell_overlap: f32,
    pub spike_drive: f32,
}

#[derive(Debug, Clone)]
pub struct SimulationMetrics {
    pub total_tension: f32,
    pub true_voids: usize,
    pub pops_this_frame: usize,
    pub shell_radii: BTreeMap<QuarkClass, f32>,
}

pub struct Simulation {
    pub config: Config,
    pub particles: Vec<Particle>,
    pub cavities: Vec<Cavity>,
    pub grid: SpatialHashGrid,
    pub step: u64,
    last_frame_pops: usize,
}

impl Simulation {
    pub fn new(config: Config) -> Self {
        let cavities = seeded_cavities(&config);
        let particles = seeded_particles(&config, &cavities);
        let mut grid = SpatialHashGrid::new(config.cell_size);
        grid.rebuild(&particles);

        Self {
            config,
            particles,
            cavities,
            grid,
            step: 0,
            last_frame_pops: 0,
        }
    }

    pub fn update(&mut self) -> SimulationMetrics {
        self.last_frame_pops = 0;

        for cavity in &mut self.cavities {
            if cavity.cooldown_steps > 0 {
                cavity.cooldown_steps -= 1;
            }
        }

        self.grid.rebuild(&self.particles);
        self.apply_force_step();

        let events = detect_breaches(&self.particles, &self.cavities, &self.grid, &self.config);
        self.last_frame_pops = events.len();
        if !events.is_empty() {
            self.apply_breach_events(&events);
        }

        self.step += 1;
        self.collect_metrics()
    }

    fn apply_force_step(&mut self) {
        let interaction_radius_sq = self.config.interaction_radius * self.config.interaction_radius;
        let dt = self.config.dt;
        let ambient = self.config.ambient_tension;
        let tension_relaxation = self.config.tension_relaxation;
        let attraction_strength = self.config.attraction_strength;
        let density_reactivity = self.config.density_reactivity;
        let density_force_gain = self.config.density_force_gain;
        let merge_density_threshold = self.config.merge_density_threshold;
        let merge_influence_gain = self.config.merge_influence_gain;
        let attraction_ramp_gain = self.config.attraction_ramp_gain;
        let micro_pop_threshold = self.config.micro_pop_threshold;
        let micro_pop_release = self.config.micro_pop_release;
        let flow_coupling = self.config.flow_coupling;
        let collapse_gain = self.config.collapse_gain;
        let shell_source_scale = self.config.shell_source_scale;
        let shell_overlap_gain = self.config.shell_overlap_gain;
        let spike_threshold = self.config.spike_threshold;
        let shell_reform_gain = self.config.shell_reform_gain;
        let void_geometry_gain = self.config.void_geometry_gain;
        let max_speed = self.config.max_speed;
        let velocity_drag = 1.0 - self.config.velocity_damping * dt;
        let particles = self.particles.clone();
        let cavities = self.cavities.clone();
        let grid = &self.grid;
        let expected_neighbors = expected_neighbor_count(&self.config).max(1.0);
        let base_spacing = substrate_spacing(&self.config).max(0.5);
        let causal_horizon = crate::physics::causal_horizon(&self.config);
        let min_shell_gap = base_spacing * 0.35;

        self.particles = (0..particles.len())
            .into_par_iter()
            .map(|index| {
                let mut particle = particles[index].clone();

                if particle.fixed {
                    particle.tension = 0.0;
                    particle.velocity = Vector2::zeros();
                    particle.last_force = Vector2::zeros();
                    particle.local_density = 0.0;
                    particle.shell_reactivity = 0.0;
                    particle.merge_factor = 0.0;
                    particle.influence_scale = 1.0;
                    particle.stored_energy = 0.0;
                    particle.shell_overlap = 0.0;
                    particle.spike_drive = 0.0;
                    return particle;
                }

                let mut neighbor_indices = Vec::with_capacity(32);
                grid.query_neighbors(particle.position, &mut neighbor_indices);

                let mut force = Vector2::zeros();
                let mut neighbor_tension_sum = 0.0;
                let mut neighbor_count = 0.0;
                let mut neighbor_velocity_sum = Vector2::zeros();
                let mut neighbor_position_sum = Vector2::zeros();
                let mut shell_overlap_sum = 0.0;

                for neighbor_index in neighbor_indices {
                    if neighbor_index == index {
                        continue;
                    }

                    let neighbor = &particles[neighbor_index];
                    let offset = neighbor.position - particle.position;
                    let distance_sq = offset.norm_squared();
                    if distance_sq <= 0.0001 || distance_sq > interaction_radius_sq {
                        continue;
                    }

                    let distance = distance_sq.sqrt();
                    if distance > causal_horizon {
                        continue;
                    }

                    let direction = offset / distance;
                    let tension_delta = (neighbor.tension - particle.tension).max(0.0);
                    let shared_influence = 0.5 * (particle.influence_scale + neighbor.influence_scale);
                    let shell_radius = base_spacing * shell_source_scale * shared_influence;
                    let overlap_ratio = ((shell_radius - distance) / shell_radius.max(min_shell_gap)).max(0.0);
                    let shell_gap = (distance - shell_radius).max(min_shell_gap);
                    let inverse_square = 1.0 / (shell_gap * shell_gap + min_shell_gap * min_shell_gap);
                    let causal_weight = 1.0 - (distance / causal_horizon);
                    let overlap_damping = 1.0 / (1.0 + shell_overlap_gain * overlap_ratio);
                    force += direction
                        * attraction_strength
                        * tension_delta
                        * inverse_square
                        * causal_weight.max(0.0)
                        * overlap_damping;

                    neighbor_tension_sum += neighbor.tension;
                    neighbor_count += 1.0;
                    neighbor_velocity_sum += neighbor.velocity;
                    neighbor_position_sum += neighbor.position;
                    shell_overlap_sum += overlap_ratio;
                }

                let (void_tension_boost, void_force, shell_proximity) =
                    void_boundary_response(particle.position, &cavities, ambient, self.config.shell_gain, void_geometry_gain, self.config.void_shell_band_scale);
                let cavity_target = ambient + void_tension_boost;
                let local_mean = if neighbor_count > 0.0 {
                    neighbor_tension_sum / neighbor_count
                } else {
                    ambient
                };
                let local_density = neighbor_count / expected_neighbors;
                let crowding = (local_density - 1.0).max(0.0);
                let shell_reactivity = crowding * shell_proximity;
                let merge_factor = ((local_density - merge_density_threshold) / merge_density_threshold.max(0.1)).max(0.0);
                let shell_overlap = if neighbor_count > 0.0 {
                    shell_overlap_sum / neighbor_count
                } else {
                    0.0
                };
                let spike_drive = ((particle.tension - local_mean) / ambient.max(0.1)).max(0.0)
                    + 0.6 * shell_overlap
                    + 0.25 * merge_factor;
                let spike_excess = (spike_drive - spike_threshold).max(0.0);
                let influence_scale = ((1.0 + merge_influence_gain * merge_factor).sqrt()
                    / (1.0 + shell_overlap_gain * shell_overlap + 0.6 * spike_excess))
                    .clamp(0.85, 3.2);

                force += void_force * (1.0 + 0.35 * crowding + 0.25 * shell_reactivity);

                if neighbor_count > 0.0 {
                    let local_velocity = neighbor_velocity_sum / neighbor_count;
                    let local_centroid = neighbor_position_sum / neighbor_count;
                    let collapse_offset = local_centroid - particle.position;
                    force += collapse_offset * collapse_gain * influence_scale * (0.2 + shell_proximity) * crowding.max(0.05);
                    particle.velocity += (local_velocity - particle.velocity) * flow_coupling * dt;
                }

                force *= (1.0 + density_force_gain * shell_reactivity + 0.22 * (influence_scale - 1.0))
                    * (1.0 - 0.35 * shell_overlap).max(0.25);

                let available_supply = (local_mean - ambient).max(0.0)
                    + 0.45 * crowding
                    + 0.2 * particle.stored_energy
                    + 0.18 * merge_factor;

                let surface_release = shell_reform_gain * (shell_overlap + spike_excess) * (0.35 + shell_proximity);

                let raw_target = cavity_target
                    + 0.08 * (local_mean - ambient)
                    + density_reactivity * shell_reactivity
                    + 0.16 * merge_factor * (0.4 + shell_proximity)
                    - surface_release
                    + particle.ejecta_timer.max(0.0) * 0.05;

                let max_target_step = attraction_ramp_gain * (0.2 + available_supply) * dt * 10.0;
                let target_tension = particle.tension + (raw_target - particle.tension).clamp(-max_target_step, max_target_step);

                particle.tension += (target_tension - particle.tension) * tension_relaxation * dt * 10.0;
                particle.tension = particle.tension.clamp(0.0, ambient * 3.2);

                particle.velocity += force * dt;
                particle.velocity *= velocity_drag.max(0.0);

                let stored_energy = (particle.stored_energy
                    + (merge_factor * (0.32 + shell_proximity) + 0.18 * crowding + 0.4 * surface_release) * dt
                    - micro_pop_release * 0.18 * dt)
                    .max(0.0);
                particle.stored_energy = stored_energy;

                if surface_release > 0.02 && neighbor_count > 0.0 {
                    let local_centroid = neighbor_position_sum / neighbor_count;
                    let outward = particle.position - local_centroid;
                    if outward.norm_squared() > 0.0001 {
                        let outward_dir = outward.normalize();
                        let tangent = Vector2::new(-outward_dir.y, outward_dir.x);
                        let phase = ((self.step as f32 * 0.11) + index as f32 * 0.37).sin();
                        particle.velocity += outward_dir * micro_pop_release * 0.22 * surface_release;
                        particle.velocity += tangent * 0.05 * surface_release * phase;
                        particle.ejecta_timer = (particle.ejecta_timer).max(0.08 + 0.12 * surface_release);
                    }
                }

                if particle.stored_energy > micro_pop_threshold {
                    let phase = ((self.step as usize + index * 17) % 19) == 0;
                    if phase {
                        let mut release_dir = if particle.last_force.norm_squared() > 0.0001 {
                            particle.last_force.normalize()
                        } else {
                            Vector2::new(1.0, 0.0)
                        };

                        if neighbor_count > 0.0 {
                            let local_centroid = neighbor_position_sum / neighbor_count;
                            let outward = particle.position - local_centroid;
                            if outward.norm_squared() > 0.0001 {
                                release_dir = outward.normalize();
                            }
                        }

                        let pop_strength = (particle.stored_energy - micro_pop_threshold + 0.08).min(1.25);
                        particle.velocity += release_dir * micro_pop_release * pop_strength;
                        particle.tension = (particle.tension - 0.16 * pop_strength).max(0.0);
                        particle.ejecta_timer = (0.12 + 0.25 * pop_strength).max(particle.ejecta_timer);
                        particle.stored_energy *= 0.45;
                    }
                }

                let speed = particle.velocity.norm();
                if speed > max_speed {
                    particle.velocity *= max_speed / speed;
                }

                particle.position += particle.velocity;
                confine_to_domain(&mut particle, self.config.world_size());

                particle.last_force = force;
                particle.base_tension = target_tension;
                particle.ejecta_timer = (particle.ejecta_timer - dt).max(0.0);
                particle.local_density = local_density;
                particle.shell_reactivity = shell_reactivity;
                particle.merge_factor = merge_factor;
                particle.influence_scale = influence_scale;
                particle.shell_overlap = shell_overlap;
                particle.spike_drive = spike_drive;
                particle
            })
            .collect();
    }

    fn apply_breach_events(&mut self, events: &[BreachEvent]) {
        for event in events {
            let cavity = &mut self.cavities[event.cavity_id];
            cavity.pop_count += 1;
            cavity.cooldown_steps = 120;

            for particle in &mut self.particles {
                if particle.fixed {
                    continue;
                }

                let offset = particle.position - cavity.center;
                let distance = offset.norm();
                if distance <= 0.0001 {
                    continue;
                }

                let radial = offset / distance;
                let alignment = radial.dot(&event.direction);
                let ring_distance = (distance - cavity.radius).abs();
                let shell_band = cavity.shell_thickness * 1.5;

                if alignment > 0.82 && distance > cavity.radius && ring_distance < shell_band {
                    particle.tension *= (1.0 - 0.35 * event.strength.min(1.0)).max(0.35);
                }

                if ring_distance < cavity.shell_thickness * 2.2 {
                    let tangent = Vector2::new(-radial.y, radial.x);
                    let wave = (alignment * PI).sin();
                    particle.velocity += tangent * wave * 0.08 * event.strength;
                }

                if alignment > 0.75 && distance < cavity.radius + cavity.shell_thickness * 1.2 {
                    particle.velocity -= radial * (0.18 + 0.15 * event.strength);
                }

                if alignment > 0.88
                    && distance >= cavity.radius + cavity.shell_thickness
                    && distance <= cavity.radius + cavity.shell_thickness * 3.8
                {
                    particle.velocity += radial * (0.24 + 0.25 * event.strength);
                    particle.tension = (particle.tension + 0.22 * event.strength).clamp(0.0, self.config.ambient_tension * 3.5);
                    particle.ejecta_timer = 0.75;
                }

                if distance <= cavity.radius + cavity.shell_thickness * 4.0 {
                    particle.velocity -= radial * 0.04 * event.strength;
                }
            }
        }
    }

    fn collect_metrics(&self) -> SimulationMetrics {
        let total_tension = self.particles.iter().map(|particle| particle.tension).sum();
        let true_voids = self.cavities.len();
        let mut shell_radii = BTreeMap::new();

        for cavity in &self.cavities {
            let mut total_radius = 0.0;
            let mut count = 0.0;

            for particle in &self.particles {
                if particle.kind == ParticleKind::TrueVoid {
                    continue;
                }

                let distance = (particle.position - cavity.center).norm();
                if distance >= cavity.radius && distance <= cavity.radius + cavity.shell_thickness * 2.5 {
                    total_radius += distance;
                    count += 1.0;
                }
            }

            let average = if count > 0.0 {
                total_radius / count
            } else {
                cavity.radius
            };

            shell_radii.insert(cavity.class, average);
        }

        SimulationMetrics {
            total_tension,
            true_voids,
            pops_this_frame: self.last_frame_pops,
            shell_radii,
        }
    }
}

fn seeded_cavities(config: &Config) -> Vec<Cavity> {
    let (width, height) = config.world_size();
    vec![
        Cavity {
            id: 0,
            center: Vector2::new(width * 0.28, height * 0.48),
            radius: config.up_core_radius,
            shell_thickness: config.up_shell_thickness,
            core_tension: 0.0,
            class: QuarkClass::Up,
            pop_count: 0,
            cooldown_steps: 0,
        },
        Cavity {
            id: 1,
            center: Vector2::new(width * 0.56, height * 0.42),
            radius: config.charm_core_radius,
            shell_thickness: config.charm_shell_thickness,
            core_tension: 0.0,
            class: QuarkClass::Charm,
            pop_count: 0,
            cooldown_steps: 0,
        },
        Cavity {
            id: 2,
            center: Vector2::new(width * 0.76, height * 0.62),
            radius: config.top_core_radius,
            shell_thickness: config.top_shell_thickness,
            core_tension: 0.0,
            class: QuarkClass::Top,
            pop_count: 0,
            cooldown_steps: 0,
        },
    ]
}

fn seeded_particles(config: &Config, cavities: &[Cavity]) -> Vec<Particle> {
    let (width, height) = config.world_size();
    let cell_area = (width * height) / config.particle_count as f32;
    let spacing = cell_area.sqrt().max(2.0);
    let mut particles = Vec::with_capacity(config.particle_count + 256);
    let mut rng = StdRng::seed_from_u64(7);

    let mut y = spacing * 0.5;
    while y < height && particles.len() < config.particle_count {
        let mut x = spacing * 0.5;
        while x < width && particles.len() < config.particle_count {
            let jitter = Vector2::new(
                rng.gen_range(-0.22 * spacing..0.22 * spacing),
                rng.gen_range(-0.22 * spacing..0.22 * spacing),
            );
            let position = Vector2::new(x, y) + jitter;

            let mut fixed = false;
            let mut kind = ParticleKind::Substrate;
            for cavity in cavities {
                if (position - cavity.center).norm() <= cavity.radius {
                    fixed = true;
                    kind = ParticleKind::TrueVoid;
                    break;
                }
            }

            let mut velocity = Vector2::new(
                rng.gen_range(-0.015..0.015),
                rng.gen_range(-0.015..0.015),
            );

            if fixed {
                velocity = Vector2::zeros();
            }

            let tension = if fixed {
                0.0
            } else {
                let (void_tension_boost, _, _) = void_boundary_response(
                    position,
                    cavities,
                    config.ambient_tension,
                    config.shell_gain,
                    config.void_geometry_gain,
                    config.void_shell_band_scale,
                );
                config.ambient_tension + void_tension_boost
            };

            particles.push(Particle {
                position,
                velocity,
                tension,
                kind,
                fixed,
                last_force: Vector2::zeros(),
                base_tension: tension,
                ejecta_timer: 0.0,
                local_density: 0.0,
                shell_reactivity: 0.0,
                merge_factor: 0.0,
                influence_scale: 1.0,
                stored_energy: 0.0,
                shell_overlap: 0.0,
                spike_drive: 0.0,
            });

            x += spacing;
        }
        y += spacing;
    }

    particles
}

fn void_boundary_response(
    position: Vector2<f32>,
    cavities: &[Cavity],
    ambient: f32,
    shell_gain: f32,
    void_geometry_gain: f32,
    shell_band_scale: f32,
) -> (f32, Vector2<f32>, f32) {
    let mut tension_boost = 0.0;
    let mut force = Vector2::zeros();
    let mut strongest_proximity: f32 = 0.0;

    for cavity in cavities {
        let offset = position - cavity.center;
        let distance = offset.norm();
        if distance <= cavity.radius {
            return (ambient, Vector2::zeros(), 1.0);
        }

        let ratio = (cavity.radius / distance.max(cavity.radius + 1.0e-4)).clamp(0.0, 1.0);
        let missing_fraction = ratio.asin() / PI;
        let band = cavity.shell_thickness.max(1.0) * shell_band_scale;
        let shell_distance = (distance - cavity.radius) / band;
        let shell_kernel = (-(shell_distance * shell_distance)).exp();
        let shell_strength = shell_gain * (ambient - cavity.core_tension) * missing_fraction * shell_kernel;

        tension_boost += shell_strength;
        strongest_proximity = strongest_proximity.max(missing_fraction * shell_kernel);

        if distance > 1.0e-4 {
            force += (offset / distance) * void_geometry_gain * missing_fraction * shell_kernel;
        }
    }

    (tension_boost, force, strongest_proximity)
}

fn expected_neighbor_count(config: &Config) -> f32 {
    let (width, height) = config.world_size();
    let particle_area = (width * height) / config.particle_count as f32;
    let influence_area = PI * config.interaction_radius * config.interaction_radius;
    influence_area / particle_area
}

pub fn substrate_spacing(config: &Config) -> f32 {
    let (width, height) = config.world_size();
    ((width * height) / config.particle_count as f32).sqrt()
}

/// Maximum distance a signal may travel in one update step.
///
/// This is the discrete causal horizon `ell_causal <= v_eff * dt` from the
/// simulation spec. It is floored at the substrate spacing so the nearest
/// neighbors stay reachable, and capped by the interaction radius used for
/// neighbor search. When `signal_speed * dt` is the binding term the horizon is
/// physical; when the cap or floor binds instead, the effective speed is set by
/// numerics rather than by `signal_speed`.
pub fn causal_horizon(config: &Config) -> f32 {
    let base_spacing = substrate_spacing(config).max(0.5);
    (config.signal_speed * config.dt)
        .min(config.interaction_radius)
        .max(base_spacing)
}

fn confine_to_domain(particle: &mut Particle, (width, height): (f32, f32)) {
    if particle.position.x < 0.0 {
        particle.position.x = 0.0;
        particle.velocity.x *= -0.4;
    } else if particle.position.x > width - 1.0 {
        particle.position.x = width - 1.0;
        particle.velocity.x *= -0.4;
    }

    if particle.position.y < 0.0 {
        particle.position.y = 0.0;
        particle.velocity.y *= -0.4;
    } else if particle.position.y > height - 1.0 {
        particle.position.y = height - 1.0;
        particle.velocity.y *= -0.4;
    }
}