use std::collections::BTreeMap;
use std::f32::consts::PI;

use nalgebra::Vector3;
use rand::{rngs::StdRng, Rng, SeedableRng};
use rayon::prelude::*;

use crate::{config::Config, neighbor_search_3d::SpatialHashGrid3d};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParticleKind3d {
    Substrate,
    TrueVoid,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum QuarkClass3d {
    Up,
    Charm,
    Top,
}

impl QuarkClass3d {
    pub fn label(self) -> &'static str {
        match self {
            Self::Up => "up",
            Self::Charm => "charm",
            Self::Top => "top",
        }
    }
}

#[derive(Debug, Clone)]
pub struct VoidSphere {
    pub center: Vector3<f32>,
    pub radius: f32,
    pub shell_thickness: f32,
    pub class: QuarkClass3d,
}

#[derive(Debug, Clone)]
pub struct Particle3d {
    pub position: Vector3<f32>,
    pub velocity: Vector3<f32>,
    pub tension: f32,
    pub kind: ParticleKind3d,
    pub fixed: bool,
    pub influence_scale: f32,
    pub local_density: f32,
    pub shell_reactivity: f32,
    pub merge_factor: f32,
    pub stored_energy: f32,
    pub ejecta_timer: f32,
    pub shell_overlap: f32,
    pub spike_drive: f32,
}

#[derive(Debug, Clone)]
pub struct SimulationMetrics3d {
    pub total_tension: f32,
    pub true_voids: usize,
    pub micro_pops: usize,
    pub shell_radii: BTreeMap<QuarkClass3d, f32>,
}

pub struct Simulation3d {
    pub config: Config,
    pub particles: Vec<Particle3d>,
    pub voids: Vec<VoidSphere>,
    pub grid: SpatialHashGrid3d,
    pub step: u64,
    last_frame_micro_pops: usize,
}

impl Simulation3d {
    pub fn new(config: Config) -> Self {
        let voids = seeded_voids(&config);
        let particles = seeded_particles(&config, &voids);
        let mut grid = SpatialHashGrid3d::new(config.cell_size);
        grid.rebuild(&particles);

        Self {
            config,
            particles,
            voids,
            grid,
            step: 0,
            last_frame_micro_pops: 0,
        }
    }

    pub fn update(&mut self) -> SimulationMetrics3d {
        self.last_frame_micro_pops = 0;
        self.grid.rebuild(&self.particles);
        self.apply_force_step();
        self.step += 1;
        self.collect_metrics()
    }

    fn apply_force_step(&mut self) {
        let interaction_radius_sq = self.config.interaction_radius * self.config.interaction_radius;
        let dt = self.config.dt;
        let ambient = self.config.ambient_tension;
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
        let signal_speed = self.config.signal_speed;
        let tension_relaxation = self.config.tension_relaxation;
        let max_speed = self.config.max_speed;
        let velocity_drag = 1.0 - self.config.velocity_damping * dt;
        let particles = self.particles.clone();
        let voids = self.voids.clone();
        let grid = &self.grid;
        let expected_neighbors = expected_neighbor_count_3d(&self.config).max(1.0);
        let base_spacing = substrate_spacing_3d(&self.config).max(0.5);
        let causal_horizon = (signal_speed * dt).min(self.config.interaction_radius).max(base_spacing);
        let min_shell_gap = base_spacing * 0.35;
        let step_seed = self.step as usize;

        let updated: Vec<(Particle3d, bool)> = (0..particles.len())
            .into_par_iter()
            .map(|index| {
                let mut particle = particles[index].clone();

                if particle.fixed {
                    particle.tension = 0.0;
                    particle.velocity = Vector3::zeros();
                    particle.local_density = 0.0;
                    particle.shell_reactivity = 0.0;
                    particle.merge_factor = 0.0;
                    particle.influence_scale = 1.0;
                    particle.stored_energy = 0.0;
                    particle.shell_overlap = 0.0;
                    particle.spike_drive = 0.0;
                    return (particle, false);
                }

                let mut neighbor_indices = Vec::with_capacity(48);
                grid.query_neighbors(particle.position, &mut neighbor_indices);

                let mut force = Vector3::zeros();
                let mut neighbor_tension_sum = 0.0;
                let mut neighbor_count = 0.0;
                let mut neighbor_velocity_sum = Vector3::zeros();
                let mut neighbor_position_sum = Vector3::zeros();
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

                let (void_tension_boost, void_force, shell_proximity) = void_boundary_response_3d(
                    particle.position,
                    &voids,
                    ambient,
                    self.config.shell_gain,
                    void_geometry_gain,
                    self.config.void_shell_band_scale,
                );
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
                    .clamp(0.85, 3.4);

                force += void_force * (1.0 + 0.35 * crowding + 0.25 * shell_reactivity);

                if neighbor_count > 0.0 {
                    let local_velocity = neighbor_velocity_sum / neighbor_count;
                    let local_centroid = neighbor_position_sum / neighbor_count;
                    let collapse_offset = local_centroid - particle.position;
                    force += collapse_offset * collapse_gain * influence_scale * (0.18 + shell_proximity) * crowding.max(0.05);
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
                    + 0.16 * merge_factor * (0.35 + shell_proximity)
                    - surface_release
                    + particle.ejecta_timer.max(0.0) * 0.05;

                let max_target_step = attraction_ramp_gain * (0.2 + available_supply) * dt * 10.0;
                let target_tension = particle.tension + (raw_target - particle.tension).clamp(-max_target_step, max_target_step);

                particle.tension += (target_tension - particle.tension) * tension_relaxation * dt * 10.0;
                particle.tension = particle.tension.clamp(0.0, ambient * 3.4);

                particle.velocity += force * dt;
                particle.velocity *= velocity_drag.max(0.0);

                let mut micro_pop = false;
                let stored_energy = (particle.stored_energy
                    + (merge_factor * (0.28 + shell_proximity) + 0.16 * crowding + 0.4 * surface_release) * dt
                    - micro_pop_release * 0.18 * dt)
                    .max(0.0);
                particle.stored_energy = stored_energy;

                if surface_release > 0.02 && neighbor_count > 0.0 {
                    let local_centroid = neighbor_position_sum / neighbor_count;
                    let outward = particle.position - local_centroid;
                    if outward.norm_squared() > 0.0001 {
                        let outward_dir = outward.normalize();
                        let axis = if outward_dir.z.abs() < 0.9 {
                            Vector3::new(0.0, 0.0, 1.0)
                        } else {
                            Vector3::new(1.0, 0.0, 0.0)
                        };
                        let tangent = outward_dir.cross(&axis).normalize();
                        let phase = ((step_seed as f32 * 0.09) + index as f32 * 0.31).sin();
                        particle.velocity += outward_dir * micro_pop_release * 0.22 * surface_release;
                        particle.velocity += tangent * 0.05 * surface_release * phase;
                        particle.ejecta_timer = particle.ejecta_timer.max(0.08 + 0.12 * surface_release);
                    }
                }

                if particle.stored_energy > micro_pop_threshold {
                    let phase = ((step_seed + index * 23) % 29) == 0;
                    if phase {
                        let mut release_dir = if neighbor_count > 0.0 {
                            let local_centroid = neighbor_position_sum / neighbor_count;
                            let outward = particle.position - local_centroid;
                            if outward.norm_squared() > 0.0001 {
                                outward.normalize()
                            } else {
                                Vector3::new(1.0, 0.0, 0.0)
                            }
                        } else {
                            Vector3::new(1.0, 0.0, 0.0)
                        };

                        if release_dir.norm_squared() <= 0.0001 {
                            release_dir = Vector3::new(1.0, 0.0, 0.0);
                        }

                        let pop_strength = (particle.stored_energy - micro_pop_threshold + 0.06).min(1.2);
                        particle.velocity += release_dir * micro_pop_release * pop_strength;
                        particle.tension = (particle.tension - 0.14 * pop_strength).max(0.0);
                        particle.ejecta_timer = (0.1 + 0.22 * pop_strength).max(particle.ejecta_timer);
                        particle.stored_energy *= 0.48;
                        micro_pop = true;
                    }
                }

                let speed = particle.velocity.norm();
                if speed > max_speed {
                    particle.velocity *= max_speed / speed;
                }

                particle.position += particle.velocity;
                confine_to_domain_3d(&mut particle, self.config.world_size_3d());
                particle.ejecta_timer = (particle.ejecta_timer - dt).max(0.0);
                particle.local_density = local_density;
                particle.shell_reactivity = shell_reactivity;
                particle.merge_factor = merge_factor;
                particle.influence_scale = influence_scale;
                particle.shell_overlap = shell_overlap;
                particle.spike_drive = spike_drive;

                (particle, micro_pop)
            })
            .collect();

        self.last_frame_micro_pops = updated.iter().filter(|(_, micro_pop)| *micro_pop).count();
        self.particles = updated.into_iter().map(|(particle, _)| particle).collect();
    }

    fn collect_metrics(&self) -> SimulationMetrics3d {
        let total_tension = self.particles.iter().map(|particle| particle.tension).sum();
        let true_voids = self.voids.len();
        let mut shell_radii = BTreeMap::new();

        for void in &self.voids {
            let mut total_radius = 0.0;
            let mut count = 0.0;

            for particle in &self.particles {
                if particle.kind == ParticleKind3d::TrueVoid {
                    continue;
                }

                let distance = (particle.position - void.center).norm();
                if distance >= void.radius && distance <= void.radius + void.shell_thickness * 2.5 {
                    total_radius += distance;
                    count += 1.0;
                }
            }

            let average = if count > 0.0 { total_radius / count } else { void.radius };
            shell_radii.insert(void.class, average);
        }

        SimulationMetrics3d {
            total_tension,
            true_voids,
            micro_pops: self.last_frame_micro_pops,
            shell_radii,
        }
    }
}

fn seeded_voids(config: &Config) -> Vec<VoidSphere> {
    let (width, height, depth) = config.world_size_3d();
    vec![
        VoidSphere {
            center: Vector3::new(width * 0.28, height * 0.48, depth * 0.42),
            radius: config.up_core_radius,
            shell_thickness: config.up_shell_thickness,
            class: QuarkClass3d::Up,
        },
        VoidSphere {
            center: Vector3::new(width * 0.56, height * 0.42, depth * 0.56),
            radius: config.charm_core_radius,
            shell_thickness: config.charm_shell_thickness,
            class: QuarkClass3d::Charm,
        },
        VoidSphere {
            center: Vector3::new(width * 0.76, height * 0.62, depth * 0.50),
            radius: config.top_core_radius,
            shell_thickness: config.top_shell_thickness,
            class: QuarkClass3d::Top,
        },
    ]
}

fn seeded_particles(config: &Config, voids: &[VoidSphere]) -> Vec<Particle3d> {
    let (width, height, depth) = config.world_size_3d();
    let total_count = config.active_particle_count();
    let cell_volume = (width * height * depth) / total_count as f32;
    let spacing = cell_volume.cbrt().max(1.6);
    let mut particles = Vec::with_capacity(total_count + 512);
    let mut rng = StdRng::seed_from_u64(11);

    let mut z = spacing * 0.5;
    while z < depth && particles.len() < total_count {
        let mut y = spacing * 0.5;
        while y < height && particles.len() < total_count {
            let mut x = spacing * 0.5;
            while x < width && particles.len() < total_count {
                let jitter = Vector3::new(
                    rng.gen_range(-0.18 * spacing..0.18 * spacing),
                    rng.gen_range(-0.18 * spacing..0.18 * spacing),
                    rng.gen_range(-0.18 * spacing..0.18 * spacing),
                );
                let position = Vector3::new(x, y, z) + jitter;

                let mut fixed = false;
                for void in voids {
                    if (position - void.center).norm() <= void.radius {
                        fixed = true;
                        break;
                    }
                }

                let velocity = if fixed {
                    Vector3::zeros()
                } else {
                    Vector3::new(
                        rng.gen_range(-0.01..0.01),
                        rng.gen_range(-0.01..0.01),
                        rng.gen_range(-0.01..0.01),
                    )
                };

                let tension = if fixed {
                    0.0
                } else {
                    let (void_tension_boost, _, _) = void_boundary_response_3d(
                        position,
                        voids,
                        config.ambient_tension,
                        config.shell_gain,
                        config.void_geometry_gain,
                        config.void_shell_band_scale,
                    );
                    config.ambient_tension + void_tension_boost
                };

                particles.push(Particle3d {
                    position,
                    velocity,
                    tension,
                    kind: if fixed { ParticleKind3d::TrueVoid } else { ParticleKind3d::Substrate },
                    fixed,
                    influence_scale: 1.0,
                    local_density: 0.0,
                    shell_reactivity: 0.0,
                    merge_factor: 0.0,
                    stored_energy: 0.0,
                    ejecta_timer: 0.0,
                    shell_overlap: 0.0,
                    spike_drive: 0.0,
                });

                x += spacing;
            }
            y += spacing;
        }
        z += spacing;
    }

    particles
}

fn void_boundary_response_3d(
    position: Vector3<f32>,
    voids: &[VoidSphere],
    ambient: f32,
    shell_gain: f32,
    void_geometry_gain: f32,
    shell_band_scale: f32,
) -> (f32, Vector3<f32>, f32) {
    let mut tension_boost = 0.0;
    let mut force = Vector3::zeros();
    let mut strongest_proximity: f32 = 0.0;

    for void in voids {
        let offset = position - void.center;
        let distance = offset.norm();
        if distance <= void.radius {
            return (ambient, Vector3::zeros(), 1.0);
        }

        let ratio = (void.radius / distance.max(void.radius + 1.0e-4)).clamp(0.0, 1.0);
        let cos_theta = (1.0 - ratio * ratio).sqrt();
        let missing_fraction = 0.5 * (1.0 - cos_theta);
        let band = void.shell_thickness.max(1.0) * shell_band_scale;
        let shell_distance = (distance - void.radius) / band;
        let shell_kernel = (-(shell_distance * shell_distance)).exp();
        let shell_strength = shell_gain * ambient * missing_fraction * shell_kernel;

        tension_boost += shell_strength;
        strongest_proximity = strongest_proximity.max(missing_fraction * shell_kernel);

        if distance > 1.0e-4 {
            force += (offset / distance) * void_geometry_gain * missing_fraction * shell_kernel;
        }
    }

    (tension_boost, force, strongest_proximity)
}

fn expected_neighbor_count_3d(config: &Config) -> f32 {
    let (width, height, depth) = config.world_size_3d();
    let particle_volume = (width * height * depth) / config.active_particle_count() as f32;
    let influence_volume = (4.0 / 3.0) * PI * config.interaction_radius.powi(3);
    influence_volume / particle_volume
}

fn substrate_spacing_3d(config: &Config) -> f32 {
    let (width, height, depth) = config.world_size_3d();
    ((width * height * depth) / config.active_particle_count() as f32).cbrt()
}

fn confine_to_domain_3d(particle: &mut Particle3d, (width, height, depth): (f32, f32, f32)) {
    if particle.position.x < 0.0 {
        particle.position.x = 0.0;
        particle.velocity.x *= -0.35;
    } else if particle.position.x > width - 1.0 {
        particle.position.x = width - 1.0;
        particle.velocity.x *= -0.35;
    }

    if particle.position.y < 0.0 {
        particle.position.y = 0.0;
        particle.velocity.y *= -0.35;
    } else if particle.position.y > height - 1.0 {
        particle.position.y = height - 1.0;
        particle.velocity.y *= -0.35;
    }

    if particle.position.z < 0.0 {
        particle.position.z = 0.0;
        particle.velocity.z *= -0.35;
    } else if particle.position.z > depth - 1.0 {
        particle.position.z = depth - 1.0;
        particle.velocity.z *= -0.35;
    }
}