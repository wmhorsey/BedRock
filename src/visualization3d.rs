use minifb::{Key, Window, WindowOptions};
use nalgebra::{Rotation3, Vector3};

use crate::physics3d::{ParticleKind3d, Simulation3d, SimulationMetrics3d};

pub struct Renderer3d {
    window: Window,
    buffer: Vec<u32>,
    depth_buffer: Vec<f32>,
    width: usize,
    height: usize,
    yaw: f32,
}

impl Renderer3d {
    pub fn new(width: usize, height: usize) -> Result<Self, minifb::Error> {
        let mut window = Window::new(
            "BedRock Tension-Field Prototype 3D",
            width,
            height,
            WindowOptions {
                resize: false,
                scale: minifb::Scale::X1,
                ..WindowOptions::default()
            },
        )?;
        window.set_target_fps(60);

        Ok(Self {
            window,
            buffer: vec![0x050B13; width * height],
            depth_buffer: vec![f32::NEG_INFINITY; width * height],
            width,
            height,
            yaw: 0.0,
        })
    }

    pub fn is_open(&self) -> bool {
        self.window.is_open() && !self.window.is_key_down(Key::Escape)
    }

    pub fn draw(&mut self, simulation: &Simulation3d, _metrics: &SimulationMetrics3d) -> Result<(), minifb::Error> {
        self.buffer.fill(0x050B13);
        self.depth_buffer.fill(f32::NEG_INFINITY);
        self.yaw += simulation.config.camera_yaw_speed;

        let rotation = Rotation3::from_euler_angles(simulation.config.camera_pitch, self.yaw, 0.0);
        let (world_w, world_h, world_d) = simulation.config.world_size_3d();
        let center = Vector3::new(world_w * 0.5, world_h * 0.5, world_d * 0.5);
        let projection_scale = simulation.config.projection_scale;
        let far_stride = simulation.config.far_field_stride_3d.max(1);
        let far_fade_cutoff = simulation.config.far_field_fade_cutoff;

        for (index, particle) in simulation.particles.iter().enumerate() {
            let local = particle.position - center;
            let rotated = rotation * local;
            let depth = rotated.z + simulation.config.camera_distance;
            if depth <= 1.0 {
                continue;
            }

            let screen_x = self.width as f32 * 0.5 + rotated.x * projection_scale / depth;
            let screen_y = self.height as f32 * 0.5 + rotated.y * projection_scale / depth;
            if screen_x < -8.0 || screen_y < -8.0 || screen_x >= self.width as f32 + 8.0 || screen_y >= self.height as f32 + 8.0 {
                continue;
            }

            let base_color = match particle.kind {
                ParticleKind3d::TrueVoid => 0x000000,
                _ if particle.ejecta_timer > 0.0 => 0xFFD27F,
                _ => tension_color(particle.tension, simulation.config.ambient_tension * 3.4),
            };

            let depth_fade = (1.0 - depth / (simulation.config.camera_distance + world_d * 1.4)).clamp(0.2, 1.0);
            if particle.kind != ParticleKind3d::TrueVoid
                && depth_fade <= far_fade_cutoff
                && particle.merge_factor < 0.05
                && particle.shell_reactivity < 0.05
                && particle.ejecta_timer <= 0.0
                && index % far_stride != 0
            {
                continue;
            }

            let radius = match particle.kind {
                ParticleKind3d::TrueVoid => 1,
                _ => {
                    let base = simulation.config.field_splat_radius as f32;
                    let halo = (particle.influence_scale - 1.0) * base;
                    ((base + halo) * (simulation.config.camera_distance / depth).sqrt().clamp(0.55, 2.4)).round() as isize
                }
            };

            draw_depth_splat(
                &mut self.buffer,
                &mut self.depth_buffer,
                self.width,
                self.height,
                screen_x,
                screen_y,
                radius.max(1),
                scale_color(base_color, depth_fade),
                depth,
            );
        }

        for void in &simulation.voids {
            let local = void.center - center;
            let rotated = rotation * local;
            let depth = rotated.z + simulation.config.camera_distance;
            if depth <= 1.0 {
                continue;
            }

            let screen_x = self.width as f32 * 0.5 + rotated.x * projection_scale / depth;
            let screen_y = self.height as f32 * 0.5 + rotated.y * projection_scale / depth;
            let radius = (void.radius * projection_scale / depth).round() as isize;
            draw_depth_splat(
                &mut self.buffer,
                &mut self.depth_buffer,
                self.width,
                self.height,
                screen_x,
                screen_y,
                radius.max(2),
                0x000000,
                depth + 10.0,
            );
        }

        self.window.update_with_buffer(&self.buffer, self.width, self.height)
    }
}

fn tension_color(value: f32, max_value: f32) -> u32 {
    let normalized = (value / max_value).clamp(0.0, 1.0);
    let red = (normalized * 255.0) as u32;
    let green = ((1.0 - (normalized - 0.5).abs() * 2.0).clamp(0.0, 1.0) * 144.0) as u32;
    let blue = ((1.0 - normalized) * 255.0) as u32;
    (red << 16) | (green << 8) | blue
}

fn draw_depth_splat(
    buffer: &mut [u32],
    depth_buffer: &mut [f32],
    width: usize,
    height: usize,
    cx: f32,
    cy: f32,
    radius: isize,
    color: u32,
    depth: f32,
) {
    let cx = cx as isize;
    let cy = cy as isize;
    let radius_sq = (radius * radius) as f32;

    for y in (cy - radius).max(0)..=(cy + radius).min(height as isize - 1) {
        for x in (cx - radius).max(0)..=(cx + radius).min(width as isize - 1) {
            let dx = (x - cx) as f32;
            let dy = (y - cy) as f32;
            let distance_sq = dx * dx + dy * dy;
            if distance_sq > radius_sq {
                continue;
            }

            let index = y as usize * width + x as usize;
            let fragment_depth = depth - 0.02 * distance_sq;
            if fragment_depth < depth_buffer[index] {
                continue;
            }

            let weight = (1.0 - distance_sq / radius_sq.max(1.0)).clamp(0.1, 1.0);
            let existing = buffer[index];
            buffer[index] = blend_pixel(existing, color, weight);
            depth_buffer[index] = fragment_depth;
        }
    }
}

fn blend_pixel(existing: u32, color: u32, weight: f32) -> u32 {
    let er = ((existing >> 16) & 0xFF) as f32;
    let eg = ((existing >> 8) & 0xFF) as f32;
    let eb = (existing & 0xFF) as f32;

    let nr = ((color >> 16) & 0xFF) as f32;
    let ng = ((color >> 8) & 0xFF) as f32;
    let nb = (color & 0xFF) as f32;

    let r = er.max(nr * weight).min(255.0) as u32;
    let g = eg.max(ng * weight).min(255.0) as u32;
    let b = eb.max(nb * weight).min(255.0) as u32;
    (r << 16) | (g << 8) | b
}

fn scale_color(color: u32, scale: f32) -> u32 {
    let red = (((color >> 16) & 0xFF) as f32 * scale).clamp(0.0, 255.0) as u32;
    let green = (((color >> 8) & 0xFF) as f32 * scale).clamp(0.0, 255.0) as u32;
    let blue = ((color & 0xFF) as f32 * scale).clamp(0.0, 255.0) as u32;
    (red << 16) | (green << 8) | blue
}