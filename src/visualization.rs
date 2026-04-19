use minifb::{Key, Window, WindowOptions};

use crate::physics::{ParticleKind, Simulation, SimulationMetrics};

pub struct Renderer {
    window: Window,
    buffer: Vec<u32>,
    width: usize,
    height: usize,
}

impl Renderer {
    pub fn new(width: usize, height: usize) -> Result<Self, minifb::Error> {
        let mut window = Window::new(
            "BedRock Tension-Field Prototype",
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
            buffer: vec![0x091018; width * height],
            width,
            height,
        })
    }

    pub fn is_open(&self) -> bool {
        self.window.is_open() && !self.window.is_key_down(Key::Escape)
    }

    pub fn draw(&mut self, simulation: &Simulation, _metrics: &SimulationMetrics) -> Result<(), minifb::Error> {
        self.buffer.fill(0x091018);

        for particle in &simulation.particles {
            let color = match particle.kind {
                ParticleKind::TrueVoid => 0x000000,
                _ if particle.ejecta_timer > 0.0 => 0xFFD27F,
                _ => tension_color(particle.tension, simulation.config.ambient_tension * 3.2),
            };

            let splat_radius = match particle.kind {
                ParticleKind::TrueVoid => 0,
                _ => {
                    let halo = ((particle.influence_scale - 1.0) * simulation.config.field_splat_radius as f32).round() as isize;
                    simulation.config.field_splat_radius as isize + halo
                }
            };

            draw_splat(
                &mut self.buffer,
                self.width,
                self.height,
                particle.position.x,
                particle.position.y,
                splat_radius,
                color,
            );

            if particle.kind != ParticleKind::TrueVoid && particle.merge_factor > 0.08 {
                let halo_radius = splat_radius + (particle.merge_factor * 2.0).round() as isize + 1;
                draw_splat(
                    &mut self.buffer,
                    self.width,
                    self.height,
                    particle.position.x,
                    particle.position.y,
                    halo_radius,
                    scale_color(color, 0.45),
                );
            }
        }

        for cavity in &simulation.cavities {
            draw_circle(
                &mut self.buffer,
                self.width,
                self.height,
                cavity.center.x as isize,
                cavity.center.y as isize,
                cavity.radius as isize,
                0x000000,
            );

            if simulation.config.show_debug_shells {
                draw_circle_outline(
                    &mut self.buffer,
                    self.width,
                    self.height,
                    cavity.center.x as isize,
                    cavity.center.y as isize,
                    (cavity.radius + cavity.shell_thickness) as isize,
                    0xE8EEF2,
                );
            }
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

fn draw_circle(buffer: &mut [u32], width: usize, height: usize, cx: isize, cy: isize, radius: isize, color: u32) {
    let radius_sq = radius * radius;
    for y in (cy - radius).max(0)..=(cy + radius).min(height as isize - 1) {
        for x in (cx - radius).max(0)..=(cx + radius).min(width as isize - 1) {
            let dx = x - cx;
            let dy = y - cy;
            if dx * dx + dy * dy <= radius_sq {
                buffer[y as usize * width + x as usize] = color;
            }
        }
    }
}

fn draw_circle_outline(buffer: &mut [u32], width: usize, height: usize, cx: isize, cy: isize, radius: isize, color: u32) {
    let radius_sq = radius * radius;
    let inner_sq = (radius - 1).max(0) * (radius - 1).max(0);
    for y in (cy - radius).max(0)..=(cy + radius).min(height as isize - 1) {
        for x in (cx - radius).max(0)..=(cx + radius).min(width as isize - 1) {
            let dx = x - cx;
            let dy = y - cy;
            let distance_sq = dx * dx + dy * dy;
            if distance_sq <= radius_sq && distance_sq >= inner_sq {
                buffer[y as usize * width + x as usize] = color;
            }
        }
    }
}

fn draw_splat(
    buffer: &mut [u32],
    width: usize,
    height: usize,
    cx: f32,
    cy: f32,
    radius: isize,
    color: u32,
) {
    let cx = cx as isize;
    let cy = cy as isize;

    if radius <= 0 {
        if cx >= 0 && cy >= 0 && cx < width as isize && cy < height as isize {
            blend_pixel(buffer, width, cx as usize, cy as usize, color, 1.0);
        }
        return;
    }

    let radius_sq = (radius * radius) as f32;
    for y in (cy - radius).max(0)..=(cy + radius).min(height as isize - 1) {
        for x in (cx - radius).max(0)..=(cx + radius).min(width as isize - 1) {
            let dx = (x - cx) as f32;
            let dy = (y - cy) as f32;
            let distance_sq = dx * dx + dy * dy;
            if distance_sq > radius_sq {
                continue;
            }

            let weight = (1.0 - distance_sq / radius_sq.max(1.0)).clamp(0.15, 1.0);
            blend_pixel(buffer, width, x as usize, y as usize, color, weight);
        }
    }
}

fn blend_pixel(buffer: &mut [u32], width: usize, x: usize, y: usize, color: u32, weight: f32) {
    let index = y * width + x;
    let existing = buffer[index];

    let er = ((existing >> 16) & 0xFF) as f32;
    let eg = ((existing >> 8) & 0xFF) as f32;
    let eb = (existing & 0xFF) as f32;

    let nr = ((color >> 16) & 0xFF) as f32;
    let ng = ((color >> 8) & 0xFF) as f32;
    let nb = (color & 0xFF) as f32;

    let blended_r = er.max(nr * weight).min(255.0) as u32;
    let blended_g = eg.max(ng * weight).min(255.0) as u32;
    let blended_b = eb.max(nb * weight).min(255.0) as u32;

    buffer[index] = (blended_r << 16) | (blended_g << 8) | blended_b;
}

fn scale_color(color: u32, scale: f32) -> u32 {
    let red = (((color >> 16) & 0xFF) as f32 * scale).clamp(0.0, 255.0) as u32;
    let green = (((color >> 8) & 0xFF) as f32 * scale).clamp(0.0, 255.0) as u32;
    let blue = ((color & 0xFF) as f32 * scale).clamp(0.0, 255.0) as u32;
    (red << 16) | (green << 8) | blue
}