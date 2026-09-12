use clap::{Parser, ValueEnum};

#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum SimulationMode {
    #[value(name = "2d")]
    TwoD,
    #[value(name = "3d")]
    ThreeD,
}

#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum ConeModel {
    /// Current relaxational particle engine (diffusive; no wave).
    Relax,
    /// Conservative acoustic field closure (real propagating wave).
    Wave,
}

#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum ConeSource {
    /// Compact disk pulse at the domain center.
    Disk,
    /// Shell pop: a tension ring released with zero flux (physical collapse).
    Shell,
}

#[derive(Parser, Debug, Clone)]
#[command(author, version, about = "BedRock single-substance tension-field prototype")]
pub struct Config {
    #[arg(long, value_enum, default_value_t = SimulationMode::ThreeD)]
    pub mode: SimulationMode,

    #[arg(long, default_value_t = 960)]
    pub width: usize,

    #[arg(long, default_value_t = 720)]
    pub height: usize,

    #[arg(long, default_value_t = 720.0)]
    pub depth: f32,

    #[arg(long, default_value_t = 48000)]
    pub particle_count: usize,

    #[arg(long, default_value_t = 120000)]
    pub particle_count_3d: usize,

    #[arg(long, default_value_t = 1.0)]
    pub ambient_tension: f32,

    #[arg(long, default_value_t = 0.72)]
    pub attraction_strength: f32,

    #[arg(long, default_value_t = 0.38)]
    pub density_reactivity: f32,

    #[arg(long, default_value_t = 0.42)]
    pub density_force_gain: f32,

    #[arg(long, default_value_t = 1.12)]
    pub merge_density_threshold: f32,

    #[arg(long, default_value_t = 0.75)]
    pub merge_influence_gain: f32,

    #[arg(long, default_value_t = 0.55)]
    pub attraction_ramp_gain: f32,

    #[arg(long, default_value_t = 0.42)]
    pub micro_pop_threshold: f32,

    #[arg(long, default_value_t = 0.34)]
    pub micro_pop_release: f32,

    #[arg(long, default_value_t = 0.28)]
    pub flow_coupling: f32,

    #[arg(long, default_value_t = 0.32)]
    pub collapse_gain: f32,

    #[arg(long, default_value_t = 1.12)]
    pub breach_threshold_base: f32,

    #[arg(long, default_value_t = 2.2)]
    pub curvature_scale: f32,

    #[arg(long, default_value_t = 16.0)]
    pub interaction_radius: f32,

    #[arg(long, default_value_t = 220.0)]
    pub signal_speed: f32,

    #[arg(long, default_value_t = 0.68)]
    pub shell_source_scale: f32,

    #[arg(long, default_value_t = 1.4)]
    pub shell_overlap_gain: f32,

    #[arg(long, default_value_t = 0.28)]
    pub spike_threshold: f32,

    #[arg(long, default_value_t = 0.52)]
    pub shell_reform_gain: f32,

    #[arg(long, default_value_t = 0.36)]
    pub void_geometry_gain: f32,

    #[arg(long, default_value_t = 1.15)]
    pub void_shell_band_scale: f32,

    #[arg(long, default_value_t = 16.0)]
    pub cell_size: f32,

    #[arg(long, default_value_t = 0.08)]
    pub dt: f32,

    #[arg(long, default_value_t = 2)]
    pub substeps_per_frame: usize,

    #[arg(long, default_value_t = 0)]
    pub headless_steps: usize,

    #[arg(long, default_value_t = false)]
    pub show_debug_shells: bool,

    #[arg(long, default_value_t = 2)]
    pub field_splat_radius: usize,

    #[arg(long, default_value_t = 0.82)]
    pub shell_gain: f32,

    #[arg(long, default_value_t = 28.0)]
    pub up_core_radius: f32,

    #[arg(long, default_value_t = 18.0)]
    pub charm_core_radius: f32,

    #[arg(long, default_value_t = 10.0)]
    pub top_core_radius: f32,

    #[arg(long, default_value_t = 14.0)]
    pub up_shell_thickness: f32,

    #[arg(long, default_value_t = 10.0)]
    pub charm_shell_thickness: f32,

    #[arg(long, default_value_t = 6.0)]
    pub top_shell_thickness: f32,

    #[arg(long, default_value_t = 0.18)]
    pub velocity_damping: f32,

    #[arg(long, default_value_t = 0.22)]
    pub tension_relaxation: f32,

    #[arg(long, default_value_t = 0.55)]
    pub max_speed: f32,

    #[arg(long, default_value_t = 920.0)]
    pub camera_distance: f32,

    #[arg(long, default_value_t = 640.0)]
    pub projection_scale: f32,

    #[arg(long, default_value_t = 0.0035)]
    pub camera_yaw_speed: f32,

    #[arg(long, default_value_t = 0.42)]
    pub camera_pitch: f32,

    #[arg(long, default_value_t = 2)]
    pub far_field_stride_3d: usize,

    #[arg(long, default_value_t = 0.45)]
    pub far_field_fade_cutoff: f32,

    /// Run the causal-front cone test instead of the normal simulation.
    #[arg(long, default_value_t = false)]
    pub cone_test: bool,

    /// Constitutive law used by the cone test.
    #[arg(long, value_enum, default_value_t = ConeModel::Wave)]
    pub cone_model: ConeModel,

    /// Excitation used by the wave cone test.
    #[arg(long, value_enum, default_value_t = ConeSource::Disk)]
    pub cone_source: ConeSource,

    /// Ring radius for the shell-pop source (wave model).
    #[arg(long, default_value_t = 120.0)]
    pub cone_shell_radius: f32,

    /// Ring half-width for the shell-pop source (wave model).
    #[arg(long, default_value_t = 18.0)]
    pub cone_shell_width: f32,

    /// Number of update steps to evolve the cone pulse.
    #[arg(long, default_value_t = 80)]
    pub cone_steps: usize,

    /// Radius of the localized tension pulse injected at the domain center.
    #[arg(long, default_value_t = 12.0)]
    pub cone_perturb_radius: f32,

    /// Peak tension added above ambient at the center of the pulse.
    #[arg(long, default_value_t = 1.5)]
    pub cone_perturb_amplitude: f32,

    /// Angular bins used to measure front radius and isotropy.
    #[arg(long, default_value_t = 24)]
    pub cone_bins: usize,

    /// CSV path for the per-step cone measurement (empty disables the file).
    #[arg(long, default_value = "cone_test.csv")]
    pub cone_csv: String,
}

impl Config {
    pub fn apply_runtime_defaults(&mut self) {
        if self.mode == SimulationMode::ThreeD {
            if self.particle_count_3d == 120000 {
                self.particle_count_3d = 140000;
            }
            if (self.interaction_radius - 16.0).abs() < f32::EPSILON {
                self.interaction_radius = 12.0;
            }
            if (self.camera_yaw_speed - 0.0035).abs() < f32::EPSILON {
                self.camera_yaw_speed = 0.004;
            }
        }
    }

    /// Cone-friendly setup: force the 2D reference model and, only where the
    /// user left global defaults, pick a domain and parameters that make the
    /// selected cone model measurable out of the box.
    pub fn apply_cone_defaults(&mut self) {
        self.mode = SimulationMode::TwoD;
        match self.cone_model {
            // Relaxational engine: keep signal_speed * dt between the substrate
            // spacing and the interaction radius so the causal horizon is
            // physical rather than clamped by the neighbor cutoff.
            ConeModel::Relax => {
                if self.width == 960 {
                    self.width = 1400;
                }
                if self.height == 720 {
                    self.height = 1400;
                }
                if self.particle_count == 48000 {
                    self.particle_count = 120000;
                }
                if (self.signal_speed - 220.0).abs() < f32::EPSILON {
                    self.signal_speed = 110.0;
                }
            }
            // Conservative wave: signal_speed is now the wave speed c. Pick a
            // sample spacing and support so the leapfrog sits near CFL ~0.3
            // (c*dt ~ 0.3 * spacing) with a well-resolved kernel (~4 spacings).
            ConeModel::Wave => {
                if self.width == 960 {
                    self.width = 1200;
                }
                if self.height == 720 {
                    self.height = 1200;
                }
                if self.particle_count == 48000 {
                    self.particle_count = 40000;
                }
                if (self.interaction_radius - 16.0).abs() < f32::EPSILON {
                    self.interaction_radius = 24.0;
                }
                if (self.cell_size - 16.0).abs() < f32::EPSILON {
                    self.cell_size = 24.0;
                }
                if (self.signal_speed - 220.0).abs() < f32::EPSILON {
                    self.signal_speed = 22.5;
                }
                if self.cone_steps == 80 {
                    self.cone_steps = 300;
                }
                if (self.cone_perturb_radius - 12.0).abs() < f32::EPSILON {
                    self.cone_perturb_radius = 18.0;
                }
            }
        }
    }

    pub fn world_size(&self) -> (f32, f32) {
        (self.width as f32, self.height as f32)
    }

    pub fn world_size_3d(&self) -> (f32, f32, f32) {
        (self.width as f32, self.height as f32, self.depth)
    }

    pub fn active_particle_count(&self) -> usize {
        match self.mode {
            SimulationMode::TwoD => self.particle_count,
            SimulationMode::ThreeD => self.particle_count_3d,
        }
    }
}