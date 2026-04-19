pub const PARTICLE_UPDATE_WGSL: &str = r#"
struct Particle {
    position : vec2<f32>,
    velocity : vec2<f32>,
    tension : f32,
    kind : u32,
};

@group(0) @binding(0)
var<storage, read> particles_in : array<Particle>;

@group(0) @binding(1)
var<storage, read_write> particles_out : array<Particle>;

@compute @workgroup_size(128)
fn accumulate_forces(@builtin(global_invocation_id) gid : vec3<u32>) {
    let index = gid.x;
    let particle = particles_in[index];
    particles_out[index] = particle;
}
"#;

pub fn acceleration_notes() -> &'static str {
    let _shader_size = PARTICLE_UPDATE_WGSL.len();
    "CPU Rayon is the default path. The intended GPU path is a wgpu compute pass that consumes packed particle buffers, hashed cell spans, and cavity descriptors so force accumulation and tension relaxation can move to the GPU without changing the high-level model."
}