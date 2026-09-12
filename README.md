# BedRock Tension-Field Prototype

This crate contains a 2D Rust prototype for a single-substance tension-field simulation driven by mesh-free particles. It starts CPU-first with Rayon parallelism, uses a spatial hash grid for neighbor lookups, and keeps a documented `wgpu` compute path ready for later offload.

## How To Run

```powershell
cargo run --release
```

The default run now opens the 3D view. To force the old 2D mode:

```powershell
cargo run --release -- --mode 2d
```

Short Cargo aliases:

```powershell
cargo sim
cargo sim2d
cargo headless3d
```

Run the higher-density 3D view:

```powershell
cargo run --release -- --mode 3d --particle-count-3d 120000 --field-splat-radius 2 --interaction-radius 12 --depth 720
```

Plain `cargo run --release` now auto-applies the tuned 3D defaults when you stay in 3D mode, including the denser particle count and tighter interaction radius.

Useful overrides:

```powershell
cargo run --release -- --particle-count 36000 --attraction-strength 0.75
cargo run --release -- --curvature-scale 2.8 --breach-threshold-base 1.05
cargo run --release -- --density-reactivity 0.55 --density-force-gain 0.65
cargo run --release -- --particle-count 56000 --field-splat-radius 3 --flow-coupling 0.38 --collapse-gain 0.42
cargo run --release -- --merge-density-threshold 1.05 --merge-influence-gain 0.9 --attraction-ramp-gain 0.45 --micro-pop-threshold 0.34
cargo run --release -- --signal-speed 220 --shell-source-scale 0.68
cargo run --release -- --shell-overlap-gain 1.4 --spike-threshold 0.28 --shell-reform-gain 0.52
cargo run --release -- --mode 3d --particle-count-3d 140000 --field-splat-radius 2 --camera-distance 920 --projection-scale 640 --far-field-stride-3d 2 --camera-yaw-speed 0.004
cargo run --release -- --headless-steps 200
cargo run --release -- --show-debug-shells
```

Controls:

- `Esc`: quit

Headless mode is useful for fast validation and parameter sweeps without opening the renderer.

## Cone Test (Causal-Front Falsification)

The cone test turns the finite-speed claim of the ontology into a measurable experiment: inject one localized tension pulse into a quiescent substrate and measure the disturbance front over time. Two constitutive laws can be tested with `--cone-model`:

- `wave` (default): the conservative acoustic closure `dA/dt = -c*div(G)`, `dG/dt = -c*grad(A)`. Energy-conserving, real propagating wave; `signal_speed` is the wave speed `c`.
- `relax`: the current relaxational particle engine (diffusive; kept as a comparison baseline).

```powershell
cargo cone
# equivalently:
cargo run --release -- --cone-test
# compare against the relaxational engine:
cargo run --release -- --cone-test --cone-model relax
```

Useful overrides:

```powershell
cargo run --release -- --cone-test --cone-steps 300 --cone-perturb-amplitude 2.0
cargo run --release -- --cone-test --signal-speed 18 --cone-csv runs/wave_a.csv
```

Wave verdicts:

- **wave speed (|A|)** vs **c** — the amplitude front should radiate at ~`c`, so `signal_speed` is a real speed, not just the neighbor-list reach.
- **finite propagation** — PASS if nothing outruns the operator stencil. A faint dispersive precursor ahead of the front is an expected lattice artifact.
- **isotropy** — PASS if the wavefront is round.
- **energy** — bounded (no damping) confirms the medium is conservative rather than dissipative.

Relax verdicts (diffusive baseline): causal-horizon classification, causality, isotropy, and the `edge/causal` ratio (~0.02, showing the front barely moves). Per-step data is written to `cone_test.csv`. See [SimulationSpec.md](SimulationSpec.md) section 9 for the full methodology.

## First Experiments

1. Run the default setup and watch whether the `top` cavity breaches first while `up` remains more stable.
2. Raise `--curvature-scale` and compare how quickly shell collapse localizes around the smallest cavity.
3. Lower `--breach-threshold-base` to study whether outward tension packets start appearing in dense shell regions earlier.
4. Increase `--particle-count` and inspect whether the same behavior survives a finer substrate discretization.

## Notes

- True voids are represented as hard masked zero-tension cores.
- The visible field dots are substrate particles: discrete samples of the tensioned medium, not sensors.
- `--mode 3d` switches to a projected 3D point-cloud simulation with spherical voids and a rotating camera.
- The 3D renderer can be made larger and faster with `--projection-scale`, `--camera-distance`, and `--far-field-stride-3d` without lowering simulation density.
- Shells are encoded as explicit ring tension profiles derived from the ambient/core ratio.
- The old white ring was only a debug outline marking the nominal shell radius. It is now hidden by default and can be shown with `--show-debug-shells`.
- Local particle crowding now feeds back into shell reactivity, so bunching near a cavity can amplify local attraction and breach likelihood.
- When local nodes coalesce, attraction growth is now ramp-limited instead of jumping instantly; their influence spreads as a wider halo rather than a permanent spike.
- Compressed clumps store attraction energy and intermittently bleed it off as low-energy micro-pops, which is the intended source of the boiling-style flicker.
- Particle-to-particle attraction is now sourced from an effective shell radius around each local clump rather than from a point center, and it is limited by a finite causal horizon set by `signal_speed`.
- Overlapping local shells are now damped and converted into reforming surface-release behavior; when a compressed clump behaves more like a spike than a depression, it sheds energy back into the field instead of holding a persistent overlapping shell.
- Void-adjacent shell formation is now driven by missing interaction geometry around the void boundary rather than only by a hand-shaped ring profile.
- The renderer now draws particles as small field splats so the substrate appears as a denser sheet rather than isolated pixels.
- Neighbor velocity coupling and local collapse coupling help the substrate move more like a continuous medium when regions compress or rebound.
- Breach events are intentionally observable: the prototype prints metrics every 100 steps and visually highlights ejecta.
- The 3D extension path is straightforward: switch `Vector2<f32>` to `Vector3<f32>`, replace the 2D spatial hash with voxel hashing, and render slices or volume projections.