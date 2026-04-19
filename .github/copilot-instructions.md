# BedRock Guidelines

## Scope
- Work in Rust unless the user explicitly asks for another language.
- Treat BedRock as a scientific simulation codebase focused on emergent field behavior, particle methods, numerical update loops, and GPU-aware architecture.

## Modeling Standard
- Prefer falsifiable simulation design over visually plausible behavior.
- Distinguish core model rules from numerical safeguards, initialization bias, rendering choices, and heuristic patches.
- Make every non-physical assumption explicit, configurable, and easy to disable.
- Do not hardcode stabilizers, symmetry corrections, damping terms, or convenience constraints unless they are explicitly requested or clearly labeled as temporary scaffolding.
- Do not "fix" surprising behavior before checking whether it is a real consequence of the model.
- When the physics is underspecified, formalize the missing assumptions before implementation instead of silently inventing them.

## Implementation
- Prefer modular Cargo projects with separate physics, diagnostics, configuration, neighbor search, breach/event logic, and visualization boundaries when those concerns exist.
- Keep data layout and update loops performance-conscious, but do not optimize ahead of correctness and observability.
- Prefer CPU-first correctness with a credible GPU acceleration path.
- Use established Rust crates when they reduce implementation risk for math, parallelism, rendering, CLI handling, or GPU compute.

## Validation
- Add instrumentation so claims about emergent behavior can be inspected rather than assumed.
- Use comparative runs, ablations, and parameter sweeps when checking whether a behavior depends on hidden constraints.
- When reviewing code, prioritize correctness risks, numerical instability, accidental imposed behavior, and missing observability before style concerns.

## Build And Test
- Use Cargo commands for validation.
- Prefer `cargo check` for fast iteration, `cargo test` when tests exist, and `cargo run --release` for simulation execution.