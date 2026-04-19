mod breach_detection;
mod config;
mod gpu;
mod neighbor_search;
mod neighbor_search_3d;
mod physics;
mod physics3d;
mod visualization;
mod visualization3d;

use clap::Parser;

use config::{Config, SimulationMode};
use physics::Simulation;
use physics3d::Simulation3d;
use visualization::Renderer;
use visualization3d::Renderer3d;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut config = Config::parse();
    config.apply_runtime_defaults();

    println!("BedRock prototype initialized.");
    println!("GPU path: {}", gpu::acceleration_notes());

    match config.mode {
        SimulationMode::TwoD => run_2d(config),
        SimulationMode::ThreeD => run_3d(config),
    }
}

fn run_2d(config: Config) -> Result<(), Box<dyn std::error::Error>> {
    let mut simulation = Simulation::new(config.clone());

    if config.headless_steps > 0 {
        for _ in 0..config.headless_steps {
            let metrics = simulation.update();
            if simulation.step % 100 == 0 || simulation.step == config.headless_steps as u64 {
                log_metrics_2d(&simulation, &metrics);
            }
        }
        return Ok(());
    }

    let mut renderer = Renderer::new(config.width, config.height)?;

    while renderer.is_open() {
        let mut metrics = None;

        for _ in 0..config.substeps_per_frame {
            metrics = Some(simulation.update());
        }

        let metrics = metrics.expect("substeps_per_frame must be greater than zero");
        renderer.draw(&simulation, &metrics)?;

        if simulation.step % 100 == 0 {
            log_metrics_2d(&simulation, &metrics);
        }
    }

    Ok(())
}

fn run_3d(config: Config) -> Result<(), Box<dyn std::error::Error>> {
    let mut simulation = Simulation3d::new(config.clone());

    if config.headless_steps > 0 {
        for _ in 0..config.headless_steps {
            let metrics = simulation.update();
            if simulation.step % 100 == 0 || simulation.step == config.headless_steps as u64 {
                log_metrics_3d(&simulation, &metrics);
            }
        }
        return Ok(());
    }

    let mut renderer = Renderer3d::new(config.width, config.height)?;

    while renderer.is_open() {
        let mut metrics = None;
        for _ in 0..config.substeps_per_frame {
            metrics = Some(simulation.update());
        }

        let metrics = metrics.expect("substeps_per_frame must be greater than zero");
        renderer.draw(&simulation, &metrics)?;

        if simulation.step % 100 == 0 {
            log_metrics_3d(&simulation, &metrics);
        }
    }

    Ok(())
}

fn log_metrics_2d(simulation: &Simulation, metrics: &physics::SimulationMetrics) {
    let shell_summary = metrics
        .shell_radii
        .iter()
        .map(|(class, radius)| format!("{}={:.2}", class.label(), radius))
        .collect::<Vec<_>>()
        .join(", ");

    println!(
        "step={} total_tension={:.3} true_voids={} pops_this_frame={} shell_radii=[{}]",
        simulation.step,
        metrics.total_tension,
        metrics.true_voids,
        metrics.pops_this_frame,
        shell_summary,
    );
}

fn log_metrics_3d(simulation: &Simulation3d, metrics: &physics3d::SimulationMetrics3d) {
    let shell_summary = metrics
        .shell_radii
        .iter()
        .map(|(class, radius)| format!("{}={:.2}", class.label(), radius))
        .collect::<Vec<_>>()
        .join(", ");

    println!(
        "step={} total_tension={:.3} true_voids={} micro_pops={} shell_radii=[{}]",
        simulation.step,
        metrics.total_tension,
        metrics.true_voids,
        metrics.micro_pops,
        shell_summary,
    );
}