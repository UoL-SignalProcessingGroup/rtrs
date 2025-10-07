mod input;
mod output;
mod rays;
mod ssp;
mod bty;
mod reflect;
mod influence;
mod utils;
mod engine;

use anyhow::Result;
use std::fs;

use input::config::SimulationConfig;
use output::write_hdf5;

fn load_config(path: &str) -> Result<SimulationConfig> {
    let text = fs::read_to_string(path)?;
    let cfg: SimulationConfig = serde_json::from_str(&text)?;
    Ok(cfg)
}

fn main() -> Result<()> {

    // Load config from JSON
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <config.json>", args.get(0).unwrap_or(&"rtrs".into()));
        return Ok(());
    }
    let in_path = &args[1];
    let config = load_config(in_path)?;

    let (ray_paths, pressure_field) = engine::core(&config);

    // Derive output path: replace extension with .h5 (append if none)
    let out_path = {
        let mut p = std::path::PathBuf::from(in_path);
        p.set_extension("h5");
        p
    };

    // Write HDF5 output
    let out_path_str = out_path.to_str().expect("Invalid output path");
    write_hdf5(out_path_str, &config, ray_paths, pressure_field)?;
    Ok(())
}





