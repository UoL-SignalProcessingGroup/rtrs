mod input;
mod output;
mod rays;
mod ssp;
mod bty;
mod reflect;
mod influence;
mod utils;

use anyhow::Result;
use std::fs;

use input::config::SimulationConfig;
use output::write_hdf5;
use rays::trace_ray;
use ssp::init_ssp;
use influence::{
    hat_beam_influence,
    gaussian_beam_influence, 
    init_pressure_field, 
    PressureField
};

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

    let (ray_paths, pressure_field) = core(&config);

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


fn core(cfg: &SimulationConfig) -> (Vec<Vec<[f32; 3]>>, PressureField) {

    // convert angles to radians
    let launch_elev_rad: Vec<f32> = cfg.source.launch_elev_deg.iter().map(|d| d.to_radians()).collect();    // "alpha" in Fortran
    let launch_azim_rad: Vec<f32> = cfg.source.launch_azim_deg.iter().map(|d| d.to_radians()).collect();    // "beta" in Fortran
    let d_elev = if launch_elev_rad.len() >= 2 {
        (launch_elev_rad[1] - launch_elev_rad[0]).abs()
    } else {
        1.0_f32.to_radians()
    };
    let d_azim = if launch_azim_rad.len() >= 2 {
        (launch_azim_rad[1] - launch_azim_rad[0]).abs()
    } else {
        1.0_f32.to_radians()
    };

    // allocate ray paths and pressure field
    let mut ray_paths = Vec::with_capacity(launch_azim_rad.len() * launch_elev_rad.len());

    // initialize environmental fields
    let ssp_field = init_ssp(cfg);
    let bty_field = bty::init_bty(cfg);

    // allocate pressure field array
    let mut pressure_field = init_pressure_field(cfg);

    // angular frequency
    let omega = 2.0 * std::f32::consts::PI * cfg.source.freq_hz;


    // loop over launch angles
    for &azim in &launch_azim_rad {
        for &elev in &launch_elev_rad {
            
            // trace rays
            let mut ray_history = trace_ray(azim, elev, cfg, &ssp_field, &bty_field);

            // beam influence
            gaussian_beam_influence(&mut ray_history, &mut pressure_field, elev, d_azim, d_elev, omega);
            // hat_beam_influence(&mut ray_history, &mut pressure_field, elev, d_azim, d_elev, omega);

            // save ray path history for output
            let path = ray_history.iter().map(|r| r.position).collect::<Vec<[f32; 3]>>();
            ray_paths.push(path);
        }
    }
    // debug print
    // println!("pressure field {}", pressure_field.pressure);
    return (ray_paths, pressure_field);
        
}


