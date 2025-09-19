mod input;
mod output;
mod rays;
mod ssp;
mod reflect;

use anyhow::Result;
use std::fs;
use input::config::SimulationConfig;
use output::write_hdf5;
use rays::{trace_ray, Ray};
use ssp::init_ssp;

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

    let ray_paths = core(&config);

    // Derive output path: replace extension with .h5 (append if none)
    let out_path = {
        let mut p = std::path::PathBuf::from(in_path);
        p.set_extension("h5");
        p
    };

    // Write HDF5 output
    let out_path_str = out_path.to_str().expect("Invalid output path");
    write_hdf5(out_path_str, &config, ray_paths)?;
    Ok(())
}


fn core(cfg: &SimulationConfig) -> Vec<Vec<[f64; 3]>> {

    // convert angles to radians
    let launch_elev_rad: Vec<f64> = cfg.source.launch_elev_deg.iter().map(|d| d.to_radians()).collect();    // "alpha" in Fortran
    let launch_azim_rad: Vec<f64> = cfg.source.launch_azim_deg.iter().map(|d| d.to_radians()).collect();    // "beta" in Fortran

    // allocate pressure field array
    let n_r_rcvr = cfg.receivers.ranges_m.as_ref().map_or(0, |v| v.len());
    let n_bear_rcvr = cfg.receivers.bearings_deg.as_ref().map_or(0, |v| v.len());
    let n_depths_rcvr = cfg.receivers.depths_m.as_ref().map_or(0, |v| v.len());
    let mut ray_paths = Vec::new();
    let mut pressure_field_real = vec![0.0; n_r_rcvr * n_bear_rcvr * n_depths_rcvr];

    let ssp_field = init_ssp(cfg);

    // loop over launch angles
    for &azim in &launch_azim_rad {
        for &elev in &launch_elev_rad {
            
            // trace rays
            let ray_history = trace_ray(azim, elev, cfg, &ssp_field);

            // beam influence

            // accumulate pressure field

            // save ray path history for output
            let path = ray_history.iter().map(|r| r.position).collect::<Vec<[f64; 3]>>();
            // println!("{}", ray_history.len());
            // println!("{:?}", path);
            ray_paths.push(path);

            // break;
        }
    }
    return ray_paths;
        
}


