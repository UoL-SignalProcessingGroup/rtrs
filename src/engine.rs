use crate::input::config::SimulationConfig;
use crate::rays::trace_ray;
use crate::ssp::init_ssp;
use crate::bty;
use crate::influence::{gaussian_beam_influence, init_pressure_field, PressureField};
use rayon::prelude::*;
use std::f32::consts::PI;

pub fn core(cfg: &SimulationConfig) -> (Vec<Vec<[f32; 3]>>, PressureField) {

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

    // initialize environmental fields
    let ssp_field = init_ssp(cfg);
    let bty_field = bty::init_bty(cfg);

    // angular frequency
    let omega: Vec<f32> = cfg.source.freq_hz.iter().map(|f| 2.0 * PI * f).collect();

    // Create all angle pairs upfront for parallel processing
    let angle_pairs: Vec<(f32, f32)> = launch_azim_rad.iter()
        .flat_map(|&azim| {
            launch_elev_rad.iter().map(move |&elev| (azim, elev))
        })
        .collect();
    
    // Process rays in parallel - each thread gets its own local pressure field
    let results: Vec<_> = angle_pairs
        .par_iter()
        .map(|&(azim, elev)| {
            // Each thread initializes its own local pressure field
            let mut local_pressure_field = init_pressure_field(cfg);
            
            // Trace ray and compute influence
            let mut ray_history = trace_ray(azim, elev, cfg, &ssp_field, &bty_field);
            gaussian_beam_influence(&mut ray_history, &mut local_pressure_field, 
                                   elev, d_azim, d_elev, &omega);
            
            // Extract ray path
            let path = ray_history.iter()
                .map(|r| r.position)
                .collect::<Vec<[f32; 3]>>();
            
            (path, local_pressure_field)
        })
        .collect();
    
    // Separate ray paths and local pressure fields
    let (ray_paths, local_fields): (Vec<_>, Vec<_>) = results.into_iter().unzip();
    
    // Merge all local pressure fields into the final one
    let mut pressure_field = init_pressure_field(cfg);
    for local_field in local_fields {
        // Element-wise addition of pressure arrays
        pressure_field.pressure += &local_field.pressure;
        
        // Merge delay and amplitude: keep the earliest arrival with its amplitude
        for idx in ndarray::indices(pressure_field.delay_s.raw_dim()) {
            let local_delay = local_field.delay_s[idx];
            if local_delay < pressure_field.delay_s[idx] {
                pressure_field.delay_s[idx] = local_delay;
                pressure_field.amplitude[idx] = local_field.amplitude[idx];
            }
        }
    }

    return (ray_paths, pressure_field);
}
