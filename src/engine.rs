use crate::input::config::SimulationConfig;
use crate::rays::trace_ray;
use crate::ssp::init_ssp;
use crate::bty;
use crate::influence::{gaussian_beam_influence, init_pressure_field, PressureField};
use ndarray::Zip;
use rayon::prelude::*;
use std::f32::consts::PI;

fn merge_pressure_fields(dst: &mut PressureField, src: &PressureField) {
    Zip::from(&mut dst.pressure)
        .and(&src.pressure)
        .for_each(|d, s| *d += *s);

    Zip::from(&mut dst.delay_s)
        .and(&mut dst.amplitude)
        .and(&src.delay_s)
        .and(&src.amplitude)
        .for_each(|d_delay, d_amp, &s_delay, &s_amp| {
            if s_delay < *d_delay {
                *d_delay = s_delay;
                *d_amp = s_amp;
            }
        });
}

pub fn core(cfg: &SimulationConfig) -> (Option<Vec<Vec<[f32; 3]>>>, PressureField) {

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
    
    // Process rays in parallel with worker-local accumulation and reduce
    let store_ray_paths = cfg.beam.store_ray_paths;

    let (mut indexed_paths, pressure_field) = angle_pairs
        .par_iter()
        .enumerate()
        .map(|(idx, &(azim, elev))| {
            let mut local_pressure_field = init_pressure_field(cfg);
            
            // Trace ray and compute influence
            let mut ray_history = trace_ray(azim, elev, cfg, &ssp_field, &bty_field);
            gaussian_beam_influence(&mut ray_history, &mut local_pressure_field, 
                                   &bty_field, elev, d_azim, d_elev, &omega);
            
            // Extract ray path
            let path = if store_ray_paths {
                Some(
                    ray_history
                        .iter()
                        .map(|r| r.position)
                        .collect::<Vec<[f32; 3]>>(),
                )
            } else {
                None
            };
            
            (idx, path, local_pressure_field)
        })
        .fold(
            || (Vec::<(usize, Vec<[f32; 3]>)>::new(), init_pressure_field(cfg)),
            |mut acc, (idx, path, local_field)| {
                if let Some(path) = path {
                    acc.0.push((idx, path));
                }
                merge_pressure_fields(&mut acc.1, &local_field);
                acc
            },
        )
        .reduce(
            || (Vec::<(usize, Vec<[f32; 3]>)>::new(), init_pressure_field(cfg)),
            |mut left, right| {
                left.0.extend(right.0);
                merge_pressure_fields(&mut left.1, &right.1);
                left
            },
        );

    let ray_paths = if store_ray_paths {
        indexed_paths.sort_by_key(|(idx, _)| *idx);
        Some(
            indexed_paths
                .into_iter()
                .map(|(_, path)| path)
                .collect::<Vec<_>>(),
        )
    } else {
        None
    };

    return (ray_paths, pressure_field);
}
