use crate::input::config::SimulationConfig;
use crate::rays::trace_ray;
use crate::ssp::init_ssp;
use crate::bty;
use crate::influence::{gaussian_beam_influence, init_pressure_field, PressureField};
use std::f32::consts::PI;
use std::time::Instant;

const TIMER: bool = false;

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

    // allocate ray paths and pressure field
    let mut ray_paths = Vec::with_capacity(launch_azim_rad.len() * launch_elev_rad.len());

    // initialize environmental fields
    let ssp_field = init_ssp(cfg);
    let bty_field = bty::init_bty(cfg);

    // allocate pressure field array
    let mut pressure_field = init_pressure_field(cfg);

    // angular frequency
    let omega: Vec<f32> = cfg.source.freq_hz.iter().map(|f| 2.0 * PI * f).collect();

    // timing containers (keep in outer scope so they're available after the loop)
    let mut trace_durations: Vec<f32> = if TIMER { Vec::with_capacity(launch_azim_rad.len() * launch_elev_rad.len()) } else { Vec::new() };
    let mut influence_durations: Vec<f32> = if TIMER { Vec::with_capacity(launch_azim_rad.len() * launch_elev_rad.len()) } else { Vec::new() };
    
    // loop over launch angles
    for &azim in &launch_azim_rad {
        for &elev in &launch_elev_rad {

            // run trace and influence; measure timings if enabled. Use a block expression
            // so we can return the populated ray history and still keep timings in scope.
            let ray_history = {
                if TIMER {
                    // trace rays (timed)
                    let t0 = Instant::now();
                    let mut rh = trace_ray(azim, elev, cfg, &ssp_field, &bty_field);
                    let t_trace = t0.elapsed().as_secs_f32();
                    trace_durations.push(t_trace);

                    // beam influence (timed)
                    let t1 = Instant::now();
                    gaussian_beam_influence(&mut rh, &mut pressure_field, elev, d_azim, d_elev, &omega);
                    let t_infl = t1.elapsed().as_secs_f32();
                    influence_durations.push(t_infl);

                    rh
                } else {
                    // no timing
                    let mut rh = trace_ray(azim, elev, cfg, &ssp_field, &bty_field);
                    gaussian_beam_influence(&mut rh, &mut pressure_field, elev, d_azim, d_elev, &omega);
                    rh
                }
            };

            // save ray path history for output
            let path = ray_history.iter().map(|r| r.position).collect::<Vec<[f32; 3]>>();
            ray_paths.push(path);
        }
    }

    if TIMER {
        // Print average timings
        if !trace_durations.is_empty() {
            let sum: f32 = trace_durations.iter().sum();
            let avg = sum / (trace_durations.len() as f32);
            println!("Average trace_ray time: {:.6} s over {} calls", avg, trace_durations.len());
            println!("  (total time {:.3} s)", sum);
        }
        if !influence_durations.is_empty() {
            let sum: f32 = influence_durations.iter().sum();
            let avg = sum / (influence_durations.len() as f32);
            println!("Average gaussian_beam_influence time: {:.6} s over {} calls", avg, influence_durations.len());
            println!("  (total time {:.3} s)", sum);
        }
    }
    

    return (ray_paths, pressure_field);
}
