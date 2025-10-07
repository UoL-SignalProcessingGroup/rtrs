use crate::input::config::SimulationConfig;
use crate::rays::trace_ray;
use crate::ssp::init_ssp;
use crate::bty;
use crate::influence::{gaussian_beam_influence, init_pressure_field, PressureField};
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

    // allocate ray paths and pressure field
    let mut ray_paths = Vec::with_capacity(launch_azim_rad.len() * launch_elev_rad.len());

    // initialize environmental fields
    let ssp_field = init_ssp(cfg);
    let bty_field = bty::init_bty(cfg);

    // allocate pressure field array
    let mut pressure_field = init_pressure_field(cfg);

    // angular frequency
    let omega: Vec<f32> = cfg.source.freq_hz.iter().map(|f| 2.0 * PI * f).collect();

    // loop over launch angles
    use std::time::Instant;

    let mut trace_durations: Vec<f64> = Vec::with_capacity(launch_azim_rad.len() * launch_elev_rad.len());
    let mut influence_durations: Vec<f64> = Vec::with_capacity(launch_azim_rad.len() * launch_elev_rad.len());

    for &azim in &launch_azim_rad {
        for &elev in &launch_elev_rad {
            // trace rays
            let t0 = Instant::now();
            let mut ray_history = trace_ray(azim, elev, cfg, &ssp_field, &bty_field);
            let t_trace = t0.elapsed().as_secs_f64();
            trace_durations.push(t_trace);

            // beam influence
            let t1 = Instant::now();
            gaussian_beam_influence(&mut ray_history, &mut pressure_field, elev, d_azim, d_elev, &omega);
            let t_infl = t1.elapsed().as_secs_f64();
            influence_durations.push(t_infl);

            // save ray path history for output
            let path = ray_history.iter().map(|r| r.position).collect::<Vec<[f32; 3]>>();
            ray_paths.push(path);
        }
    }

    // Print average timings
    if !trace_durations.is_empty() {
        let sum: f64 = trace_durations.iter().sum();
        let avg = sum / (trace_durations.len() as f64);
        println!("Average trace_ray time: {:.6} s over {} calls", avg, trace_durations.len());
        println!("  (total time {:.3} s)", sum);
    }
    if !influence_durations.is_empty() {
        let sum: f64 = influence_durations.iter().sum();
        let avg = sum / (influence_durations.len() as f64);
        println!("Average gaussian_beam_influence time: {:.6} s over {} calls", avg, influence_durations.len());
        println!("  (total time {:.3} s)", sum);
    }

    return (ray_paths, pressure_field);
}
