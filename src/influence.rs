use crate::input::config::SimulationConfig;
use crate::rays::{Ray};
use crate::bty::BTYfield;
use crate::reflect::compute_bottom_reflection_coefficient;
use crate::utils::{
    dot,
    sub,
};

use std::f32::consts::PI;
use ndarray::{Array4, Array3};
use num_complex::Complex32;
use std::cmp::Ordering;


pub struct PressureField {
    // 3D pressure field array [x, y, z] with complex values
    // now 4D: [freq, x, y, z]
    pub pressure: Array4<Complex32>,    // [nfreq, x,y,z]
    pub x_m: Vec<f32>,
    pub y_m: Vec<f32>,
    pub z_m: Vec<f32>,
    // when true, receivers are an explicit list of positions in `receiver_positions_m`
    pub is_array: bool,
    pub receiver_positions_m: Option<Vec<[f32; 3]>>,
    // per-receiver earliest arrival delay (s) and corresponding amplitude
    pub delay_s: Array3<f32>, // [x,y,z]
    pub amplitude: Array3<f32>, // [x,y,z]
}

pub fn init_pressure_field(config: &SimulationConfig) -> PressureField {

    // initialize empty pressure field
    let nfreq = config.source.freq_hz.len();

    // support two receiver configuration modes:
    // - "grid" (default): x_rcvr_m/y_rcvr_m/z_rcvr_m are axes -> produce 4D grid (nx,ny,nz)
    // - "array": x_rcvr_m/y_rcvr_m/z_rcvr_m are explicit per-receiver positions of equal length
    let x_m: Vec<f32>;
    let y_m: Vec<f32>;
    let z_m: Vec<f32>;
    let mut is_array = false;
    let mut receiver_positions_m: Option<Vec<[f32; 3]>> = None;
    let pressure;
    let delay_s;
    let amplitude;

    if config.receivers.config_type.to_lowercase() == "array" {
        // validate lengths match
        let nx = config.receivers.x_rcvr_m.len();
        let ny = config.receivers.y_rcvr_m.len();
        let nz = config.receivers.z_rcvr_m.len();
        assert!(nx == ny && ny == nz, "array receiver coordinates must have equal lengths");
        let nrec = nx;
        // create receiver positions vector
        let mut recs: Vec<[f32; 3]> = Vec::with_capacity(nrec);
        for i in 0..nrec {
            recs.push([config.receivers.x_rcvr_m[i], config.receivers.y_rcvr_m[i], config.receivers.z_rcvr_m[i]]);
        }
        receiver_positions_m = Some(recs);
        is_array = true;

        // keep x_m/y_m/z_m empty placeholders for array mode
        x_m = Vec::new(); y_m = Vec::new(); z_m = Vec::new();

        let shape = (nfreq, nrec, 1usize, 1usize);
        pressure = Array4::<Complex32>::zeros(shape);
        delay_s = Array3::<f32>::from_elem((nrec, 1, 1), f32::INFINITY);
        amplitude = Array3::<f32>::from_elem((nrec, 1, 1), 0.0_f32);
    } else {
        // grid mode (existing behavior)
        x_m = config.receivers.x_rcvr_m.clone();
        y_m = config.receivers.y_rcvr_m.clone();
        z_m = config.receivers.z_rcvr_m.clone();
        let shape = (nfreq, x_m.len(), y_m.len(), z_m.len());
        pressure = Array4::<Complex32>::zeros(shape);
        // per-receiver arrays (nx, ny, nz)
        let nx = x_m.len();
        let ny = y_m.len();
        let nz = z_m.len();
        delay_s = Array3::<f32>::from_elem((nx, ny, nz), f32::INFINITY);
        amplitude = Array3::<f32>::from_elem((nx, ny, nz), 0.0_f32);
    }

    // init struct 
    PressureField {
        pressure,
        x_m,
        y_m,
        z_m,
        is_array,
        receiver_positions_m,
        delay_s,
        amplitude,
    }
}

fn scale_beam(elev: f32, d_elev: f32, d_azim: f32, ray: &mut [Ray]) {
    let c0 = ray[0].c;
    let ratio1 = (elev.cos().abs().sqrt()) * (d_elev*d_azim).sqrt() / c0 / (2.0*PI);
    for r in ray.iter_mut() { r.amplitude *= ratio1 * r.c; }
    // DetQ with ORIGINAL q’s (then scale q’s)
    for r in ray.iter_mut() {
        r.det_q = r.q_tilde[0]*r.q_hat[1] - r.q_tilde[1]*r.q_hat[0];
    }
    let s_a = d_elev / c0;
    let s_b = elev.cos().abs() * d_azim / c0;
    for r in ray.iter_mut() {
        r.q_tilde[0] *= s_a; r.q_tilde[1] *= s_a;
        r.q_hat  [0] *= s_b; r.q_hat  [1] *= s_b;
    }
    // ratio1
}


pub fn gaussian_beam_influence(
    ray_history: &mut Vec<Ray>, 
    pressure_field: &mut PressureField,
    bty_field: &BTYfield,
    elev: f32,
    d_azim: f32,
    d_elev: f32,
    omega: &Vec<f32>
) {
    let n_steps = ray_history.len();
    if n_steps < 2 { return; }

    // Constants and beam window similar to Bellhop implementation
    let beam_window = 4.0_f32; // kills beams outside exp(-0.5 * BeamWindow^2)

    // Scale the beam (note: scale_beam applies geometric scaling)
    scale_beam(elev, d_elev, d_azim, ray_history);

    // Compute KMAH phase shifts at caustics
    let mut kmah_phase = vec![0.0f32; n_steps];
    for is in 1..n_steps {
        kmah_phase[is] = kmah_phase[is - 1];
        let det_q_curr = ray_history[is].det_q;
        let det_q_prev = ray_history[is - 1].det_q;
        if (det_q_curr <= 0.0 && det_q_prev > 0.0) || (det_q_curr >= 0.0 && det_q_prev < 0.0) {
            kmah_phase[is] = kmah_phase[is - 1] + PI / 2.0;
        }
    }

    // Pre-calculate tangents and normals for each ray step
    let mut e1_g = vec![[0.0; 3]; n_steps];
    let mut e2_g = vec![[0.0; 3]; n_steps];
    for is in 0..n_steps {
        let ray = &ray_history[is];
        (e1_g[is], e2_g[is]) = crate::rays::ray_normal(ray.direction, ray.phi, ray.c);
    }

    let mut running_bottom_reflection = vec![Complex32::new(1.0, 0.0); omega.len()];

    // Pre-calc maximum gaussian radii for quick rejection (approx from q elements)
    let mut max_radius_a = vec![0.0f32; n_steps.saturating_sub(1)];
    let mut max_radius_b = vec![0.0f32; n_steps.saturating_sub(1)];
    for is in 0..n_steps - 1 {
        let q_t_a = ray_history[is].q_tilde;
        let q_t_b = ray_history[is + 1].q_tilde;
        let q_h_a = ray_history[is].q_hat;
        let q_h_b = ray_history[is + 1].q_hat;
        // compute vector norms L1 = ||q_tilde|| and L2 = ||q_hat|| (Euclidean)
        let l2_a = (q_h_a[0]*q_h_a[0] + q_h_a[1]*q_h_a[1]).sqrt();
        let l2_b = (q_h_b[0]*q_h_b[0] + q_h_b[1]*q_h_b[1]).sqrt();
        let l1_a = (q_t_a[0]*q_t_a[0] + q_t_a[1]*q_t_a[1]).sqrt();
        let l1_b = (q_t_b[0]*q_t_b[0] + q_t_b[1]*q_t_b[1]).sqrt();
        max_radius_a[is] = beam_window * l2_a.max(l2_b);
        max_radius_b[is] = beam_window * l1_a.max(l1_b);
    }

    // Process each receiver in the 3D grid
    // Iterate over ray segments first, then only visit receivers inside a per-segment AABB.
    // If pressure_field.is_array, treat receivers as explicit list in receiver_positions_m
    if pressure_field.is_array {
        let recs = pressure_field.receiver_positions_m.as_ref().expect("receiver_positions_m must be Some in array mode");
        let nrec = recs.len();

        for is in 1..n_steps {
            // Precompute per-segment values
            let ray_start = ray_history[is - 1].position;
            let ray_end = ray_history[is].position;
            let ray_vec = sub(&ray_end, &ray_start);
            let ray_length_sq = dot(&ray_vec, &ray_vec);
            if ray_length_sq >= 1e-12 {
                let seg_len = ray_length_sq.sqrt();
                if seg_len > 1e-4 {
                    // beam radius estimate (use precomputed maxima)
                    let mut beam_radius = 0.0_f32;
                    if is - 1 < max_radius_a.len() && is - 1 < max_radius_b.len() {
                        beam_radius = max_radius_a[is - 1].max(max_radius_b[is - 1]);
                    }
                    let half_len = 0.5 * seg_len;
                    let search_radius = beam_radius + half_len + 1e-6_f32;

                    // Use precomputed normals
                    let e1 = e1_g[is];
                    let e2 = e2_g[is];

                    for irec in 0..nrec {
                        let receiver_pos = recs[irec];

                        // Find closest point on ray segment to receiver
                        let to_receiver = sub(&receiver_pos, &ray_start);
                        let t_raw = dot(&to_receiver, &ray_vec) / ray_length_sq;
                        if t_raw < 0.0 || t_raw > 1.0 { continue; }
                        let t = t_raw;
                        let closest_point = [
                            ray_start[0] + t * ray_vec[0],
                            ray_start[1] + t * ray_vec[1],
                            ray_start[2] + t * ray_vec[2],
                        ];

                        let dx = receiver_pos[0] - closest_point[0];
                        let dy = receiver_pos[1] - closest_point[1];
                        let dz = receiver_pos[2] - closest_point[2];
                        let dist_sq = dx*dx + dy*dy + dz*dz;
                        if dist_sq > (search_radius * search_radius) { continue; }

                        // Linear interpolation of q's at closest point
                        let q_tilde = [
                            ray_history[is - 1].q_tilde[0] + t * (ray_history[is].q_tilde[0] - ray_history[is - 1].q_tilde[0]),
                            ray_history[is - 1].q_tilde[1] + t * (ray_history[is].q_tilde[1] - ray_history[is - 1].q_tilde[1]),
                        ];
                        let q_hat = [
                            ray_history[is - 1].q_hat[0] + t * (ray_history[is].q_hat[0] - ray_history[is - 1].q_hat[0]),
                            ray_history[is - 1].q_hat[1] + t * (ray_history[is].q_hat[1] - ray_history[is - 1].q_hat[1]),
                        ];

                        let det_q_int = q_tilde[0] * q_hat[1] - q_hat[0] * q_tilde[1];
                        if det_q_int.abs() < 1e-12 { continue; }

                        let offset = sub(&receiver_pos, &closest_point);
                        let n = dot(&offset, &e1).abs();
                        let m = dot(&offset, &e2).abs();

                        let a = if q_hat[1].abs() > 1e-12 {
                            ((- q_hat[0] * m + q_hat[1] * n) / det_q_int).abs()
                        } else {
                            (m / det_q_int.abs()).abs()
                        };
                        let b = if q_tilde[0].abs() > 1e-12 {
                            ((q_tilde[0] * m - q_tilde[1] * n) / det_q_int).abs()
                        } else {
                            (n / det_q_int.abs()).abs()
                        };

                        if a + b > beam_window { continue; }
                        let w = (-0.5_f32 * (a * a + b * b)).exp();
                        let delay = ray_history[is - 1].travel_time + t * (ray_history[is].travel_time - ray_history[is - 1].travel_time);
                        let const_amp = ray_history[is].amplitude / det_q_int.abs().sqrt();
                        let amp = const_amp * w;

                        let mut phase_int = ray_history[is - 1].phase + kmah_phase[is - 1];
                        let det_q_prev = ray_history[is - 1].det_q;
                        if (det_q_int <= 0.0 && det_q_prev > 0.0) || (det_q_int >= 0.0 && det_q_prev < 0.0) {
                            phase_int += PI / 2.0;
                        }

                        // update earliest arrival for this receiver index
                        {
                            let cur_delay = pressure_field.delay_s[(irec, 0, 0)];
                            let cur_amp = pressure_field.amplitude[(irec, 0, 0)];
                            if delay < cur_delay || (delay == cur_delay && amp.abs() > cur_amp.abs()) {
                                pressure_field.delay_s[(irec, 0, 0)] = delay;
                                pressure_field.amplitude[(irec, 0, 0)] = amp;
                            }
                        }

                        // accumulate pressure for each frequency into (nfreq, irec, 0, 0)
                        for (ifreq, &om) in omega.iter().enumerate() {
                            let phase = om * delay - phase_int;
                            let (s, c) = phase.sin_cos();
                            let base_contribution = Complex32::new(amp * c, amp * s);
                            let contribution = running_bottom_reflection[ifreq] * base_contribution;
                            pressure_field.pressure[[ifreq, irec, 0, 0]] += contribution;
                        }
                    }
                }
            }

            if let Some(bottom_bounce) = ray_history[is].bottom_bounce {
                for (ifreq, &angular_frequency_rad_s) in omega.iter().enumerate() {
                    let reflection = compute_bottom_reflection_coefficient(
                        &bty_field.bottom_model,
                        bty_field.water_density_g_cm3,
                        bottom_bounce.water_sound_speed_m_s,
                        bottom_bounce.incident_slowness,
                        bottom_bounce.boundary_normal,
                        angular_frequency_rad_s,
                    );
                    running_bottom_reflection[ifreq] *= reflection;
                }
            }
        }

        return; // finished array-mode processing
    }
    let find_index_range = |arr: &Vec<f32>, min_v: f32, max_v: f32| -> Option<(usize, usize)> {
        if min_v > max_v { return None; }
        // lower bound
        let lo = match arr.binary_search_by(|v| v.partial_cmp(&min_v).unwrap_or(Ordering::Equal)) {
            Ok(i) => i,
            Err(i) => i,
        };
        // upper bound -> find last index <= max_v
        let hi = match arr.binary_search_by(|v| v.partial_cmp(&max_v).unwrap_or(Ordering::Equal)) {
            Ok(i) => i,
            Err(i) => { if i == 0 { return None; } else { i - 1 } }
        };
        if lo > hi || lo >= arr.len() { return None; }
        Some((lo, hi))
    };

    // grid mode
    for is in 1..n_steps {
        // Precompute per-segment values
        let ray_start = ray_history[is - 1].position;
        let ray_end = ray_history[is].position;
        let ray_vec = sub(&ray_end, &ray_start);
        let ray_length_sq = dot(&ray_vec, &ray_vec);
        // skip degenerate/very short segments
        if ray_length_sq >= 1e-12 {
            let seg_len = ray_length_sq.sqrt();
            if seg_len > 1e-4 {
                // beam radius estimate (use precomputed maxima)
                let mut beam_radius = 0.0_f32;
                if is - 1 < max_radius_a.len() && is - 1 < max_radius_b.len() {
                    beam_radius = max_radius_a[is - 1].max(max_radius_b[is - 1]);
                }
                // Use a conservative radius including half-segment length
                let half_len = 0.5 * seg_len;
                let search_radius = beam_radius + half_len + 1e-6_f32;

                // AABB of the segment inflated by search_radius
                let min_x = ray_start[0].min(ray_end[0]) - search_radius;
                let max_x = ray_start[0].max(ray_end[0]) + search_radius;
                let min_y = ray_start[1].min(ray_end[1]) - search_radius;
                let max_y = ray_start[1].max(ray_end[1]) + search_radius;
                let min_z = ray_start[2].min(ray_end[2]) - search_radius;
                let max_z = ray_start[2].max(ray_end[2]) + search_radius;

                let x_range = find_index_range(&pressure_field.x_m, min_x, max_x);
                let y_range = find_index_range(&pressure_field.y_m, min_y, max_y);
                let z_range = find_index_range(&pressure_field.z_m, min_z, max_z);

                if let (Some(x_range), Some(y_range), Some(z_range)) = (x_range, y_range, z_range) {
                    // Use precomputed normals
                    let e1 = e1_g[is];
                    let e2 = e2_g[is];

                    // iterate only candidate receivers
                    for ix in x_range.0..=x_range.1 {
                        let x_rcvr = pressure_field.x_m[ix];
                        for iy in y_range.0..=y_range.1 {
                            let y_rcvr = pressure_field.y_m[iy];
                            for iz in z_range.0..=z_range.1 {
                                let z_rcvr = pressure_field.z_m[iz];

                                let receiver_pos = [x_rcvr, y_rcvr, z_rcvr];

                                // Find closest point on ray segment to receiver
                                let to_receiver = sub(&receiver_pos, &ray_start);
                                let t_raw = dot(&to_receiver, &ray_vec) / ray_length_sq;
                                // reject contributions from projections outside the segment (no backward extrapolation)
                                if t_raw < 0.0 || t_raw > 1.0 { continue; }
                                let t = t_raw;
                                let closest_point = [
                                    ray_start[0] + t * ray_vec[0],
                                    ray_start[1] + t * ray_vec[1],
                                    ray_start[2] + t * ray_vec[2],
                                ];

                                // squared distance quick reject
                                let dx = receiver_pos[0] - closest_point[0];
                                let dy = receiver_pos[1] - closest_point[1];
                                let dz = receiver_pos[2] - closest_point[2];
                                let dist_sq = dx*dx + dy*dy + dz*dz;
                                if dist_sq > (search_radius * search_radius) { continue; }

                                // Linear interpolation of q's at closest point
                                let q_tilde = [
                                    ray_history[is - 1].q_tilde[0] + t * (ray_history[is].q_tilde[0] - ray_history[is - 1].q_tilde[0]),
                                    ray_history[is - 1].q_tilde[1] + t * (ray_history[is].q_tilde[1] - ray_history[is - 1].q_tilde[1]),
                                ];
                                let q_hat = [
                                    ray_history[is - 1].q_hat[0] + t * (ray_history[is].q_hat[0] - ray_history[is - 1].q_hat[0]),
                                    ray_history[is - 1].q_hat[1] + t * (ray_history[is].q_hat[1] - ray_history[is - 1].q_hat[1]),
                                ];

                                // Determinant of the ray tube
                                let det_q_int = q_tilde[0] * q_hat[1] - q_hat[0] * q_tilde[1];
                                if det_q_int.abs() < 1e-12 { continue; }

                                // Beam coordinates using local normals
                                let offset = sub(&receiver_pos, &closest_point);
                                let n = dot(&offset, &e1).abs();
                                let m = dot(&offset, &e2).abs();

                                let a = if q_hat[1].abs() > 1e-12 {
                                    ((- q_hat[0] * m + q_hat[1] * n) / det_q_int).abs()
                                } else {
                                    (m / det_q_int.abs()).abs()
                                };
                                let b = if q_tilde[0].abs() > 1e-12 {
                                    ((q_tilde[0] * m - q_tilde[1] * n) / det_q_int).abs()
                                } else {
                                    (n / det_q_int.abs()).abs()
                                };

                                if a + b > beam_window { continue; }

                                // Gaussian weight
                                let w = (-0.5_f32 * (a * a + b * b)).exp();

                                // Travel time to closest point
                                let delay = ray_history[is - 1].travel_time + t * (ray_history[is].travel_time - ray_history[is - 1].travel_time);

                                let const_amp = ray_history[is].amplitude / det_q_int.abs().sqrt();
                                let amp = const_amp * w;

                                // Phase shift at caustics
                                let mut phase_int = ray_history[is - 1].phase + kmah_phase[is - 1];
                                let det_q_prev = ray_history[is - 1].det_q;
                                if (det_q_int <= 0.0 && det_q_prev > 0.0) || (det_q_int >= 0.0 && det_q_prev < 0.0) {
                                    phase_int += PI / 2.0;
                                }

                                // Record per-receiver earliest-arrival delay and corresponding amplitude
                                // If this contribution arrives earlier than stored delay, replace; if equal delay, keep larger amplitude
                                {
                                    let cur_delay = pressure_field.delay_s[(ix, iy, iz)];
                                    let cur_amp = pressure_field.amplitude[(ix, iy, iz)];
                                    // choose update when delay is smaller (earlier) or if equal delay but amplitude larger
                                    if delay < cur_delay || (delay == cur_delay && amp.abs() > cur_amp.abs()) {
                                        pressure_field.delay_s[(ix, iy, iz)] = delay;
                                        pressure_field.amplitude[(ix, iy, iz)] = amp;
                                    }
                                }

                                // Loop over frequencies and accumulate into 4D pressure field
                                for (ifreq, &om) in omega.iter().enumerate() {
                                    let phase = om * delay - phase_int;
                                    let (s, c) = phase.sin_cos();
                                    let base_contribution = Complex32::new(amp * c, amp * s);
                                    let contribution = running_bottom_reflection[ifreq] * base_contribution;
                                    pressure_field.pressure[[ifreq, ix, iy, iz]] += contribution;
                                }
                            }
                        }
                    }
                }
            }
        }

        if let Some(bottom_bounce) = ray_history[is].bottom_bounce {
            for (ifreq, &angular_frequency_rad_s) in omega.iter().enumerate() {
                let reflection = compute_bottom_reflection_coefficient(
                    &bty_field.bottom_model,
                    bty_field.water_density_g_cm3,
                    bottom_bounce.water_sound_speed_m_s,
                    bottom_bounce.incident_slowness,
                    bottom_bounce.boundary_normal,
                    angular_frequency_rad_s,
                );
                running_bottom_reflection[ifreq] *= reflection;
            }
        }
    }
}