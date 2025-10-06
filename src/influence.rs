use crate::input::config::SimulationConfig;
use crate::rays::{Ray};
use crate::utils::{
    dot,
    sub,
    
};

use std::f32::consts::PI;
use ndarray::Array4;
use num_complex::Complex32;
use std::cmp::Ordering;


pub struct PressureField {
    // 3D pressure field array [x, y, z] with complex values
    // now 4D: [freq, x, y, z]
    pub pressure: Array4<Complex32>,    // [nfreq, x,y,z]
    pub x_m: Vec<f32>,
    pub y_m: Vec<f32>,
    pub z_m: Vec<f32>,
}

pub fn init_pressure_field(config: &SimulationConfig) -> PressureField {

    // initialize empty pressure field
    let x_m = config.receivers.x_rcvr_m.clone();
    let y_m = config.receivers.y_rcvr_m.clone();
    let z_m = config.receivers.z_rcvr_m.clone();
    let nfreq = config.source.freq_hz.len();
    let shape = (nfreq, x_m.len(), y_m.len(), z_m.len());
    let pressure = Array4::<Complex32>::zeros(shape);

    // init struct 
    PressureField {
        pressure,
        x_m,
        y_m,
        z_m,
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
    elev: f32,
    d_azim: f32,
    d_elev: f32,
    omega: &Vec<f32>
) {
    let n_steps = ray_history.len();
    if n_steps < 2 { return; }

    // Constants and beam window similar to Fortran implementation
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

    // Pre-calc maximum gaussian radii for quick rejection (approx from q elements)
    let mut max_radius_a = vec![0.0f32; n_steps.saturating_sub(1)];
    let mut max_radius_b = vec![0.0f32; n_steps.saturating_sub(1)];
    for is in 0..n_steps - 1 {
        let q_t_a = ray_history[is].q_tilde;
        let q_t_b = ray_history[is + 1].q_tilde;
        let q_h_a = ray_history[is].q_hat;
        let q_h_b = ray_history[is + 1].q_hat;
        max_radius_a[is] = beam_window * q_h_a[1].abs().max(q_h_b[1].abs());
        max_radius_b[is] = beam_window * q_t_a[0].abs().max(q_t_b[0].abs());
    }

    // Process each receiver in the 3D grid
    // Iterate over ray segments first, then only visit receivers inside a per-segment AABB.
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

    for is in 1..n_steps {
        // Precompute per-segment values
        let ray_start = ray_history[is - 1].position;
        let ray_end = ray_history[is].position;
        let ray_vec = sub(&ray_end, &ray_start);
        let ray_length_sq = dot(&ray_vec, &ray_vec);
        // skip degenerate/very short segments
        if ray_length_sq < 1e-12 { continue; }
        let seg_len = ray_length_sq.sqrt();
        if seg_len <= 1e-4 { continue; }

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

        let x_range = match find_index_range(&pressure_field.x_m, min_x, max_x) { None => continue, Some(r) => r };
        let y_range = match find_index_range(&pressure_field.y_m, min_y, max_y) { None => continue, Some(r) => r };
        let z_range = match find_index_range(&pressure_field.z_m, min_z, max_z) { None => continue, Some(r) => r };

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
                    let mut t = dot(&to_receiver, &ray_vec) / ray_length_sq;
                    t = t.max(0.0).min(1.0);
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
                    let m = dot(&offset, &e1).abs();
                    let n = dot(&offset, &e2).abs();

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

                    // Loop over frequencies and accumulate into 4D pressure field
                    for (ifreq, &om) in omega.iter().enumerate() {
                        let phase = om * delay - phase_int;
                        let (s, c) = phase.sin_cos();
                        let contribution = Complex32::new(amp * c, amp * s);
                        pressure_field.pressure[[ifreq, ix, iy, iz]] += contribution;
                    }
                }
            }
        }
    }
}




// old
// pub fn hat_beam_influence(
//     ray_history: &mut Vec<Ray>, 
//     pressure_field: &mut PressureField,
//     elev: f32,
//     d_azim: f32,
//     d_elev: f32,
//     omega: f32
// ) {
//     let n_steps = ray_history.len();
//     if n_steps < 2 { return; }

//     // Scale the beam
//     scale_beam(elev, d_elev, d_azim, ray_history);
    
//     // Calculate KMAH phase shifts at caustics
//     let mut kmah_phase = vec![0.0; n_steps];
//     for is in 1..n_steps {
//         kmah_phase[is] = kmah_phase[is - 1];
//         let det_q_curr = ray_history[is].det_q;
//         let det_q_prev = ray_history[is - 1].det_q;
        
//         if (det_q_curr <= 0.0 && det_q_prev > 0.0) || (det_q_curr >= 0.0 && det_q_prev < 0.0) {
//             kmah_phase[is] = kmah_phase[is - 1] + PI / 2.0;
//         }
//     }

//     // Pre-calculate tangents and normals for each ray step
//     let mut e1_g = vec![[0.0; 3]; n_steps];
//     let mut e2_g = vec![[0.0; 3]; n_steps];
    
//     for is in 0..n_steps {
//         let ray = &ray_history[is];
//         (e1_g[is], e2_g[is]) = crate::rays::ray_normal(ray.direction, ray.phi, ray.c);
//     }

//     // Process each receiver in the 3D grid
//     for (ix, &x_rcvr) in pressure_field.x_m.iter().enumerate() {
//         for (iy, &y_rcvr) in pressure_field.y_m.iter().enumerate() {
//             for (iz, &z_rcvr) in pressure_field.z_m.iter().enumerate() {
                
//                 // Vector from receiver to each ray step
//                 let mut xt = vec![[0.0; 3]; n_steps];
//                 for is in 0..n_steps {
//                     xt[is] = [
//                         ray_history[is].position[0] - x_rcvr,
//                         ray_history[is].position[1] - y_rcvr,
//                         ray_history[is].position[2] - z_rcvr,
//                     ];
//                 }

//                 // Cross products for each step
//                 let mut xtxe1 = vec![[0.0; 3]; n_steps];
//                 let mut xtxe2 = vec![[0.0; 3]; n_steps];
//                 for is in 0..n_steps {
//                     xtxe1[is] = cross_product(&xt[is], &e1_g[is]);
//                     xtxe2[is] = cross_product(&xt[is], &e2_g[is]);
//                 }

//                 // Step through the ray
//                 for is in 1..n_steps {
//                     // Skip duplicate points (boundary reflections)
//                     let ray_segment_length = norm(&sub(&ray_history[is].position, &ray_history[is - 1].position));
//                     if ray_segment_length <= 1e-4 {
//                         continue;
//                     }

//                     // Check if the ray passes close to this receiver
//                     let ray_start = ray_history[is - 1].position;
//                     let ray_end = ray_history[is].position;
//                     let receiver_pos = [x_rcvr, y_rcvr, z_rcvr];
                    
//                     // Find closest point on ray segment to receiver
//                     let ray_vec = sub(&ray_end, &ray_start);
//                     let to_receiver = sub(&receiver_pos, &ray_start);
//                     let ray_length_sq = dot(&ray_vec, &ray_vec);
                    
//                     if ray_length_sq < 1e-12 { continue; } // Degenerate segment
                    
//                     let t = dot(&to_receiver, &ray_vec) / ray_length_sq;
//                     let t = t.max(0.0).min(1.0); // Clamp to segment
                    
//                     // Closest point on ray segment
//                     let closest_point = [
//                         ray_start[0] + t * ray_vec[0],
//                         ray_start[1] + t * ray_vec[1],
//                         ray_start[2] + t * ray_vec[2],
//                     ];
                    
//                     // Distance from receiver to ray
//                     let dist_to_ray = norm(&sub(&receiver_pos, &closest_point));
                    
//                     // Only process if receiver is within reasonable distance
//                     // (this is a simplified beam width check)
//                     if dist_to_ray > 100.0 { continue; } // 100m beam width for now
                    
//                     // Linear interpolation of ray properties at closest point
//                     let q_tilde = [
//                         ray_history[is - 1].q_tilde[0] + t * (ray_history[is].q_tilde[0] - ray_history[is - 1].q_tilde[0]),
//                         ray_history[is - 1].q_tilde[1] + t * (ray_history[is].q_tilde[1] - ray_history[is - 1].q_tilde[1]),
//                     ];
//                     let q_hat = [
//                         ray_history[is - 1].q_hat[0] + t * (ray_history[is].q_hat[0] - ray_history[is - 1].q_hat[0]),
//                         ray_history[is - 1].q_hat[1] + t * (ray_history[is].q_hat[1] - ray_history[is - 1].q_hat[1]),
//                     ];
                    
//                     // Determinant of the ray tube
//                     let det_q_int = q_tilde[0] * q_hat[1] - q_hat[0] * q_tilde[1];
//                     if det_q_int.abs() < 1e-12 {
//                         continue; // Receiver outside beam or degenerate beam
//                     }
                    
//                     // Calculate beam coordinates
//                     // Project distance onto the two beam directions
//                     let (e1, e2) = crate::rays::ray_normal(ray_history[is].direction, ray_history[is].phi, ray_history[is].c);
//                     let offset = sub(&receiver_pos, &closest_point);
//                     let m = dot(&offset, &e1).abs();
//                     let n = dot(&offset, &e2).abs();

//                     // Beam coordinates in paraxial system
//                     let a = if q_hat[1].abs() > 1e-12 { 
//                         ((- q_hat[0] * m + q_hat[1] * n) / det_q_int).abs() 
//                     } else { 
//                         (m / det_q_int.abs()).abs() 
//                     };
//                     let b = if q_tilde[0].abs() > 1e-12 { 
//                         ((q_tilde[0] * m - q_tilde[1] * n) / det_q_int).abs() 
//                     } else { 
//                         (n / det_q_int.abs()).abs() 
//                     };
                    
//                     if a.max(b) > 1.0 {
//                         continue; // Receiver outside beam
//                     }

//                     // Hat function beam pattern
//                     let w = (1.0 - a) * (1.0 - b);
                    
//                     // Travel time to closest point
//                     let delay = ray_history[is - 1].travel_time + t * (ray_history[is].travel_time - ray_history[is - 1].travel_time);

//                     let const_amp = ray_history[is].amplitude / det_q_int.abs().sqrt();
//                     let amp = const_amp * w;

//                     // Phase shift at caustics
//                     let mut phase_int = ray_history[is - 1].phase + kmah_phase[is - 1];
//                     let det_q_prev = ray_history[is - 1].det_q;
//                     if (det_q_int <= 0.0 && det_q_prev > 0.0) || (det_q_int >= 0.0 && det_q_prev < 0.0) {
//                         phase_int += PI / 2.0;
//                     }

//                     // Apply coherent contribution
//                     let phase = omega * delay - phase_int;
//                     let contribution = Complex32::new(
//                         amp * phase.cos(),
//                         amp * phase.sin()
//                     );

//                     pressure_field.pressure[[ix, iy, iz]] += contribution;
//                 }
//             }
//         }
//     }
// }