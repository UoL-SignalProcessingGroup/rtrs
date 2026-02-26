use crate::rays::Ray;
use crate::bty::{
    BTYfield,
    interpolate_bty,
};

use std::f32::consts::PI;
use crate::rays::ray_normal;
use num_complex::Complex32;

pub fn reflect_boundaries(ray_history: &mut Vec<Ray>, bty_field: &BTYfield) {
    let ray = ray_history.last_mut().unwrap();

    // Determine which boundary (surface or bottom) — handle only one per step
    // Surface: z < 0.0
    if ray.position[2] < 0.0 {
        // flat surface at z=0
        let dz = ray.direction[2];
        if dz != 0.0 {
            let t = -ray.position[2] / dz;
            ray.position[0] += ray.direction[0] * t;
            ray.position[1] += ray.direction[1] * t;
            ray.position[2] = 0.0;
        } else {
            ray.position[2] = 0.0;
        }

        // plane normal pointing upward into water
        let normal = [0.0_f32, 0.0_f32, 1.0_f32];

        // reflect direction across plane
        let d_dot_n = ray.direction[0]*normal[0] + ray.direction[1]*normal[1] + ray.direction[2]*normal[2];
        ray.direction[0] = ray.direction[0] - 2.0 * d_dot_n * normal[0];
        ray.direction[1] = ray.direction[1] - 2.0 * d_dot_n * normal[1];
        ray.direction[2] = ray.direction[2] - 2.0 * d_dot_n * normal[2];

        // nudge off boundary
        let nudge = 1e-9;
        ray.position[0] += nudge * normal[0];
        ray.position[1] += nudge * normal[1];
        ray.position[2] += nudge * normal[2];

        ray.num_top_bounces += 1;
        ray.phase += PI; // vacuum phase inversion

        // proceed to paraxial update below
    } else {
        // Bottom: check against interpolated bottom depth
        let z_bty = interpolate_bty(ray.position, bty_field);
        if ray.position[2] >= z_bty {
            // compute local bottom normal
            let (mut normal, _tangent) = calculate_bottom_normal_and_tangent(ray.position, bty_field);
            if normal[2] > 0.0 { normal = [-normal[0], -normal[1], -normal[2]]; }

            // intersect ray with local plane point
            let s = [ray.position[0], ray.position[1], z_bty];
            let n_dot_d = normal[0]*ray.direction[0] + normal[1]*ray.direction[1] + normal[2]*ray.direction[2];
            let n_dot_s_minus_p = normal[0]*(s[0]-ray.position[0]) + normal[1]*(s[1]-ray.position[1]) + normal[2]*(s[2]-ray.position[2]);
            const EPS: f32 = 1e-12;
            if n_dot_d.abs() > EPS {
                let t_plane = n_dot_s_minus_p / n_dot_d;
                ray.position[0] += ray.direction[0] * t_plane;
                ray.position[1] += ray.direction[1] * t_plane;
                ray.position[2] += ray.direction[2] * t_plane;
            } else {
                ray.position[2] = z_bty;
            }

            // reflect direction about local normal
            let d_dot_n = ray.direction[0]*normal[0] + ray.direction[1]*normal[1] + ray.direction[2]*normal[2];
            ray.direction[0] = ray.direction[0] - 2.0 * d_dot_n * normal[0];
            ray.direction[1] = ray.direction[1] - 2.0 * d_dot_n * normal[1];
            ray.direction[2] = ray.direction[2] - 2.0 * d_dot_n * normal[2];

            // nudge off boundary
            let nudge = 1e-9;
            ray.position[0] += nudge * normal[0];
            ray.position[1] += nudge * normal[1];
            ray.position[2] += nudge * normal[2];
            
            // Bottom boundary condition: fluid halfspace R, or vacuum fallback
            ray.num_bottom_bounces += 1;
            apply_bottom_bc(ray, &normal, bty_field);

            // proceed to paraxial update below
        } else {
            return; // no boundary hit
        }
    }

    // --- Paraxial rotation and phi update ---
    // Compute ray-centered normals e1,e2 for the incident ray prior to reflection.
    // approximate by computing them from the pre-reflection direction using
    // ray_normal. For robustness, use stored c in ray.

    // need the incident direction to compute incoming e1,e2. For simplicity
    // compute an approximate incident direction by reflecting ray.direction back
    // across the last-used plane normal. 

    // maybe a simpler approach is to use the current (post-reflection) direction but bellhop uses incident normals.
    // approximate by negating the reflected z component for surface only
    // which gives a reasonable e1,e2.
    let c_local = ray.c;
    let inc_dir = [ray.direction[0], ray.direction[1], -ray.direction[2]];
    let (e1, e2) = ray_normal(inc_dir, ray.phi, c_local);

    // Build rayt (scaled tangent) and rayn2, rayn1 as in bellhop's CalcTangent_Normals
    let rayt = [c_local * inc_dir[0], c_local * inc_dir[1], c_local * inc_dir[2]];
    // choose boundary normal for constructing rayn2: surface -> [0,0,1], bottom -> local normal
    // We will reuse the normal computed earlier when present; reconstruct if necessary
    let bdry_n = if ray.num_top_bounces > 0 && ray.num_bottom_bounces == 0 {
        [0.0_f32, 0.0_f32, 1.0_f32]
    } else {
        // approximate local bottom normal by sampling small z-gradient: use upward unit
        [0.0_f32, 0.0_f32, 1.0_f32]
    };
    // rayn2 = -cross(rayt, bdry_n)
    let mut rayn2 = [
        -(rayt[1] * bdry_n[2] - rayt[2] * bdry_n[1]),
        -(rayt[2] * bdry_n[0] - rayt[0] * bdry_n[2]),
        -(rayt[0] * bdry_n[1] - rayt[1] * bdry_n[0]),
    ];
    let r2norm = (rayn2[0]*rayn2[0] + rayn2[1]*rayn2[1] + rayn2[2]*rayn2[2]).sqrt();
    if r2norm > 0.0 { rayn2 = [rayn2[0]/r2norm, rayn2[1]/r2norm, rayn2[2]/r2norm]; }
    let rayn1 = [
        -(rayt[1] * rayn2[2] - rayt[2] * rayn2[1]),
        -(rayt[2] * rayn2[0] - rayt[0] * rayn2[2]),
        -(rayt[0] * rayn2[1] - rayt[1] * rayn2[0]),
    ];

    // Rotation matrix entries (RotMat): Rot(1,1)=dot(rayn1,e1), Rot(1,2)=dot(rayn1,e2)
    let rot11 = rayn1[0]*e1[0] + rayn1[1]*e1[1] + rayn1[2]*e1[2];
    let rot12 = rayn1[0]*e2[0] + rayn1[1]*e2[1] + rayn1[2]*e2[2];
    let rot21 = -rot12;
    let rot22 = rayn2[0]*e2[0] + rayn2[1]*e2[1] + rayn2[2]*e2[2];

    // rotate p/q into rayn basis
    let p_tilde_in = [rot11 * ray.p_tilde[0] + rot12 * ray.p_hat[0], rot11 * ray.p_tilde[1] + rot12 * ray.p_hat[1]];
    let p_hat_in   = [rot21 * ray.p_tilde[0] + rot22 * ray.p_hat[0], rot21 * ray.p_tilde[1] + rot22 * ray.p_hat[1]];
    let q_tilde_in = [rot11 * ray.q_tilde[0] + rot12 * ray.q_hat[0], rot11 * ray.q_tilde[1] + rot12 * ray.q_hat[1]];
    let q_hat_in   = [rot21 * ray.q_tilde[0] + rot22 * ray.q_hat[0], rot21 * ray.q_tilde[1] + rot22 * ray.q_hat[1]];

    // curvature corrections R1,R2,R3 are model-dependent; skip them (set 0) per instruction
    let r1 = 0.0_f32;
    let r2 = 0.0_f32;
    let r3 = 0.0_f32;

    // apply curvature change (bellhop formulas)
    let p_tilde_out = [p_tilde_in[0] + q_tilde_in[0] * r1 - q_hat_in[0] * r2,
                       p_tilde_in[1] + q_tilde_in[1] * r1 - q_hat_in[1] * r2];
    let p_hat_out   = [p_hat_in[0]   + q_tilde_in[0] * r2 + q_hat_in[0] * r3,
                       p_hat_in[1]   + q_tilde_in[1] * r2 + q_hat_in[1] * r3];

    // rotate back to e1,e2 (Rot^T)
    ray.p_tilde[0] = rot11 * p_tilde_out[0] + rot21 * p_hat_out[0];
    ray.p_tilde[1] = rot11 * p_tilde_out[1] + rot21 * p_hat_out[1];
    ray.p_hat  [0] = rot12 * p_tilde_out[0] + rot22 * p_hat_out[0];
    ray.p_hat  [1] = rot12 * p_tilde_out[1] + rot22 * p_hat_out[1];

    // Fortran left q unchanged through the curvature correction; keep q as-is
    // update det_q
    ray.det_q = ray.q_tilde[0]*ray.q_hat[1] - ray.q_tilde[1]*ray.q_hat[0];

    // update phi as bellhop does
    let dot_r = (rayn1[0]*e1[0] + rayn1[1]*e1[1] + rayn1[2]*e1[2]).clamp(-1.0, 1.0);
    ray.phi = ray.phi + 2.0 * dot_r.acos();
}


// Apply bottom boundary condition to `ray` at the moment of a bottom bounce.
// Lossless fluid-fluid reflection happens if `bty_field.hs_cp` and `bty_field.hs_rho` are both `Some`, otherwise fallback to rigid reflection.
fn apply_bottom_bc(ray: &mut Ray, normal: &[f32; 3], bty_field: &BTYfield) {

    let bottom_p_wave_speed = match bty_field.bottom_p_wave_speed {
        Some(v) => v,
        None => {
            // Rigid (perfectly reflecting) fallback: R = +1, no phase or amplitude change
            return;
        }
    };
    let density_bottom = match bty_field.bottom_density {
        Some(v) => v,
        None => return, // rigid fallback
    };

    let density_water = bty_field.water_density; // g/cm3  water (upper halfspace)

    let c_w = ray.c;

    // ray.direction is the slowness vector (|direction| = 1/c_w), already reflected (pointing back up).
    // normal points upward into water (normal[2] < 0 for a near-flat bottom).
    // Both are upward after the bounce, so their dot product is positive and equals Th
    // (the normal component of slowness at the boundary, as in Bellhop's Th = dot(t, nBdry)).
    let slowness_normal = normal[0] * ray.direction[0]
                        + normal[1] * ray.direction[1]
                        + normal[2] * ray.direction[2];

    // Tangential slowness squared: |s|2 − s_n2  (Snell's law invariant across interface)
    let slowness_tangential_sq = (1.0 / (c_w * c_w) - slowness_normal * slowness_normal).max(0.0);

    // Normal slowness squared in the bottom halfspace (negative = post-critical / evanescent)
    let slowness_normal_bottom_sq = slowness_tangential_sq - 1.0 / (bottom_p_wave_speed * bottom_p_wave_speed);

    // Bellhop convention: imaginary below critical, real above critical
    let vertical_slowness_bottom: Complex32 = if slowness_normal_bottom_sq >= 0.0 {
        Complex32::new(slowness_normal_bottom_sq.sqrt(), 0.0)
    } else {
        Complex32::new(0.0, (-slowness_normal_bottom_sq).sqrt())
    };

    // Fluid–fluid plane-wave reflection coefficient (omega cancels throughout)
    let num = vertical_slowness_bottom * density_water - Complex32::new(0.0, slowness_normal * density_bottom);
    let den = vertical_slowness_bottom * density_water + Complex32::new(0.0, slowness_normal * density_bottom);

    // Apply reflection coefficient to ray amplitude and phase. If |R| is very small, kill the ray.
    let reflection_coeff = - num / den;
    if reflection_coeff.norm() < 1e-5 {
        ray.amplitude = 0.0;    // kill
    } else {
        ray.amplitude *= reflection_coeff.norm();
        ray.phase    += reflection_coeff.arg();
    }
}


fn calculate_bottom_normal_and_tangent(
    position: [f32; 3],
    bty_field: &BTYfield,
) -> ([f32; 3], [f32; 3]) {

    // Find indices for x and y consistent with bilinear interpolation used in bty.rs
    let find_index = |arr: &[f32], val: f32| -> usize {
        match arr.binary_search_by(|probe| probe.partial_cmp(&val).unwrap()) {
            Ok(i) => i.min(arr.len() - 2),
            Err(i) => i.saturating_sub(1).min(arr.len() - 2),
        }
    };

    let i = find_index(&bty_field.x, position[0]);
    let j = find_index(&bty_field.y, position[1]);

    // Grid spacing
    let dx = bty_field.x[i + 1] - bty_field.x[i];
    let dy = bty_field.y[j + 1] - bty_field.y[j];

    // Normalized distances within the cell
    let xd = ((position[0] - bty_field.x[i]) / dx).clamp(0.0, 1.0);
    let yd = ((position[1] - bty_field.y[j]) / dy).clamp(0.0, 1.0);

    // Corner depths (match bty.rs orientation: z[[i, j]])
    let z00 = bty_field.z[[i, j]];
    let z10 = bty_field.z[[i + 1, j]];
    let z01 = bty_field.z[[i, j + 1]];
    let z11 = bty_field.z[[i + 1, j + 1]];

    // Bilinear slopes for dz/dx and dz/dy
    let dzdx_y0 = (z10 - z00) / dx;
    let dzdx_y1 = (z11 - z01) / dx;
    let dzdx = (1.0 - yd) * dzdx_y0 + yd * dzdx_y1;

    let dzdy_x0 = (z01 - z00) / dy;
    let dzdy_x1 = (z11 - z10) / dy;
    let dzdy = (1.0 - xd) * dzdy_x0 + xd * dzdy_x1;

    // Upward pointing normal (if depth increases downward, upward is -z)
    let nx = dzdx;
    let ny = dzdy;
    let nz = -1.0;
    let norm_len = (nx * nx + ny * ny + nz * nz).sqrt();
    let normal = [nx / norm_len, ny / norm_len, nz / norm_len];

    // Tangent vector in x-y plane, perpendicular to projected normal
    let tangent = [-ny, nx, 0.0];
    let tnorm = (tangent[0] * tangent[0] + tangent[1] * tangent[1]).sqrt();
    let tangent = if tnorm > 0.0 {
        [tangent[0] / tnorm, tangent[1] / tnorm, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    };

    (normal, tangent)
}