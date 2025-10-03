use crate::rays::Ray;
use crate::bty::{
    BTYfield,
    interpolate_bty,
};

use std::f32::consts::PI;
use num_complex::Complex;

pub fn surface_reflection(ray_history: &mut Vec<Ray>) {
    let ray = ray_history.last_mut().unwrap();
    if ray.position[2] < 0.0 {
        let dz = ray.direction[2];
        if dz != 0.0 {
            let t = -ray.position[2] / dz;
            ray.position[0] += ray.direction[0] * t;
            ray.position[1] += ray.direction[1] * t;
            ray.position[2] = 0.0;
        }
        ray.direction[2] = -ray.direction[2];
        ray.num_top_bounces += 1;
        ray.phase += PI; // phase shift on reflection
        // amplitude stays the same
        // ray.travel_time += t;
    }
}

pub fn bottom_reflections(ray_history: &mut Vec<Ray>, bty_field: &BTYfield) {
    // acoustic vacuum pressure release bottom reflections 

    let ray = ray_history.last_mut().unwrap();

    // interpolate bottom depth at current (x,y)
    let z_bty = interpolate_bty(ray.position, bty_field);
    if ray.position[2] >= z_bty {
        // Compute local bottom normal (pointing upward, into water)
        let (mut normal, _tangent) = calculate_bottom_normal_and_tangent(ray.position, bty_field);
        // Ensure the normal points upward (negative z if depth increases downward)
        if normal[2] > 0.0 { normal = [-normal[0], -normal[1], -normal[2]]; }

        // Define a point on the plane at current (x,y)
        let s = [ray.position[0], ray.position[1], z_bty];

        // Ray-plane intersection: t = n·(s - p) / n·d
        let n_dot_d = normal[0] * ray.direction[0]
            + normal[1] * ray.direction[1]
            + normal[2] * ray.direction[2];
        let n_dot_s_minus_p = normal[0] * (s[0] - ray.position[0])
            + normal[1] * (s[1] - ray.position[1])
            + normal[2] * (s[2] - ray.position[2]);

        const EPS: f32 = 1e-12;
        if n_dot_d.abs() > EPS {
            let t_plane = n_dot_s_minus_p / n_dot_d;
            ray.position[0] += ray.direction[0] * t_plane;
            ray.position[1] += ray.direction[1] * t_plane;
            ray.position[2] += ray.direction[2] * t_plane;
        } else {
            // Nearly parallel to plane: snap to surface vertically at current (x,y)
            ray.position[2] = z_bty;
        }

        // Reflect direction across plane: d' = d - 2(d·n)n
        let d_dot_n = ray.direction[0] * normal[0]
            + ray.direction[1] * normal[1]
            + ray.direction[2] * normal[2];
        ray.direction[0] = ray.direction[0] - 2.0 * d_dot_n * normal[0];
        ray.direction[1] = ray.direction[1] - 2.0 * d_dot_n * normal[1];
        ray.direction[2] = ray.direction[2] - 2.0 * d_dot_n * normal[2];

        // Nudge the ray slightly off the boundary along the normal to avoid re-hitting
        let nudge = 1e-9;
        ray.position[0] += nudge * normal[0];
        ray.position[1] += nudge * normal[1];
        ray.position[2] += nudge * normal[2];

        // pressure release reflection
        // ray.num_bottom_bounces += 1;
        // ray.phase += PI; // phase shift on reflection

        // Fresnel reflection coefficient for fluid-fluid interface
        let u = {
            let mag = (ray.direction[0]*ray.direction[0]
            + ray.direction[1]*ray.direction[1]
            + ray.direction[2]*ray.direction[2]).sqrt();
            if mag > 0.0 {
            [ray.direction[0]/mag, ray.direction[1]/mag, ray.direction[2]/mag]
            } else {
            [0.0, 0.0, 0.0]
            }
        };
        let cos_th1 = (u[0]*normal[0] + u[1]*normal[1] + u[2]*normal[2]).abs();
        let sin_th1 = (1.0_f32 - cos_th1 * cos_th1).max(0.0).sqrt();

        let sin_th2 = (ray.c / bty_field.c) * sin_th1;
        let rho_ocean = 1.0;

        let refl_c = if sin_th2 > 1.0 {
            // total internal reflection -> magnitude 1, zero imaginary part here
            Complex::new(1.0_f32, 0.0_f32)
        } else {
            let cos_th2 = (1.0_f32 - sin_th2 * sin_th2).max(0.0).sqrt();
            let z1 = rho_ocean * ray.c;
            let z2 = bty_field.density * bty_field.c;
            let re = (z2 * cos_th1 - z1 * cos_th2) / (z2 * cos_th1 + z1 * cos_th2);
            Complex::new(re, 0.0_f32)
        };


        // update amplitude and phase
        ray.amplitude *= refl_c.norm();
        ray.phase += refl_c.arg();
        // println!("Bottom reflection: refl_c = {}, new amp = {}, new phase = {}", refl_c, ray.amplitude, ray.phase);
        
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