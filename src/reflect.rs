use crate::bty::{
    bottom_normal_at_with_cursor,
    interpolate_bty_with_cursor,
    BottomBoundaryRuntimeModel,
    BTYCursor,
    BTYfield,
};
use crate::rays::{ray_normal, BottomBounceMetadata, Ray};
use num_complex::Complex32;
use std::f32::consts::PI;

pub fn reflect_boundaries(ray_history: &mut Vec<Ray>, bty_field: &BTYfield, bty_cursor: &mut BTYCursor) {
    let ray = ray_history.last_mut().unwrap();
    ray.bottom_bounce = None;

    let boundary_normal_for_paraxial: [f32; 3];
    let incident_slowness_for_paraxial: [f32; 3];

    // Determine which boundary (surface or bottom) — handle only one per step
    // Surface: z < 0.0
    if ray.position[2] <= 0.0 && ray.direction[2] < 0.0 {
        incident_slowness_for_paraxial = ray.direction;

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

        // outward normal from water at the top boundary (positive-down convention)
        let normal = [0.0_f32, 0.0_f32, -1.0_f32];
        boundary_normal_for_paraxial = normal;

        // reflect direction across plane
        let d_dot_n = ray.direction[0]*normal[0] + ray.direction[1]*normal[1] + ray.direction[2]*normal[2];
        ray.direction[0] = ray.direction[0] - 2.0 * d_dot_n * normal[0];
        ray.direction[1] = ray.direction[1] - 2.0 * d_dot_n * normal[1];
        ray.direction[2] = ray.direction[2] - 2.0 * d_dot_n * normal[2];

        // nudge off boundary
        let nudge = 1e-9;
        ray.position[0] -= nudge * normal[0];
        ray.position[1] -= nudge * normal[1];
        ray.position[2] -= nudge * normal[2];

        ray.num_top_bounces += 1;
        ray.phase += PI; // vacuum phase inversion

        // proceed to paraxial update below
    } else {
        // Bottom: check against interpolated bottom depth
        let z_bty = interpolate_bty_with_cursor(ray.position, bty_field, bty_cursor);
        if ray.position[2] >= z_bty && ray.direction[2] > 0.0 {
            let incident_slowness = ray.direction;
            incident_slowness_for_paraxial = incident_slowness;

            // compute local bottom normal
            let (mut normal, _tangent) = bottom_normal_at_with_cursor(ray.position, bty_field, bty_cursor);
            if normal[2] < 0.0 { normal = [-normal[0], -normal[1], -normal[2]]; }
            boundary_normal_for_paraxial = normal;

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
            ray.position[0] -= nudge * normal[0];
            ray.position[1] -= nudge * normal[1];
            ray.position[2] -= nudge * normal[2];
            
            ray.num_bottom_bounces += 1;
            ray.bottom_bounce = Some(BottomBounceMetadata {
                boundary_normal: normal,
                incident_slowness,
                water_sound_speed_m_s: ray.c,
            });

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
    // use the true incident slowness direction captured before reflection.
    let c_local = ray.c;
    let inc_dir = incident_slowness_for_paraxial;
    let (e1, e2) = ray_normal(inc_dir, ray.phi, c_local);

    // Build rayt (scaled tangent) and rayn2, rayn1 as in bellhop's CalcTangent_Normals
    let rayt = [c_local * inc_dir[0], c_local * inc_dir[1], c_local * inc_dir[2]];
    // choose boundary normal for constructing rayn2: surface -> [0,0,1], bottom -> local normal
    // We will reuse the normal computed earlier when present; reconstruct if necessary
    let bdry_n = boundary_normal_for_paraxial;
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

pub fn compute_bottom_reflection_coefficient(
    bottom_model: &BottomBoundaryRuntimeModel,
    water_density_g_cm3: f32,
    water_sound_speed_m_s: f32,
    incident_slowness: [f32; 3],
    boundary_normal: [f32; 3],
    angular_frequency_rad_s: f32,
) -> Complex32 {
    // No frequency means no dispersive/lossy boundary correction.
    if angular_frequency_rad_s <= 0.0 {
        return Complex32::new(1.0, 0.0);
    }

    // Rigid boundary: unity reflection with zero phase shift.
    if matches!(bottom_model, BottomBoundaryRuntimeModel::Rigid) {
        return Complex32::new(1.0, 0.0);
    }

    let normal_length = (boundary_normal[0] * boundary_normal[0]
        + boundary_normal[1] * boundary_normal[1]
        + boundary_normal[2] * boundary_normal[2])
        .sqrt();
    if normal_length <= 1e-12 {
        return Complex32::new(1.0, 0.0);
    }
    let unit_normal = [
        boundary_normal[0] / normal_length,
        boundary_normal[1] / normal_length,
        boundary_normal[2] / normal_length,
    ];

    // Re-normalize slowness to 1/c to limit integration drift.
    let expected_slowness_magnitude = 1.0 / water_sound_speed_m_s;
    let incident_slowness_magnitude = (incident_slowness[0] * incident_slowness[0]
        + incident_slowness[1] * incident_slowness[1]
        + incident_slowness[2] * incident_slowness[2])
        .sqrt();
    let slowness_scale = if incident_slowness_magnitude > 1e-12 {
        expected_slowness_magnitude / incident_slowness_magnitude
    } else {
        1.0
    };
    let normalized_incident_slowness = [
        incident_slowness[0] * slowness_scale,
        incident_slowness[1] * slowness_scale,
        incident_slowness[2] * slowness_scale,
    ];

    let normal_slowness = normalized_incident_slowness[0] * unit_normal[0]
        + normalized_incident_slowness[1] * unit_normal[1]
        + normalized_incident_slowness[2] * unit_normal[2];

    // Tangential slowness is Snell-invariant and determines kx.
    let tangential_slowness = [
        normalized_incident_slowness[0] - normal_slowness * unit_normal[0],
        normalized_incident_slowness[1] - normal_slowness * unit_normal[1],
        normalized_incident_slowness[2] - normal_slowness * unit_normal[2],
    ];
    let tangential_slowness_magnitude = (tangential_slowness[0] * tangential_slowness[0]
        + tangential_slowness[1] * tangential_slowness[1]
        + tangential_slowness[2] * tangential_slowness[2])
        .sqrt();

    let horizontal_wavenumber = angular_frequency_rad_s * tangential_slowness_magnitude;
    let vertical_wavenumber_water = angular_frequency_rad_s * normal_slowness;
    let imaginary_unit = Complex32::new(0.0, 1.0);

    let reflection_coefficient = match bottom_model {
        BottomBoundaryRuntimeModel::Rigid => Complex32::new(1.0, 0.0),
        BottomBoundaryRuntimeModel::Acoustic {
            compressional_speed_m_s,
            density_g_cm3,
            compressional_attenuation_db_per_wavelength,
        } => {
            // Convert user attenuation units (dB/lambda) to Np/m at this frequency.
            let compressional_attenuation_np_per_m =
                attenuation_db_per_wavelength_to_np_per_m(
                    *compressional_attenuation_db_per_wavelength,
                    *compressional_speed_m_s,
                    angular_frequency_rad_s,
                );
            let bottom_compressional_speed = complex_wave_speed_from_attenuation(
                *compressional_speed_m_s,
                compressional_attenuation_np_per_m,
                angular_frequency_rad_s,
            );
            let vertical_wavenumber_bottom_p = stable_complex_sqrt(
                Complex32::new(horizontal_wavenumber * horizontal_wavenumber, 0.0)
                    - complex_omega_over_speed_sq(angular_frequency_rad_s, bottom_compressional_speed),
            );

            // Bellhop fluid-halfspace form using f (vertical wavenumber) and g (density).
            let f_term = vertical_wavenumber_bottom_p;
            let g_term = Complex32::new(*density_g_cm3, 0.0);

            let numerator = Complex32::new(water_density_g_cm3, 0.0) * f_term
                - imaginary_unit * vertical_wavenumber_water * g_term;
            let denominator = Complex32::new(water_density_g_cm3, 0.0) * f_term
                + imaginary_unit * vertical_wavenumber_water * g_term;

            if denominator.norm() <= 1e-12 {
                Complex32::new(1.0, 0.0)
            } else {
                -numerator / denominator
            }
        }
        BottomBoundaryRuntimeModel::Elastic {
            compressional_speed_m_s,
            shear_speed_m_s,
            density_g_cm3,
            compressional_attenuation_db_per_wavelength,
            shear_attenuation_db_per_wavelength,
        } => {
            // Convert both P and S attenuation values to Np/m.
            let compressional_attenuation_np_per_m =
                attenuation_db_per_wavelength_to_np_per_m(
                    *compressional_attenuation_db_per_wavelength,
                    *compressional_speed_m_s,
                    angular_frequency_rad_s,
                );
            let shear_attenuation_np_per_m = attenuation_db_per_wavelength_to_np_per_m(
                *shear_attenuation_db_per_wavelength,
                *shear_speed_m_s,
                angular_frequency_rad_s,
            );

            let bottom_compressional_speed = complex_wave_speed_from_attenuation(
                *compressional_speed_m_s,
                compressional_attenuation_np_per_m,
                angular_frequency_rad_s,
            );
            let bottom_shear_speed = complex_wave_speed_from_attenuation(
                *shear_speed_m_s,
                shear_attenuation_np_per_m,
                angular_frequency_rad_s,
            );

            let horizontal_wavenumber_sq = Complex32::new(horizontal_wavenumber * horizontal_wavenumber, 0.0);
            let vertical_wavenumber_bottom_s_sq =
                horizontal_wavenumber_sq - complex_omega_over_speed_sq(angular_frequency_rad_s, bottom_shear_speed);
            let vertical_wavenumber_bottom_p_sq =
                horizontal_wavenumber_sq - complex_omega_over_speed_sq(angular_frequency_rad_s, bottom_compressional_speed);

            let vertical_wavenumber_bottom_s = stable_complex_sqrt(vertical_wavenumber_bottom_s_sq);
            let vertical_wavenumber_bottom_p = stable_complex_sqrt(vertical_wavenumber_bottom_p_sq);

            // Bellhop elastic halfspace terms (f, g) with shear-coupled boundary conditions.
            let shear_modulus = Complex32::new(*density_g_cm3, 0.0) * bottom_shear_speed * bottom_shear_speed;

            let y2 = ((vertical_wavenumber_bottom_s_sq + horizontal_wavenumber_sq)
                * (vertical_wavenumber_bottom_s_sq + horizontal_wavenumber_sq)
                - Complex32::new(4.0, 0.0)
                    * vertical_wavenumber_bottom_s
                    * vertical_wavenumber_bottom_p
                    * horizontal_wavenumber_sq)
                * shear_modulus;
            let y4 = vertical_wavenumber_bottom_p
                * (horizontal_wavenumber_sq - vertical_wavenumber_bottom_s_sq);

            let f_term = Complex32::new(angular_frequency_rad_s * angular_frequency_rad_s, 0.0) * y4;
            let g_term = y2;

            let numerator = Complex32::new(water_density_g_cm3, 0.0) * f_term
                - imaginary_unit * vertical_wavenumber_water * g_term;
            let denominator = Complex32::new(water_density_g_cm3, 0.0) * f_term
                + imaginary_unit * vertical_wavenumber_water * g_term;

            if denominator.norm() <= 1e-12 {
                Complex32::new(1.0, 0.0)
            } else {
                -numerator / denominator
            }
        }
    };

    if reflection_coefficient.norm() < 1e-12 {
        Complex32::new(0.0, 0.0)
    } else {
        reflection_coefficient
    }
}

fn attenuation_db_per_wavelength_to_np_per_m(
    attenuation_db_per_wavelength: f32,
    phase_speed_m_s: f32,
    angular_frequency_rad_s: f32,
) -> f32 {
    if attenuation_db_per_wavelength <= 0.0 || phase_speed_m_s <= 0.0 || angular_frequency_rad_s <= 0.0 {
        return 0.0;
    }
    let frequency_hz = angular_frequency_rad_s / (2.0 * PI);
    let wavelength_m = phase_speed_m_s / frequency_hz;
    let nepers_per_wavelength = attenuation_db_per_wavelength * (10.0_f32.ln() / 20.0);
    nepers_per_wavelength / wavelength_m
}

fn complex_omega_over_speed_sq(angular_frequency_rad_s: f32, wave_speed_m_s: Complex32) -> Complex32 {
    // Convenience helper for (omega / c)^2 in complex arithmetic.
    let omega_complex = Complex32::new(angular_frequency_rad_s, 0.0);
    let ratio = omega_complex / wave_speed_m_s;
    ratio * ratio
}

fn complex_wave_speed_from_attenuation(
    real_speed_m_s: f32,
    attenuation_np_per_m: f32,
    angular_frequency_rad_s: f32,
) -> Complex32 {
    // Build complex phase speed from real speed and attenuation (imaginary wavenumber part).
    if attenuation_np_per_m <= 0.0 || angular_frequency_rad_s <= 0.0 {
        return Complex32::new(real_speed_m_s, 0.0);
    }
    let real_wavenumber = angular_frequency_rad_s / real_speed_m_s;
    let complex_wavenumber = Complex32::new(real_wavenumber, attenuation_np_per_m);
    Complex32::new(angular_frequency_rad_s, 0.0) / complex_wavenumber
}

fn stable_complex_sqrt(value: Complex32) -> Complex32 {
    // Pick the branch with non-negative imaginary part (decaying/physical branch).
    let mut root = value.sqrt();
    if root.im < 0.0 || (root.im == 0.0 && root.re < 0.0) {
        root = -root;
    }
    root
}

