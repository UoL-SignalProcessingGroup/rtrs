use crate::bty;
use crate::input::config::{IntegrationMethod, SimulationConfig};
use crate::ssp::{
    init_ssp_cursor,
    interpolate_c_with_cursor,
    interpolate_grad_c_with_cursor,
    interpolate_partials_c_with_cursor,
    reduce_step_to_ssp_interfaces,
    update_ssp_cursor,
    calculate_ray_partials_c,
    SSPCursor,
    SSPFields
};

use crate::reflect::{
    reflect_boundaries,
};

#[derive(Clone, Copy)]
pub struct BottomBounceMetadata {
    pub boundary_normal: [f32; 3],
    pub incident_slowness: [f32; 3],
    pub water_sound_speed_m_s: f32,
}

#[derive(Clone)]
pub struct Ray {
    pub position: [f32; 3],
    pub direction: [f32; 3],
    pub phi: f32, // ray torsion 
    pub travel_time: f32,
    pub amplitude: f32,
    pub phase: f32,
    pub num_top_bounces: u32,
    pub num_bottom_bounces: u32,
    pub p_tilde: [f32; 2], // paraxial vectors
    pub q_tilde: [f32; 2],
    pub p_hat: [f32; 2],
    pub q_hat: [f32; 2],
    pub det_q: f32,
    pub c: f32, // local sound speed
    pub bottom_bounce: Option<BottomBounceMetadata>,
}

struct RayDerivatives {
    d_position: [f32; 3],
    d_direction: [f32; 3],
    d_phi: f32,
    d_p_tilde: [f32; 2],
    d_q_tilde: [f32; 2],
    d_p_hat: [f32; 2],
    d_q_hat: [f32; 2],
    c: f32,
}


pub fn trace_ray(
    azim: f32, 
    elev: f32, 
    config: &SimulationConfig, 
    ssp_field: &SSPFields,
    bty_field: &bty::BTYfield
) -> Vec<Ray> {
    // ray tracing main loop

    let max_n_steps = config.beam.max_steps;

    let mut ssp_cursor = init_ssp_cursor(config.source.position, ssp_field);
    let c = interpolate_c_with_cursor(config.source.position, ssp_field, &mut ssp_cursor);
    let mut bty_cursor = bty::init_bty_cursor(config.source.position, bty_field);

    let ray_init = Ray {
        position: config.source.position,
        direction: [ elev.cos() * azim.sin() / c, elev.cos() * azim.cos() / c, elev.sin()/ c],
        phi: 0.0,
        travel_time: 0.0,
        amplitude: 1.0,
        phase: 0.0,
        num_top_bounces: 0,
        num_bottom_bounces: 0,
        p_tilde: [1.0, 0.0],
        q_tilde: [0.0, 0.0],
        p_hat: [0.0, 1.0],
        q_hat: [0.0, 0.0],
        det_q: 0.0,
        c: c,
        bottom_bounce: None,
    };

    let mut ray_history: Vec<Ray> = Vec::with_capacity(max_n_steps);
    ray_history.push(ray_init);

    for step in 0..max_n_steps {

        match config.beam.integration_method {
            IntegrationMethod::Euler => {
                euler_step_ray(
                    &mut ray_history,
                    config.beam.step_m,
                    step,
                    ssp_field,
                    &mut ssp_cursor,
                    bty_field,
                    &mut bty_cursor,
                );
            }
            IntegrationMethod::Rk2 => {
                midpoint_rk2_step_ray(
                    &mut ray_history,
                    config.beam.step_m,
                    step,
                    ssp_field,
                    &mut ssp_cursor,
                    bty_field,
                    &mut bty_cursor,
                );
            }
        }

        // check for boundary reflections
        reflect_boundaries(&mut ray_history, bty_field, &mut bty_cursor);

        // check for max range termination conditions
        if check_max_range(&mut ray_history, config.beam.max_range_m, config.source.position) {
            break;
        }
        
        // check for number of bottom bounces
        // if ray_history.last().unwrap().num_bottom_bounces >= 10 {
        //     break;
        // }
    }

    return ray_history;
}


fn reduce_step_to_boundaries(
    position: [f32; 3],
    unit_direction: [f32; 3],
    step: f32,
    ssp: &SSPFields,
    ssp_cursor: &SSPCursor,
    bty_field: &bty::BTYfield,
    bty_cursor: &mut bty::BTYCursor,
) -> f32 {
    let eps = 1.0e-8_f32;
    let mut h = step;

    h = reduce_step_to_ssp_interfaces(position, unit_direction, h, ssp, ssp_cursor);
    h = bty::reduce_step_to_bty_segments(position, unit_direction, h, bty_field, bty_cursor);

    if unit_direction[2] < -eps {
        let h_top = (0.0 - position[2]) / unit_direction[2];
        if h_top > eps && h_top < h {
            h = h_top;
        }
    }

    let z_bty0 = bty::interpolate_bty_with_cursor(position, bty_field, bty_cursor);
    let (normal, _) = bty::bottom_normal_at_with_cursor(position, bty_field, bty_cursor);

    let boundary_point = [position[0], position[1], z_bty0];
    let n_dot_u = normal[0] * unit_direction[0] + normal[1] * unit_direction[1] + normal[2] * unit_direction[2];
    if n_dot_u.abs() > eps {
        let n_dot_s_minus_p =
            normal[0] * (boundary_point[0] - position[0]) +
            normal[1] * (boundary_point[1] - position[1]) +
            normal[2] * (boundary_point[2] - position[2]);
        let h_bot = n_dot_s_minus_p / n_dot_u;
        if h_bot > eps && h_bot < h {
            h = h_bot;
        }
    }

    let min_step = 1.0e-6_f32 * step;
    if h < min_step { min_step } else { h }
}


fn euler_step_ray(
    ray_history: &mut Vec<Ray>,
    ds: f32,
    step: usize,
    ssp: &SSPFields,
    ssp_cursor: &mut SSPCursor,
    bty_field: &bty::BTYfield,
    bty_cursor: &mut bty::BTYCursor,
) {
    let ray0 = &ray_history[step];
    let derivatives_at_start = compute_ray_derivatives(ray0, ssp, ssp_cursor);

    let unit_direction = derivatives_at_start.d_position;

    let h = reduce_step_to_boundaries(
        ray0.position,
        unit_direction,
        ds,
        ssp,
        ssp_cursor,
        bty_field,
        bty_cursor,
    );

    // first half step position update
    let position = [
        ray0.position[0] + h * unit_direction[0],
        ray0.position[1] + h * unit_direction[1],
        ray0.position[2] + h * unit_direction[2],
    ];

    // direction update
    let direction = [
        ray0.direction[0] + h * derivatives_at_start.d_direction[0],
        ray0.direction[1] + h * derivatives_at_start.d_direction[1],
        ray0.direction[2] + h * derivatives_at_start.d_direction[2],
    ];

    let phi = ray0.phi + h * derivatives_at_start.d_phi;

    let p_tilde = [
        ray0.p_tilde[0] + h * derivatives_at_start.d_p_tilde[0],
        ray0.p_tilde[1] + h * derivatives_at_start.d_p_tilde[1],
    ];

    let q_tilde = [
        ray0.q_tilde[0] + h * derivatives_at_start.d_q_tilde[0],
        ray0.q_tilde[1] + h * derivatives_at_start.d_q_tilde[1],
    ];

    let p_hat = [
        ray0.p_hat[0] + h * derivatives_at_start.d_p_hat[0],
        ray0.p_hat[1] + h * derivatives_at_start.d_p_hat[1],
    ];

    let q_hat = [
        ray0.q_hat[0] + h * derivatives_at_start.d_q_hat[0],
        ray0.q_hat[1] + h * derivatives_at_start.d_q_hat[1],
    ];

    let ray1 = Ray {
        position,
        direction,
        phi,
        travel_time: ray0.travel_time + h / derivatives_at_start.c,
        amplitude: ray0.amplitude,
        phase: ray0.phase,
        num_top_bounces: ray0.num_top_bounces,
        num_bottom_bounces: ray0.num_bottom_bounces,
        p_tilde,
        q_tilde,
        p_hat,
        q_hat,
        det_q: ray0.det_q,
        c: derivatives_at_start.c,
        bottom_bounce: None,
    };

    // append new ray state
    ray_history.push(ray1);

    let last_pos = ray_history.last().unwrap().position;
    update_ssp_cursor(last_pos, ssp, ssp_cursor);
    bty::update_bty_cursor(last_pos, bty_field, bty_cursor);

}

fn midpoint_rk2_step_ray(
    ray_history: &mut Vec<Ray>,
    ds: f32,
    step: usize,
    ssp: &SSPFields,
    ssp_cursor: &mut SSPCursor,
    bty_field: &bty::BTYfield,
    bty_cursor: &mut bty::BTYCursor,
) {
    let ray0 = &ray_history[step];
    let derivatives_at_start = compute_ray_derivatives(ray0, ssp, ssp_cursor);

    let h = reduce_step_to_boundaries(
        ray0.position,
        derivatives_at_start.d_position,
        ds,
        ssp,
        ssp_cursor,
        bty_field,
        bty_cursor,
    );

    let half_h = 0.5 * h;
    let midpoint_ray = Ray {
        position: [
            ray0.position[0] + half_h * derivatives_at_start.d_position[0],
            ray0.position[1] + half_h * derivatives_at_start.d_position[1],
            ray0.position[2] + half_h * derivatives_at_start.d_position[2],
        ],
        direction: [
            ray0.direction[0] + half_h * derivatives_at_start.d_direction[0],
            ray0.direction[1] + half_h * derivatives_at_start.d_direction[1],
            ray0.direction[2] + half_h * derivatives_at_start.d_direction[2],
        ],
        phi: ray0.phi + half_h * derivatives_at_start.d_phi,
        travel_time: ray0.travel_time + half_h / derivatives_at_start.c,
        amplitude: ray0.amplitude,
        phase: ray0.phase,
        num_top_bounces: ray0.num_top_bounces,
        num_bottom_bounces: ray0.num_bottom_bounces,
        p_tilde: [
            ray0.p_tilde[0] + half_h * derivatives_at_start.d_p_tilde[0],
            ray0.p_tilde[1] + half_h * derivatives_at_start.d_p_tilde[1],
        ],
        q_tilde: [
            ray0.q_tilde[0] + half_h * derivatives_at_start.d_q_tilde[0],
            ray0.q_tilde[1] + half_h * derivatives_at_start.d_q_tilde[1],
        ],
        p_hat: [
            ray0.p_hat[0] + half_h * derivatives_at_start.d_p_hat[0],
            ray0.p_hat[1] + half_h * derivatives_at_start.d_p_hat[1],
        ],
        q_hat: [
            ray0.q_hat[0] + half_h * derivatives_at_start.d_q_hat[0],
            ray0.q_hat[1] + half_h * derivatives_at_start.d_q_hat[1],
        ],
        det_q: ray0.det_q,
        c: derivatives_at_start.c,
        bottom_bounce: None,
    };

    let mut midpoint_cursor = *ssp_cursor;
    update_ssp_cursor(midpoint_ray.position, ssp, &mut midpoint_cursor);
    let derivatives_at_midpoint = compute_ray_derivatives(&midpoint_ray, ssp, &mut midpoint_cursor);

    let ray1 = Ray {
        position: [
            ray0.position[0] + h * derivatives_at_midpoint.d_position[0],
            ray0.position[1] + h * derivatives_at_midpoint.d_position[1],
            ray0.position[2] + h * derivatives_at_midpoint.d_position[2],
        ],
        direction: [
            ray0.direction[0] + h * derivatives_at_midpoint.d_direction[0],
            ray0.direction[1] + h * derivatives_at_midpoint.d_direction[1],
            ray0.direction[2] + h * derivatives_at_midpoint.d_direction[2],
        ],
        phi: ray0.phi + h * derivatives_at_midpoint.d_phi,
        travel_time: ray0.travel_time + h / derivatives_at_midpoint.c,
        amplitude: ray0.amplitude,
        phase: ray0.phase,
        num_top_bounces: ray0.num_top_bounces,
        num_bottom_bounces: ray0.num_bottom_bounces,
        p_tilde: [
            ray0.p_tilde[0] + h * derivatives_at_midpoint.d_p_tilde[0],
            ray0.p_tilde[1] + h * derivatives_at_midpoint.d_p_tilde[1],
        ],
        q_tilde: [
            ray0.q_tilde[0] + h * derivatives_at_midpoint.d_q_tilde[0],
            ray0.q_tilde[1] + h * derivatives_at_midpoint.d_q_tilde[1],
        ],
        p_hat: [
            ray0.p_hat[0] + h * derivatives_at_midpoint.d_p_hat[0],
            ray0.p_hat[1] + h * derivatives_at_midpoint.d_p_hat[1],
        ],
        q_hat: [
            ray0.q_hat[0] + h * derivatives_at_midpoint.d_q_hat[0],
            ray0.q_hat[1] + h * derivatives_at_midpoint.d_q_hat[1],
        ],
        det_q: ray0.det_q,
        c: derivatives_at_midpoint.c,
        bottom_bounce: None,
    };

    ray_history.push(ray1);

    let last_pos = ray_history.last().unwrap().position;
    update_ssp_cursor(last_pos, ssp, ssp_cursor);
    bty::update_bty_cursor(last_pos, bty_field, bty_cursor);
}

fn compute_ray_derivatives(
    ray: &Ray,
    ssp: &SSPFields,
    ssp_cursor: &mut SSPCursor,
) -> RayDerivatives {
    let c = interpolate_c_with_cursor(ray.position, ssp, ssp_cursor);
    let grad_c = interpolate_grad_c_with_cursor(ray.position, ssp, ssp_cursor);
    let partial_c = interpolate_partials_c_with_cursor(ray.position, ssp, ssp_cursor);

    let (e1, e2) = ray_normal(ray.direction, ray.phi, c);
    let [cnn, cmn, cmm] = calculate_ray_partials_c(
        partial_c[0], partial_c[3], partial_c[4],
        partial_c[1], partial_c[5], partial_c[2],
        e1, e2,
    );

    let d_position = [
        c * ray.direction[0],
        c * ray.direction[1],
        c * ray.direction[2],
    ];

    let c_pow2 = c * c;
    let d_direction = [
        -grad_c[0] / c_pow2,
        -grad_c[1] / c_pow2,
        -grad_c[2] / c_pow2,
    ];

    let denom = ray.direction[0].powi(2) + ray.direction[1].powi(2);
    let d_phi = if denom > 1.0e-12 {
        c.recip() * ray.direction[2]
            * (ray.direction[1] * grad_c[0] - ray.direction[0] * grad_c[1])
            / denom
    } else {
        0.0
    };

    let c11 = -cnn / c_pow2;
    let c12 = -cmn / c_pow2;
    let c22 = -cmm / c_pow2;

    let d_p_tilde = [
        c11 * ray.q_tilde[0] + c12 * ray.q_tilde[1],
        c12 * ray.q_tilde[0] + c22 * ray.q_tilde[1],
    ];
    let d_q_tilde = [
        c * ray.p_tilde[0],
        c * ray.p_tilde[1],
    ];
    let d_p_hat = [
        c11 * ray.q_hat[0] + c12 * ray.q_hat[1],
        c12 * ray.q_hat[0] + c22 * ray.q_hat[1],
    ];
    let d_q_hat = [
        c * ray.p_hat[0],
        c * ray.p_hat[1],
    ];

    RayDerivatives {
        d_position,
        d_direction,
        d_phi,
        d_p_tilde,
        d_q_tilde,
        d_p_hat,
        d_q_hat,
        c,
    }
}


pub fn ray_normal(direction: [f32; 3], phi: f32, c: f32) -> ([f32; 3], [f32; 3]) {
    // compute the ray normal vector e1, e2

    let mut e1 = [0.0; 3];
    let mut e2 = [0.0; 3];

    let rl = (direction[0].powi(2) + direction[1].powi(2)).sqrt();

    if phi != 0.0 {
        // e1
        e1[0] = (c * direction[0] * direction[2] * phi.cos() + direction[1] * phi.sin()) / rl;
        e1[1] = (c * direction[1] * direction[2] * phi.cos() - direction[0] * phi.sin()) / rl;
        e1[2] = -c * rl * phi.cos();

        // e2
        e2[0] = (c * direction[0] * direction[2] * phi.sin() - direction[1] * phi.cos()) / rl;
        e2[1] = (c * direction[1] * direction[2] * phi.sin() + direction[0] * phi.cos()) / rl;
        e2[2] = -c * rl * phi.sin();
    } else {
        // e1
        e1[0] = c * direction[0] * direction[2] / rl;
        e1[1] = c * direction[1] * direction[2] / rl;
        e1[2] = -c * rl;

        // e2
        e2[0] = -direction[1] / rl;
        e2[1] = direction[0] / rl;
        e2[2] = 0.0;
    }
    return (e1, e2);
}


fn check_max_range(ray_history: &mut Vec<Ray>, max_range: f32, source_position: [f32; 3]) -> bool {
    let ray = ray_history.last().unwrap();
    let dx = ray.position[0] - source_position[0];
    let dy = ray.position[1] - source_position[1];
    let range = (dx.powi(2) + dy.powi(2)).sqrt();
    range >= max_range
}