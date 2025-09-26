use ndarray::array;

use crate::bty;
use crate::input::config::SimulationConfig;
use crate::ssp::{
    interpolate_c, 
    interpolate_grad_c,
    interpolate_partials_c,
    calculate_ray_partials_c,
    SSPFields
};

use crate::reflect::{
    bottom_reflections, surface_reflection
};

#[derive(Clone)]
pub struct Ray {
    pub position: [f64; 3],
    pub direction: [f64; 3],
    pub phi: f64, // ray torsion 
    pub travel_time: f64,
    pub amplitude: f64,
    pub phase: f64,
    pub num_top_bounces: u32,
    pub num_bottom_bounces: u32,
    pub p_tilde: [f64; 2], // paraxial vectors
    pub q_tilde: [f64; 2],
    pub p_hat: [f64; 2],
    pub q_hat: [f64; 2],
    pub det_q: f64,
    pub c: f64, // local sound speed
}


pub fn trace_ray(
    azim: f64, 
    elev: f64, 
    config: &SimulationConfig, 
    ssp_field: &SSPFields,
    bty_field: &bty::BTYfield
) -> Vec<Ray> {
    // ray tracing main loop

    // println!("Tracing ray at elev {:.2} deg, azim {:.2} deg", elev.to_degrees(), azim.to_degrees());

    let max_n_steps = config.beam.max_steps;

    let c = interpolate_c(config.source.position, &ssp_field);


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
    };

    let mut ray_history: Vec<Ray> = Vec::with_capacity(max_n_steps);
    ray_history.push(ray_init);

    for step in 0..max_n_steps {

        // perform Euler step
        euler_step_ray(&mut ray_history, config.beam.step_m, step, ssp_field);

        // check for boundary reflections
        surface_reflection(&mut ray_history);
        bottom_reflections(&mut ray_history, bty_field);

        // check for termination conditions
        if check_max_range(&mut ray_history, config.beam.max_range_m, config.source.position) {
            break;
        }
        
        
    }

    return ray_history;
}


fn euler_step_ray(ray_history: &mut Vec<Ray>, ds: f64, step: usize, ssp: &SSPFields) {
    // Euler Method

    // allocate ray states for step
    let ray0 = ray_history[step].clone();   // start state
    let mut ray1 = ray0.clone();    // step state

    // find local sound speed and gradients
    let c0 = interpolate_c(ray0.position, ssp);
    let grad_c0 = interpolate_grad_c(ray0.position, ssp);
    let partial_c0 = interpolate_partials_c(ray0.position, ssp);

    // update local sound speed
    ray1.c = c0;

    // find ray normals
    let (e1, e2) = ray_normal(ray0.direction, ray0.phi, c0);

    // calculate ray-centered partials of c
    let [cnn, cmn, cmm] = calculate_ray_partials_c(
        partial_c0[0], partial_c0[3], partial_c0[4],
        partial_c0[1], partial_c0[5], partial_c0[2],
        e1, e2);
    
    // unit direction scaled by local c
    let unit_direction = [
        c0 * ray0.direction[0],
        c0 * ray0.direction[1],
        c0 * ray0.direction[2],
    ];

    // first half step position update
    ray1.position[0] += ds * unit_direction[0];
    ray1.position[1] += ds * unit_direction[1];
    ray1.position[2] += ds * unit_direction[2];

    // direction half update
    ray1.direction[0] -= ds * grad_c0[0] / c0.powi(2);
    ray1.direction[1] -= ds * grad_c0[1] / c0.powi(2);
    ray1.direction[2] -= ds * grad_c0[2] / c0.powi(2);

    // phi half update
    ray1.phi += ds * c0.recip() * ray0.direction[2] 
    * (ray0.direction[1] * grad_c0[0] - ray0.direction[0] * grad_c0[1])
    / (ray0.direction[0].powi(2) + ray0.direction[1].powi(2));

    // Construct the 2x2 matrix of ray-centered sound speed second derivatives (c_mat0)
    let c_mat0 = array![
        [ -cnn / c0.powi(2), -cmn / c0.powi(2) ],
        [ -cmn / c0.powi(2), -cmm / c0.powi(2) ],
    ];

    // Update paraxial vectors (p_tilde, q_tilde, p_hat, q_hat) using the c_mat0 matrix and c0 scalar
    // p_tilde = p_tilde + half_ds * c_mat0 * q_tilde
    let q_tilde_vec = ndarray::arr1(&ray0.q_tilde);
    let p_tilde_update = c_mat0.dot(&q_tilde_vec);
    ray1.p_tilde[0] = ray0.p_tilde[0] + ds * p_tilde_update[0];
    ray1.p_tilde[1] = ray0.p_tilde[1] + ds * p_tilde_update[1];

    // q_tilde = q_tilde + half_ds * c0 * p_tilde
    ray1.q_tilde[0] = ray0.q_tilde[0] + ds * c0 * ray0.p_tilde[0];
    ray1.q_tilde[1] = ray0.q_tilde[1] + ds * c0 * ray0.p_tilde[1];

    // p_hat = p_hat + half_ds * c_mat0 * q_hat
    let q_hat_vec = ndarray::arr1(&ray0.q_hat);
    let p_hat_update = c_mat0.dot(&q_hat_vec);
    ray1.p_hat[0] = ray0.p_hat[0] + ds * p_hat_update[0];
    ray1.p_hat[1] = ray0.p_hat[1] + ds * p_hat_update[1];

    // q_hat = q_hat + half_ds * c0 * p_hat
    ray1.q_hat[0] = ray0.q_hat[0] + ds * c0 * ray0.p_hat[0];
    ray1.q_hat[1] = ray0.q_hat[1] + ds * c0 * ray0.p_hat[1];

    // travel time
    ray1.travel_time += ds /c0;
    
    ray1.amplitude = ray0.amplitude;
    ray1.phase = ray0.phase;
    ray1.num_top_bounces = ray0.num_top_bounces;
    ray1.num_bottom_bounces = ray0.num_bottom_bounces;

    // append new ray state
    ray_history.push(ray1);

}


pub fn ray_normal(direction: [f64; 3], phi: f64, c: f64) -> ([f64; 3], [f64; 3]) {
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


fn check_max_range(ray_history: &mut Vec<Ray>, max_range: f64, source_position: [f64; 3]) -> bool {
    let ray = ray_history.last().unwrap();
    let dx = ray.position[0] - source_position[0];
    let dy = ray.position[1] - source_position[1];
    let range = (dx.powi(2) + dy.powi(2)).sqrt();
    range >= max_range
}