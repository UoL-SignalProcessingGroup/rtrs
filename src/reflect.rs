use crate::rays::Ray;
use std::f64::consts::PI;

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
    }
}