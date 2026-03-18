use crate::input::config::SimulationConfig;
use ndarray::Array3;

#[derive(Clone, Copy, Debug)]
/// Cursor that tracks the current sound-speed cell index.
pub struct SSPCursor {
    pub i: usize,
    pub j: usize,
    pub k: usize,
}

/// Sound-speed field and precomputed first/second derivatives.
pub struct SSPFields {
    x: Vec<f32>,
    y: Vec<f32>,
    z: Vec<f32>,
    c: Array3<f32>,
    cx: Array3<f32>,
    cy: Array3<f32>,
    cz: Array3<f32>,
    cxx: Array3<f32>,
    cyy: Array3<f32>,
    czz: Array3<f32>,
    cxy: Array3<f32>,
    cxz: Array3<f32>,
    cyz: Array3<f32>,
}

fn find_cell_index(arr: &[f32], val: f32) -> usize {
    match arr.binary_search_by(|probe| probe.partial_cmp(&val).unwrap()) {
        Ok(i) => i.min(arr.len() - 2),
        Err(i) => i.saturating_sub(1).min(arr.len() - 2),
    }
}

fn march_cell_index(arr: &[f32], val: f32, idx: &mut usize) {
    if arr.len() < 2 {
        *idx = 0;
        return;
    }

    let max_idx = arr.len() - 2;

    if val <= arr[0] {
        *idx = 0;
        return;
    }
    if val >= arr[arr.len() - 1] {
        *idx = max_idx;
        return;
    }

    while *idx < max_idx && val > arr[*idx + 1] {
        *idx += 1;
    }
    while *idx > 0 && val < arr[*idx] {
        *idx -= 1;
    }
}

/// Initialize an SSP cursor at a world position.
pub fn init_ssp_cursor(position: [f32; 3], ssp: &SSPFields) -> SSPCursor {
    SSPCursor {
        i: find_cell_index(&ssp.x, position[0]),
        j: find_cell_index(&ssp.y, position[1]),
        k: find_cell_index(&ssp.z, position[2]),
    }
}

pub fn update_ssp_cursor(position: [f32; 3], ssp: &SSPFields, cursor: &mut SSPCursor) {
    march_cell_index(&ssp.x, position[0], &mut cursor.i);
    march_cell_index(&ssp.y, position[1], &mut cursor.j);
    march_cell_index(&ssp.z, position[2], &mut cursor.k);
}

pub fn reduce_step_to_ssp_interfaces(
    position: [f32; 3],
    unit_direction: [f32; 3],
    step: f32,
    ssp: &SSPFields,
    cursor: &SSPCursor,
) -> f32 {
    let eps = 1.0e-9_f32;
    let mut h = step;

    let try_reduce = |arr: &[f32], pos: f32, vel: f32, idx: usize, h_cur: &mut f32| {
        if vel.abs() <= eps {
            return;
        }

        let boundary = if vel > 0.0 { arr[idx + 1] } else { arr[idx] };

        let h_cross = (boundary - pos) / vel;
        if h_cross > eps && h_cross < *h_cur {
            *h_cur = h_cross;
        }
    };

    try_reduce(&ssp.x, position[0], unit_direction[0], cursor.i, &mut h);
    try_reduce(&ssp.y, position[1], unit_direction[1], cursor.j, &mut h);
    try_reduce(&ssp.z, position[2], unit_direction[2], cursor.k, &mut h);

    h
}

/// Build SSP fields and derivatives from validated simulation input.
pub fn init_ssp(config: &SimulationConfig) -> SSPFields {
    let (nx, ny, nz) = (
        config.ssp.x_ssp_m.len(),
        config.ssp.y_ssp_m.len(),
        config.ssp.z_ssp_m.len(),
    );
    let c = Array3::from_shape_vec((nx, ny, nz), config.ssp.c_m_s.clone())
        .expect("c_m_s does not match grid dimensions");
    // println!("c: {}", c);

    let x = &config.ssp.x_ssp_m;
    let y = &config.ssp.y_ssp_m;
    let z = &config.ssp.z_ssp_m;

    let dx = if x.len() >= 2 { x[1] - x[0] } else { 1.0 };
    let dy = if y.len() >= 2 { y[1] - y[0] } else { 1.0 };
    let dz = if z.len() >= 2 { z[1] - z[0] } else { 1.0 };

    let cx = partial_x(&c, dx);
    let cy = partial_y(&c, dy);
    let cz = partial_z(&c, dz);
    let cxx = partial_x(&cx, dx);
    let cyy = partial_y(&cy, dy);
    let czz = partial_z(&cz, dz);
    let cxy = partial_y(&cx, dy);
    let cxz = partial_z(&cx, dz);
    let cyz = partial_z(&cy, dz);

    let ssp_field = SSPFields {
        c,
        x: x.clone(),
        y: y.clone(),
        z: z.clone(),
        cx: cx,
        cy: cy,
        cz: cz,
        cxx: cxx,
        cyy: cyy,
        czz: czz,
        cxy: cxy,
        cxz: cxz,
        cyz: cyz,
    };

    return ssp_field;
}

/// Interpolate local sound speed at a position using the supplied cursor.
pub fn interpolate_c_with_cursor(
    position: [f32; 3],
    ssp: &SSPFields,
    cursor: &mut SSPCursor,
) -> f32 {
    trilinear_interpolation_with_cursor(position, &ssp.c, &ssp.x, &ssp.y, &ssp.z, cursor)
}

/// Interpolate sound speed, gradient, and Hessian components at a position.
pub fn interpolate_all_with_cursor(
    position: [f32; 3],
    ssp: &SSPFields,
    cursor: &mut SSPCursor,
) -> (f32, [f32; 3], [f32; 6]) {
    let (i, j, k, wx, wy, wz) = update_cursor_and_weights(position, &ssp.x, &ssp.y, &ssp.z, cursor);

    let mut c = 0.0_f32;
    let mut grad = [0.0_f32; 3];
    let mut partials = [0.0_f32; 6];

    for dx in 0..=1 {
        for dy in 0..=1 {
            for dz in 0..=1 {
                let weight = wx[dx] * wy[dy] * wz[dz];
                let ix = i + dx;
                let iy = j + dy;
                let iz = k + dz;

                c += ssp.c[[ix, iy, iz]] * weight;
                grad[0] += ssp.cx[[ix, iy, iz]] * weight;
                grad[1] += ssp.cy[[ix, iy, iz]] * weight;
                grad[2] += ssp.cz[[ix, iy, iz]] * weight;

                partials[0] += ssp.cxx[[ix, iy, iz]] * weight;
                partials[1] += ssp.cyy[[ix, iy, iz]] * weight;
                partials[2] += ssp.czz[[ix, iy, iz]] * weight;
                partials[3] += ssp.cxy[[ix, iy, iz]] * weight;
                partials[4] += ssp.cxz[[ix, iy, iz]] * weight;
                partials[5] += ssp.cyz[[ix, iy, iz]] * weight;
            }
        }
    }

    (c, grad, partials)
}

pub fn calculate_ray_partials_c(
    cxx: f32,
    cxy: f32,
    cxz: f32,
    cyy: f32,
    cyz: f32,
    czz: f32,
    e1: [f32; 3],
    e2: [f32; 3],
) -> [f32; 3] {
    // calculate cnn cmn cmm curvature of sound speed (ray-centered)
    let cnn = cxx * e1[0].powi(2)
        + cyy * e1[1].powi(2)
        + czz * e1[2].powi(2)
        + 2.0 * cxy * e1[0] * e1[1]
        + 2.0 * cxz * e1[0] * e1[2]
        + 2.0 * cyz * e1[1] * e1[2];

    let cmn = cxx * e1[0] * e2[0]
        + cyy * e1[1] * e2[1]
        + czz * e1[2] * e2[2]
        + cxy * (e1[0] * e2[1] + e2[0] * e1[1])
        + cxz * (e1[0] * e2[2] + e2[0] * e1[2])
        + cyz * (e1[1] * e2[2] + e2[1] * e1[2]);

    let cmm = cxx * e2[0].powi(2)
        + cyy * e2[1].powi(2)
        + czz * e2[2].powi(2)
        + 2.0 * cxy * e2[0] * e2[1]
        + 2.0 * cxz * e2[0] * e2[2]
        + 2.0 * cyz * e2[1] * e2[2];

    [cnn, cmn, cmm]
}

fn partial_x(c: &Array3<f32>, dx: f32) -> Array3<f32> {
    let (nx, ny, nz) = c.dim();
    assert!(nx >= 2, "Need nx >= 2 for finite differences along x");
    let mut d = Array3::<f32>::zeros((nx, ny, nz));

    // i = 0 (forward)
    for j in 0..ny {
        for k in 0..nz {
            d[[0, j, k]] = (c[[1, j, k]] - c[[0, j, k]]) / dx;
        }
    }

    // interior (central)
    for i in 1..nx - 1 {
        for j in 0..ny {
            for k in 0..nz {
                d[[i, j, k]] = (c[[i + 1, j, k]] - c[[i - 1, j, k]]) / (2.0 * dx);
            }
        }
    }

    // i = nx-1 (backward)
    for j in 0..ny {
        for k in 0..nz {
            d[[nx - 1, j, k]] = (c[[nx - 1, j, k]] - c[[nx - 2, j, k]]) / dx;
        }
    }

    d
}

fn partial_y(c: &Array3<f32>, dy: f32) -> Array3<f32> {
    let (nx, ny, nz) = c.dim();
    assert!(ny >= 2, "Need ny >= 2 for finite differences along y");
    let mut d = Array3::<f32>::zeros((nx, ny, nz));

    // j = 0 (forward)
    for i in 0..nx {
        for k in 0..nz {
            d[[i, 0, k]] = (c[[i, 1, k]] - c[[i, 0, k]]) / dy;
        }
    }

    // interior (central)
    for j in 1..ny - 1 {
        for i in 0..nx {
            for k in 0..nz {
                d[[i, j, k]] = (c[[i, j + 1, k]] - c[[i, j - 1, k]]) / (2.0 * dy);
            }
        }
    }

    // j = ny-1 (backward)
    for i in 0..nx {
        for k in 0..nz {
            d[[i, ny - 1, k]] = (c[[i, ny - 1, k]] - c[[i, ny - 2, k]]) / dy;
        }
    }

    d
}

fn partial_z(c: &Array3<f32>, dz: f32) -> Array3<f32> {
    let (nx, ny, nz) = c.dim();
    assert!(nz >= 2, "Need nz >= 2 for finite differences along z");
    let mut d = Array3::<f32>::zeros((nx, ny, nz));

    // k = 0 (forward)
    for i in 0..nx {
        for j in 0..ny {
            d[[i, j, 0]] = (c[[i, j, 1]] - c[[i, j, 0]]) / dz;
        }
    }

    // interior (central)
    for k in 1..nz - 1 {
        for i in 0..nx {
            for j in 0..ny {
                d[[i, j, k]] = (c[[i, j, k + 1]] - c[[i, j, k - 1]]) / (2.0 * dz);
            }
        }
    }

    // k = nz-1 (backward)
    for i in 0..nx {
        for j in 0..ny {
            d[[i, j, nz - 1]] = (c[[i, j, nz - 1]] - c[[i, j, nz - 2]]) / dz;
        }
    }

    d
}

// trilinear interpolation with cursor
fn trilinear_interpolation_with_cursor(
    position: [f32; 3],
    field: &Array3<f32>,
    x: &[f32],
    y: &[f32],
    z: &[f32],
    cursor: &mut SSPCursor,
) -> f32 {
    let (i, j, k, wx, wy, wz) = update_cursor_and_weights(position, x, y, z, cursor);

    let mut c = 0.0;
    for dx in 0..=1 {
        for dy in 0..=1 {
            for dz in 0..=1 {
                let weight = wx[dx] * wy[dy] * wz[dz];
                c += field[[i + dx, j + dy, k + dz]] * weight;
            }
        }
    }

    c
}

fn update_cursor_and_weights(
    position: [f32; 3],
    x: &[f32],
    y: &[f32],
    z: &[f32],
    cursor: &mut SSPCursor,
) -> (usize, usize, usize, [f32; 2], [f32; 2], [f32; 2]) {
    march_cell_index(x, position[0], &mut cursor.i);
    march_cell_index(y, position[1], &mut cursor.j);
    march_cell_index(z, position[2], &mut cursor.k);

    let i = cursor.i;
    let j = cursor.j;
    let k = cursor.k;

    let xd = ((position[0] - x[i]) / (x[i + 1] - x[i])).clamp(0.0, 1.0);
    let yd = ((position[1] - y[j]) / (y[j + 1] - y[j])).clamp(0.0, 1.0);
    let zd = ((position[2] - z[k]) / (z[k + 1] - z[k])).clamp(0.0, 1.0);

    let wx = [1.0 - xd, xd];
    let wy = [1.0 - yd, yd];
    let wz = [1.0 - zd, zd];

    (i, j, k, wx, wy, wz)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::config::{
        Bathymetry, BeamSettings, BottomBoundaryModel, IntegrationMethod, Receivers, SimulationConfig,
        SoundSpeed, Source,
    };

    fn base_config_with_ssp(x: Vec<f32>, y: Vec<f32>, z: Vec<f32>, c_m_s: Vec<f32>) -> SimulationConfig {
        SimulationConfig {
            ssp: SoundSpeed {
                x_ssp_m: x,
                y_ssp_m: y,
                z_ssp_m: z,
                c_m_s,
            },
            bathymetry: Bathymetry {
                x_bty_m: vec![0.0, 10.0],
                y_bty_m: vec![0.0, 10.0],
                z_bty_m: vec![100.0; 4],
                water_density_g_cm3: Some(1.0),
                bottom_model: BottomBoundaryModel::Rigid,
            },
            source: Source {
                position: [0.0, 0.0, 5.0],
                freq_hz: vec![100.0],
                launch_elev_deg: vec![0.0],
                launch_azim_deg: vec![0.0],
            },
            receivers: Receivers {
                config_type: "grid".to_string(),
                x_rcvr_m: vec![0.0],
                y_rcvr_m: vec![0.0],
                z_rcvr_m: vec![5.0],
            },
            beam: BeamSettings {
                step_m: 1.0,
                max_steps: 10,
                max_range_m: 50.0,
                store_ray_paths: false,
                show_progress: false,
                atomic_progress_counter: false,
                integration_method: IntegrationMethod::Euler,
            },
        }
    }

    #[test]
    fn init_ssp_builds_expected_shape_and_indexing() {
        let x = vec![0.0, 1.0];
        let y = vec![0.0, 1.0];
        let z = vec![0.0, 1.0];
        let mut c_m_s = Vec::new();
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    c_m_s.push(1500.0 + 100.0 * i as f32 + 10.0 * j as f32 + k as f32);
                }
            }
        }
        let cfg = base_config_with_ssp(x, y, z, c_m_s);
        let ssp = init_ssp(&cfg);

        assert_eq!(ssp.c.dim(), (2, 2, 2));
        assert!((ssp.c[[1, 0, 1]] - 1601.0).abs() < 1.0e-6);
    }

    #[test]
    fn interpolate_all_with_cursor_is_stable_on_corners_and_boundaries() {
        let x = vec![0.0, 1.0];
        let y = vec![0.0, 1.0];
        let z = vec![0.0, 1.0];
        let mut c_m_s = Vec::new();
        for &xi in &x {
            for &yi in &y {
                for &zi in &z {
                    c_m_s.push(1500.0 + xi + 2.0 * yi + 3.0 * zi);
                }
            }
        }

        let cfg = base_config_with_ssp(x, y, z, c_m_s);
        let ssp = init_ssp(&cfg);
        let mut cursor = init_ssp_cursor([0.0, 0.0, 0.0], &ssp);

        let (c_corner, grad_corner, _) = interpolate_all_with_cursor([0.0, 0.0, 0.0], &ssp, &mut cursor);
        assert!((c_corner - 1500.0).abs() < 1.0e-5);
        assert!((grad_corner[0] - 1.0).abs() < 1.0e-5);
        assert!((grad_corner[1] - 2.0).abs() < 1.0e-5);
        assert!((grad_corner[2] - 3.0).abs() < 1.0e-5);

        let (c_boundary, grad_boundary, _) =
            interpolate_all_with_cursor([0.5, 1.0, 0.5], &ssp, &mut cursor);
        assert!((c_boundary - 1504.0).abs() < 1.0e-5);
        assert!(grad_boundary.iter().all(|g| g.is_finite()));
    }

    #[test]
    fn reduce_step_to_ssp_interfaces_returns_smallest_positive_crossing() {
        let x = vec![0.0, 10.0];
        let y = vec![0.0, 10.0];
        let z = vec![0.0, 10.0];
        let cfg = base_config_with_ssp(x, y, z, vec![1500.0; 8]);
        let ssp = init_ssp(&cfg);

        let pos = [2.0, 3.0, 4.0];
        let cursor = init_ssp_cursor(pos, &ssp);

        let h_forward = reduce_step_to_ssp_interfaces(pos, [1.0, 0.1, 0.2], 100.0, &ssp, &cursor);
        assert!((h_forward - 8.0).abs() < 1.0e-6);
        assert!(h_forward > 0.0 && h_forward <= 100.0);

        let h_backward =
            reduce_step_to_ssp_interfaces(pos, [-1.0, -0.1, -0.2], 100.0, &ssp, &cursor);
        assert!((h_backward - 2.0).abs() < 1.0e-6);
        assert!(h_backward > 0.0 && h_backward <= 100.0);
    }
}
