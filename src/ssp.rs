use crate::input::config::SimulationConfig;
use ndarray::Array3;

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

    let dx = if x.len() >= 2 {
        x[1] - x[0]
    } else {
        1.0
    };
    let dy = if y.len() >= 2 {
        y[1] - y[0]
    } else {
        1.0
    };
    let dz = if z.len() >= 2 {
        z[1] - z[0]
    } else {
        1.0
    };

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
        c: c.clone(),
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


pub fn interpolate_c(position: [f32; 3], ssp: &SSPFields) -> f32 {
    let c = trilinear_interpolation(position, &ssp.c, &ssp.x, &ssp.y, &ssp.z);
    return c;
}

pub fn interpolate_grad_c(position: [f32; 3], ssp: &SSPFields) -> [f32; 3] {
    // interpolate grad c at position
    let cx = trilinear_interpolation(position, &ssp.cx, &ssp.x, &ssp.y, &ssp.z);
    let cy = trilinear_interpolation(position, &ssp.cy, &ssp.x, &ssp.y, &ssp.z);
    let cz = trilinear_interpolation(position, &ssp.cz, &ssp.x, &ssp.y, &ssp.z);
    return [cx, cy, cz];
}

pub fn interpolate_partials_c(position: [f32; 3], ssp: &SSPFields) -> [f32; 6] {
    // interpolate cxx, cyy, czz, cxy, cxz, cyz at position
    let cxx = trilinear_interpolation(position, &ssp.cxx, &ssp.x, &ssp.y, &ssp.z);
    let cyy = trilinear_interpolation(position, &ssp.cyy, &ssp.x, &ssp.y, &ssp.z);
    let czz = trilinear_interpolation(position, &ssp.czz, &ssp.x, &ssp.y, &ssp.z);
    let cxy = trilinear_interpolation(position, &ssp.cxy, &ssp.x, &ssp.y, &ssp.z);
    let cxz = trilinear_interpolation(position, &ssp.cxz, &ssp.x, &ssp.y, &ssp.z);
    let cyz = trilinear_interpolation(position, &ssp.cyz, &ssp.x, &ssp.y, &ssp.z);
    return [cxx, cyy, czz, cxy, cxz, cyz];
}

pub fn calculate_ray_partials_c(
    cxx: f32, cxy: f32, cxz: f32, cyy: f32, cyz: f32, czz: f32, e1: [f32; 3], e2: [f32; 3]) -> [f32; 3] {
    // calculate cnn cmn cmm curvature of sound speed (ray-centered)
    let cnn = cxx * e1[0].powi(2) + cyy * e1[1].powi(2) + czz * e1[2].powi(2)
        + 2.0 * cxy * e1[0] * e1[1]
        + 2.0 * cxz * e1[0] * e1[2]
        + 2.0 * cyz * e1[1] * e1[2];

    let cmn = cxx * e1[0] * e2[0] + cyy * e1[1] * e2[1] + czz * e1[2] * e2[2]
        + cxy * (e1[0] * e2[1] + e2[0] * e1[1])
        + cxz * (e1[0] * e2[2] + e2[0] * e1[2])
        + cyz * (e1[1] * e2[2] + e2[1] * e1[2]);

    let cmm = cxx * e2[0].powi(2) + cyy * e2[1].powi(2) + czz * e2[2].powi(2)
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


fn trilinear_interpolation(
    position: [f32; 3], 
    field: &Array3<f32>, 
    x: &[f32], 
    y: &[f32], 
    z: &[f32]
) -> f32 {
    // Trilinear interpolation of field at position [x, y, z]

    // Find indices function  i, j, k such that x[i] <= position[0] < x[i+1], etc.
    // use binary search as x,y,z arrays are sorted
    let find_index = |arr: &[f32], val: f32| -> usize {
        match arr.binary_search_by(|probe| probe.partial_cmp(&val).unwrap()) {
            Ok(i) => i.min(arr.len() - 2), // exact match
            Err(i) => i.saturating_sub(1).min(arr.len() - 2), // interval before insert position
        }
    };

    // find indiceis
    let i = find_index(x, position[0]);
    let j = find_index(y, position[1]);
    let k = find_index(z, position[2]);

    // Compute normalized distances
    let xd = ((position[0] - x[i]) / (x[i + 1] - x[i])).clamp(0.0, 1.0);
    let yd = ((position[1] - y[j]) / (y[j + 1] - y[j])).clamp(0.0, 1.0);
    let zd = ((position[2] - z[k]) / (z[k + 1] - z[k])).clamp(0.0, 1.0);

    // Compact trilinear interpolation
    let mut c = 0.0;
    for dx in 0..=1 {
        for dy in 0..=1 {
            for dz in 0..=1 {
                let weight = 
                    (if dx == 0 { 1.0 - xd } else { xd }) *
                    (if dy == 0 { 1.0 - yd } else { yd }) *
                    (if dz == 0 { 1.0 - zd } else { zd });

                c += field[[i+dx, j+dy, k+dz]] * weight;
            }
        }
    }

    c
}


