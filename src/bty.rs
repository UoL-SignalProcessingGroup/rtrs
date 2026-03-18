use crate::input::config::{BottomBoundaryModel as BottomBoundaryInputModel, SimulationConfig};
use ndarray::Array2;

#[derive(Clone, Copy, Debug)]
/// Cursor that tracks the current bathymetry cell index.
pub struct BTYCursor {
    pub i: usize,
    pub j: usize,
}

#[derive(Clone)]
/// Runtime bottom model with defaults resolved from input.
pub enum BottomBoundaryRuntimeModel {
    Rigid,
    Acoustic {
        compressional_speed_m_s: f32,
        density_g_cm3: f32,
        compressional_attenuation_db_per_wavelength: f32,
    },
    Elastic {
        compressional_speed_m_s: f32,
        shear_speed_m_s: f32,
        density_g_cm3: f32,
        compressional_attenuation_db_per_wavelength: f32,
        shear_attenuation_db_per_wavelength: f32,
    },
}

/// Bathymetry grid and runtime bottom parameters.
pub struct BTYfield {
    pub x: Vec<f32>,
    pub y: Vec<f32>,
    pub z: Array2<f32>,
    pub bottom_model: BottomBoundaryRuntimeModel,
    pub water_density_g_cm3: f32,
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

/// Initialize a bathymetry cursor at a world position.
pub fn init_bty_cursor(position: [f32; 3], bty: &BTYfield) -> BTYCursor {
    BTYCursor {
        i: find_cell_index(&bty.x, position[0]),
        j: find_cell_index(&bty.y, position[1]),
    }
}

pub fn update_bty_cursor(position: [f32; 3], bty: &BTYfield, cursor: &mut BTYCursor) {
    march_cell_index(&bty.x, position[0], &mut cursor.i);
    march_cell_index(&bty.y, position[1], &mut cursor.j);
}

pub fn reduce_step_to_bty_segments(
    position: [f32; 3],
    unit_direction: [f32; 3],
    step: f32,
    bty: &BTYfield,
    cursor: &BTYCursor,
) -> f32 {
    let eps = 1.0e-9_f32;
    let mut h = step;

    if unit_direction[0].abs() > eps {
        let boundary_x = if unit_direction[0] > 0.0 {
            bty.x[cursor.i + 1]
        } else {
            bty.x[cursor.i]
        };
        let h_cross_x = (boundary_x - position[0]) / unit_direction[0];
        if h_cross_x > eps && h_cross_x < h {
            h = h_cross_x;
        }
    }

    if unit_direction[1].abs() > eps {
        let boundary_y = if unit_direction[1] > 0.0 {
            bty.y[cursor.j + 1]
        } else {
            bty.y[cursor.j]
        };
        let h_cross_y = (boundary_y - position[1]) / unit_direction[1];
        if h_cross_y > eps && h_cross_y < h {
            h = h_cross_y;
        }
    }

    h
}

/// Build bathymetry runtime fields from validated simulation input.
pub fn init_bty(confg: &SimulationConfig) -> BTYfield {
    let (nx, ny) = (
        confg.bathymetry.x_bty_m.len(),
        confg.bathymetry.y_bty_m.len(),
    );

    let z_bty =
        Array2::from_shape_vec((nx, ny), confg.bathymetry.z_bty_m.clone()).expect(&format!(
            "z_bty_m does not match grid dimensions: nx = {}, ny = {}, bty_len = {}",
            nx,
            ny,
            confg.bathymetry.z_bty_m.len()
        ));

    let bottom_model = match &confg.bathymetry.bottom_model {
        BottomBoundaryInputModel::Rigid => BottomBoundaryRuntimeModel::Rigid,
        BottomBoundaryInputModel::Acoustic {
            compressional_speed_m_s,
            density_g_cm3,
            compressional_attenuation_db_per_wavelength,
        } => BottomBoundaryRuntimeModel::Acoustic {
            compressional_speed_m_s: *compressional_speed_m_s,
            density_g_cm3: *density_g_cm3,
            compressional_attenuation_db_per_wavelength:
                compressional_attenuation_db_per_wavelength.unwrap_or(0.0),
        },
        BottomBoundaryInputModel::Elastic {
            compressional_speed_m_s,
            shear_speed_m_s,
            density_g_cm3,
            compressional_attenuation_db_per_wavelength,
            shear_attenuation_db_per_wavelength,
        } => BottomBoundaryRuntimeModel::Elastic {
            compressional_speed_m_s: *compressional_speed_m_s,
            shear_speed_m_s: *shear_speed_m_s,
            density_g_cm3: *density_g_cm3,
            compressional_attenuation_db_per_wavelength:
                compressional_attenuation_db_per_wavelength.unwrap_or(0.0),
            shear_attenuation_db_per_wavelength: shear_attenuation_db_per_wavelength.unwrap_or(0.0),
        },
    };

    let bty_field = BTYfield {
        x: confg.bathymetry.x_bty_m.clone(),
        y: confg.bathymetry.y_bty_m.clone(),
        z: z_bty,
        bottom_model,
        water_density_g_cm3: confg.bathymetry.water_density_g_cm3.unwrap_or(1.0),
    };

    return bty_field;
}

pub fn interpolate_bty_from_cursor(
    position: [f32; 3],
    bty_field: &BTYfield,
    cursor: &BTYCursor,
) -> f32 {
    bilinear_interpolation_with_indices(
        position,
        &bty_field.z,
        &bty_field.x,
        &bty_field.y,
        cursor.i,
        cursor.j,
    )
}

pub fn bottom_normal_from_cursor(
    position: [f32; 3],
    bty_field: &BTYfield,
    cursor: &BTYCursor,
) -> ([f32; 3], [f32; 3]) {
    let i = cursor.i;
    let j = cursor.j;

    let dx = bty_field.x[i + 1] - bty_field.x[i];
    let dy = bty_field.y[j + 1] - bty_field.y[j];

    let xd = ((position[0] - bty_field.x[i]) / dx).clamp(0.0, 1.0);
    let yd = ((position[1] - bty_field.y[j]) / dy).clamp(0.0, 1.0);

    let z00 = bty_field.z[[i, j]];
    let z10 = bty_field.z[[i + 1, j]];
    let z01 = bty_field.z[[i, j + 1]];
    let z11 = bty_field.z[[i + 1, j + 1]];

    let dzdx_y0 = (z10 - z00) / dx;
    let dzdx_y1 = (z11 - z01) / dx;
    let dzdx = (1.0 - yd) * dzdx_y0 + yd * dzdx_y1;

    let dzdy_x0 = (z01 - z00) / dy;
    let dzdy_x1 = (z11 - z10) / dy;
    let dzdy = (1.0 - xd) * dzdy_x0 + xd * dzdy_x1;

    let nx = -dzdx;
    let ny = -dzdy;
    let nz = 1.0;
    let norm_len = (nx * nx + ny * ny + nz * nz).sqrt();
    let normal = [nx / norm_len, ny / norm_len, nz / norm_len];

    let tangent = [-ny, nx, 0.0];
    let tnorm = (tangent[0] * tangent[0] + tangent[1] * tangent[1]).sqrt();
    let tangent = if tnorm > 0.0 {
        [tangent[0] / tnorm, tangent[1] / tnorm, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    };

    (normal, tangent)
}

fn bilinear_interpolation_with_indices(
    position: [f32; 3],
    field: &Array2<f32>,
    x: &[f32],
    y: &[f32],
    i: usize,
    j: usize,
) -> f32 {
    let xd = ((position[0] - x[i]) / (x[i + 1] - x[i])).clamp(0.0, 1.0);
    let yd = ((position[1] - y[j]) / (y[j + 1] - y[j])).clamp(0.0, 1.0);

    let mut c = 0.0;
    for dx in 0..=1 {
        for dy in 0..=1 {
            let weight =
                (if dx == 0 { 1.0 - xd } else { xd }) * (if dy == 0 { 1.0 - yd } else { yd });

            c += field[[i + dx, j + dy]] * weight;
        }
    }

    c
}
