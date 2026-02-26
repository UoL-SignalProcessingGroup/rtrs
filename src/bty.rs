use crate::input::config::{BottomBoundaryModel as BottomBoundaryInputModel, SimulationConfig};
use ndarray::Array2;

#[derive(Clone)]
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

pub struct BTYfield {
    pub x: Vec<f32>,
    pub y: Vec<f32>,
    pub z: Array2<f32>,
    pub bottom_model: BottomBoundaryRuntimeModel,
    pub water_density_g_cm3: f32,
}

pub fn init_bty(confg: &SimulationConfig) -> BTYfield {

    let (nx, ny) = (
        confg.bathymetry.x_bty_m.len(),
        confg.bathymetry.y_bty_m.len(),
    );

    let z_bty = Array2::from_shape_vec((nx, ny), confg.bathymetry.z_bty_m.clone())
        .expect(&format!(
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
            compressional_attenuation_db_per_wavelength: compressional_attenuation_db_per_wavelength.unwrap_or(0.0),
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
            compressional_attenuation_db_per_wavelength: compressional_attenuation_db_per_wavelength.unwrap_or(0.0),
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

pub fn interpolate_bty(position: [f32; 3], bty_field: &BTYfield) -> f32 {
    // Bilinear interpolation of bottom depth at (x,y)
    bilinear_interpolation(position, &bty_field.z, &bty_field.x, &bty_field.y)
}

fn bilinear_interpolation(
    position: [f32; 3],
    field: &Array2<f32>,
    x: &[f32],
    y: &[f32],
) -> f32 {
    // Find index function: i, j such that x[i] <= position[0] < x[i+1], etc.
    let find_index = |arr: &[f32], val: f32| -> usize {
        match arr.binary_search_by(|probe| probe.partial_cmp(&val).unwrap()) {
            Ok(i) => i.min(arr.len() - 2), // exact match
            Err(i) => i.saturating_sub(1).min(arr.len() - 2), // interval before insert position
        }
    };

    let i = find_index(x, position[0]);
    let j = find_index(y, position[1]);

    // Normalized distances within the cell
    let xd = ((position[0] - x[i]) / (x[i + 1] - x[i])).clamp(0.0, 1.0);
    let yd = ((position[1] - y[j]) / (y[j + 1] - y[j])).clamp(0.0, 1.0);

    // Compact bilinear interpolation
    let mut c = 0.0;
    for dx in 0..=1 {
        for dy in 0..=1 {
            let weight =
                (if dx == 0 { 1.0 - xd } else { xd }) *
                (if dy == 0 { 1.0 - yd } else { yd });

            c += field[[i + dx, j + dy]] * weight;
        }
    }

    c
}
