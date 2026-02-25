use crate::input::config::SimulationConfig;
use ndarray::Array2;

pub struct BTYfield {
    pub x: Vec<f32>,
    pub y: Vec<f32>,
    pub z: Array2<f32>,
    pub bottom_p_wave_speed: Option<f32>, // Bottom halfspace P-wave speed (m/s). `None` gives rigid fallback.
    pub bottom_density: Option<f32>,      // Bottom halfspace density (g/cm3). Only meaningful when `bottom_p_wave_speed` is `Some`.
    pub water_density: f32,               // Water (upper halfspace) density (g/cm3). Defaults to 1.0.
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

    let bty_field = BTYfield {
        x: confg.bathymetry.x_bty_m.clone(),
        y: confg.bathymetry.y_bty_m.clone(),
        z: z_bty,
        bottom_p_wave_speed: confg.bathymetry.bottom_p_wave_speed_m_s,
        bottom_density: confg.bathymetry.bottom_density_g_cm3,
        water_density: confg.bathymetry.water_density_g_cm3.unwrap_or(1.0),
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
