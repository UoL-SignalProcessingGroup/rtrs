use anyhow::Result;
use serde_json::{json, Value};
use std::fs;

use crate::input::config::SimulationConfig;
use crate::influence::PressureField;

pub fn write_json(file_path: &str, simulation_config: &SimulationConfig, ray_paths: Vec<Vec<[f32; 3]>>, pressure_field: PressureField) -> Result<()> {

    // Build ray_paths object: { "ray_0": [[x,y,z], ...], ... }
    let mut ray_paths_obj = serde_json::Map::new();
    for (i, path) in ray_paths.iter().enumerate() {
        let key = format!("ray_{}", i);
        let pts: Vec<Value> = path.iter().map(|pt| json!([pt[0], pt[1], pt[2]])).collect();
        ray_paths_obj.insert(key, Value::Array(pts));
    }

    // Flatten pressure field real/imag parts from Array4<Complex32> (nfreq, nx, ny, nz)
    let shape = pressure_field.pressure.dim();
    let (nfreq, nx, ny, nz) = (shape.0, shape.1, shape.2, shape.3);
    let mut re_flat: Vec<f32> = Vec::with_capacity(nfreq * nx * ny * nz);
    let mut im_flat: Vec<f32> = Vec::with_capacity(nfreq * nx * ny * nz);
    for ifreq in 0..nfreq {
        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    let v = pressure_field.pressure[(ifreq, ix, iy, iz)];
                    re_flat.push(v.re);
                    im_flat.push(v.im);
                }
            }
        }
    }

    // Flatten delay_s and amplitude from Array3<f32> (nx, ny, nz)
    let dshape = pressure_field.delay_s.dim();
    let (dnx, dny, dnz) = (dshape.0, dshape.1, dshape.2);
    let mut delay_flat: Vec<f32> = Vec::with_capacity(dnx * dny * dnz);
    let mut amp_flat: Vec<f32> = Vec::with_capacity(dnx * dny * dnz);
    for ix in 0..dnx {
        for iy in 0..dny {
            for iz in 0..dnz {
                delay_flat.push(pressure_field.delay_s[(ix, iy, iz)]);
                amp_flat.push(pressure_field.amplitude[(ix, iy, iz)]);
            }
        }
    }

    // Build pressure_field section
    let mut pf_obj = serde_json::Map::new();
    pf_obj.insert("frequency_hz".into(), json!(simulation_config.source.freq_hz));
    if pressure_field.is_array {
        let recs = pressure_field.receiver_positions_m.as_ref()
            .expect("receiver_positions_m present for array mode");
        let recs_val: Vec<Value> = recs.iter().map(|r| json!([r[0], r[1], r[2]])).collect();
        pf_obj.insert("receiver_positions_m".into(), Value::Array(recs_val));
    } else {
        pf_obj.insert("x_m".into(), json!(pressure_field.x_m));
        pf_obj.insert("y_m".into(), json!(pressure_field.y_m));
        pf_obj.insert("z_m".into(), json!(pressure_field.z_m));
    }
    pf_obj.insert("delay_s".into(), json!({ "shape": [dnx, dny, dnz], "data": delay_flat }));
    pf_obj.insert("amplitude".into(), json!({ "shape": [dnx, dny, dnz], "data": amp_flat }));
    pf_obj.insert("pressure_re".into(), json!({ "shape": [nfreq, nx, ny, nz], "data": re_flat }));
    pf_obj.insert("pressure_im".into(), json!({ "shape": [nfreq, nx, ny, nz], "data": im_flat }));

    let output = json!({
        "src": {
            "frequency_hz": simulation_config.source.freq_hz,
            "source_position_m": simulation_config.source.position,
            "launch_elev_deg": simulation_config.source.launch_elev_deg,
            "launch_azim_deg": simulation_config.source.launch_azim_deg
        },
        "ray_paths": Value::Object(ray_paths_obj),
        "bty": {
            "x_bty_m": simulation_config.bathymetry.x_bty_m,
            "y_bty_m": simulation_config.bathymetry.y_bty_m,
            "z_bty_m": simulation_config.bathymetry.z_bty_m
        },
        "pressure_field": Value::Object(pf_obj)
    });

    fs::write(file_path, serde_json::to_string(&output)?)?;
    Ok(())
}
