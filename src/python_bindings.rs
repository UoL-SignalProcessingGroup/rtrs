#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::PyAny;
#[cfg(feature = "python")]
use pyo3::exceptions::PyRuntimeError;

#[cfg(feature = "python")]
use crate::engine;
#[cfg(feature = "python")]
use crate::input::config::SimulationConfig;

#[cfg(feature = "python")]
#[pyfunction]
fn run_simulation(py: Python, py_cfg: &PyAny) -> PyResult<PyObject> {
    let json_mod = py.import("json").map_err(|e| PyRuntimeError::new_err(format!("failed to import json: {}", e)))?;
    let json_str: String = json_mod.call_method1("dumps", (py_cfg,))
        .map_err(|e| PyRuntimeError::new_err(format!("json.dumps failed: {}", e)))?
        .extract()
        .map_err(|e| PyRuntimeError::new_err(format!("extract json string failed: {}", e)))?;

    let cfg: SimulationConfig = serde_json::from_str(&json_str)
        .map_err(|e| PyRuntimeError::new_err(format!("failed to parse config JSON: {}", e)))?;

    let (ray_paths, pressure_field) = engine::core(&cfg);

    #[derive(serde::Serialize)]
    struct PressureOut<'a> {
    x_m: &'a Vec<f32>,
    y_m: &'a Vec<f32>,
    z_m: &'a Vec<f32>,
    // for array-mode receivers, explicit positions (N x 3) will be returned instead
    receiver_positions: Option<&'a Vec<[f32; 3]>>,
        frequency_hz: &'a Vec<f32>,
        shape: (usize, usize, usize, usize),
        pressure_re: Vec<f32>,
        pressure_im: Vec<f32>,
        delay_s: Vec<f32>,
        amplitude: Vec<f32>,
    }

    #[derive(serde::Serialize)]
    struct Out<'a> {
        ray_paths: &'a Vec<Vec<[f32; 3]>>,
        pressure_field: PressureOut<'a>,
    }

    let shape = pressure_field.pressure.dim();
    let (nfreq, nx, ny, nz) = (shape.0, shape.1, shape.2, shape.3);
    let mut re_flat: Vec<f32> = Vec::with_capacity(nfreq * nx * ny * nz);
    let mut im_flat: Vec<f32> = Vec::with_capacity(nfreq * nx * ny * nz);
    // flatten delay/amplitude 3D arrays
    let mut delay_flat: Vec<f32> = Vec::with_capacity(nx * ny * nz);
    let mut amp_flat: Vec<f32> = Vec::with_capacity(nx * ny * nz);
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

    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                delay_flat.push(pressure_field.delay_s[(ix, iy, iz)]);
                amp_flat.push(pressure_field.amplitude[(ix, iy, iz)]);
            }
        }
    }

    let p_out = PressureOut{
        x_m: &pressure_field.x_m,
        y_m: &pressure_field.y_m,
        z_m: &pressure_field.z_m,
        receiver_positions: pressure_field.receiver_positions.as_ref(),
        frequency_hz: &cfg.source.freq_hz,
        shape: (nfreq, nx, ny, nz),
        pressure_re: re_flat,
        pressure_im: im_flat,
        delay_s: delay_flat,
        amplitude: amp_flat,
    };

    let out = Out { ray_paths: &ray_paths, pressure_field: p_out };

    let out_json = serde_json::to_string(&out).map_err(|e| PyRuntimeError::new_err(format!("failed to serialize output: {}", e)))?;
    let py_obj = json_mod.call_method1("loads", (out_json,))
        .map_err(|e| PyRuntimeError::new_err(format!("json.loads failed: {}", e)))?;

    Ok(py_obj.into())
}

#[cfg(feature = "python")]
#[pymodule]
fn rtrs(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(run_simulation, m)?)?;
    Ok(())
}
