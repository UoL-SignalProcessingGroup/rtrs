#![allow(unsafe_op_in_unsafe_fn)]
//! Python bindings for running rtrs simulations from Python.

#[cfg(feature = "python")]
use pyo3::exceptions::PyRuntimeError;
#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::PyAny;
#[cfg(feature = "python")]
use pyo3::types::PyDict;

#[cfg(feature = "python")]
use crate::engine;
#[cfg(feature = "python")]
use crate::input::config::SimulationConfig;

#[cfg(feature = "python")]
#[pyfunction(text_signature = "(py_cfg, /)")]
/// Run an rtrs simulation from a Python configuration object.
///
/// Parameters
/// ----------
/// py_cfg : dict-like
///     Simulation configuration following the rtrs input schema.
///
/// Returns
/// -------
/// dict
///     Output dictionary containing `pressure_field` and optional `ray_paths`.
///
/// Notes
/// -----
/// Unknown config keys are ignored and emitted as Python warnings.
pub fn run_simulation(py: Python, py_cfg: &PyAny) -> PyResult<PyObject> {
    let json_mod = py
        .import("json")
        .map_err(|e| PyRuntimeError::new_err(format!("failed to import json: {}", e)))?;
    let json_str: String = json_mod
        .call_method1("dumps", (py_cfg,))
        .map_err(|e| PyRuntimeError::new_err(format!("json.dumps failed: {}", e)))?
        .extract()
        .map_err(|e| PyRuntimeError::new_err(format!("extract json string failed: {}", e)))?;

    let mut unknown_fields = Vec::new();
    let mut de = serde_json::Deserializer::from_str(&json_str);
    let mut cfg: SimulationConfig = serde_ignored::deserialize(&mut de, |path| {
        unknown_fields.push(path.to_string());
    })
    .map_err(|e| PyRuntimeError::new_err(format!("failed to parse config JSON: {}", e)))?;

    let mut warnings = cfg
        .validate()
        .map_err(|e| PyRuntimeError::new_err(format!("invalid simulation config: {}", e)))?;

    warnings.extend(
        unknown_fields
            .into_iter()
            .map(|field| format!("unknown input key is ignored: config.{}", field)),
    );

    // Surface validation warnings as Python warnings so they are visible regardless of environment
    if !warnings.is_empty() {
        let warnings_mod = py
            .import("warnings")
            .map_err(|e| PyRuntimeError::new_err(format!("failed to import warnings: {}", e)))?;
        for w in &warnings {
            warnings_mod
                .call_method1("warn", (format!("rtrs: {}", w),))
                .map_err(|e| PyRuntimeError::new_err(format!("warnings.warn failed: {}", e)))?;
        }
    }

    let (ray_paths, pressure_field) = engine::core(&cfg);

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

    let out = PyDict::new(py);
    if let Some(ray_paths) = ray_paths.as_ref() {
        out.set_item("ray_paths", ray_paths)
            .map_err(|e| PyRuntimeError::new_err(format!("failed to set ray_paths: {}", e)))?;
    }

    let p_out = PyDict::new(py);
    if pressure_field.is_array {
        let receiver_positions = pressure_field
            .receiver_positions_m
            .as_ref()
            .expect("receiver_positions_m present for array mode");
        p_out
            .set_item("receiver_positions_m", receiver_positions)
            .map_err(|e| {
                PyRuntimeError::new_err(format!("failed to set receiver_positions_m: {}", e))
            })?;
    } else {
        p_out
            .set_item("x_m", &pressure_field.x_m)
            .map_err(|e| PyRuntimeError::new_err(format!("failed to set x_m: {}", e)))?;
        p_out
            .set_item("y_m", &pressure_field.y_m)
            .map_err(|e| PyRuntimeError::new_err(format!("failed to set y_m: {}", e)))?;
        p_out
            .set_item("z_m", &pressure_field.z_m)
            .map_err(|e| PyRuntimeError::new_err(format!("failed to set z_m: {}", e)))?;
    }

    p_out
        .set_item("frequency_hz", &cfg.source.freq_hz)
        .map_err(|e| PyRuntimeError::new_err(format!("failed to set frequency_hz: {}", e)))?;
    p_out
        .set_item("shape", (nfreq, nx, ny, nz))
        .map_err(|e| PyRuntimeError::new_err(format!("failed to set shape: {}", e)))?;
    p_out
        .set_item("pressure_re", re_flat)
        .map_err(|e| PyRuntimeError::new_err(format!("failed to set pressure_re: {}", e)))?;
    p_out
        .set_item("pressure_im", im_flat)
        .map_err(|e| PyRuntimeError::new_err(format!("failed to set pressure_im: {}", e)))?;
    p_out
        .set_item("delay_s", delay_flat)
        .map_err(|e| PyRuntimeError::new_err(format!("failed to set delay_s: {}", e)))?;
    p_out
        .set_item("amplitude", amp_flat)
        .map_err(|e| PyRuntimeError::new_err(format!("failed to set amplitude: {}", e)))?;

    out.set_item("pressure_field", p_out)
        .map_err(|e| PyRuntimeError::new_err(format!("failed to set pressure_field: {}", e)))?;

    Ok(out.into())
}

#[cfg(feature = "python")]
/// Python extension module exposing rtrs simulation functions.
#[pymodule]
fn rtrs(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add(
        "__doc__",
        "Python bindings for the rtrs underwater acoustic ray tracer.",
    )?;
    m.add_function(wrap_pyfunction!(run_simulation, m)?)?;
    Ok(())
}
