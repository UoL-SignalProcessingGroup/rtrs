#![allow(dead_code)]

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::PyAny;
#[cfg(feature = "python")]
use pyo3::exceptions::PyRuntimeError;

// declare crate modules at library root so both the CLI (main.rs) and the python
// extension see the same module tree
mod input;
mod output;
mod rays;
mod ssp;
mod bty;
mod reflect;
mod influence;
mod utils;
mod engine;

#[cfg(feature = "python")]
mod python_bindings;

// library root intentionally minimal; python bindings are feature-gated in
// `python_bindings.rs` so normal cargo builds (without `--features python`)
// don't pull in pyo3/linker requirements.
