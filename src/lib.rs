#![allow(dead_code)]

// declare crate modules at library root so both the CLI (main.rs) and the python
// extension see the same module tree
mod bty;
mod engine;
mod influence;
mod input;
mod rays;
mod reflect;
mod ssp;
mod utils;

#[cfg(feature = "python")]
mod python_bindings;

// library root intentionally minimal; python bindings are feature-gated in
// `python_bindings.rs` so normal cargo builds (without `--features python`)
// don't pull in pyo3/linker requirements.
