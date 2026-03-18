#![allow(dead_code)]
//! 3D Underwater acoustic ray tracing with arbitrary Cartesian receiver 
//! coordinates, geometry independent broadband beam tracing, and Python bindings.
//!
//! The docs here primarily are intended for IO reference not implementation 
//! details and internal APIs. They cover the Python API and JSON input/output formats.
//!
//! # Cargo Docs Focus
//!
//! The `cargo doc` pages are organized around the same user-facing workflow as
//! the Markdown guides:
//!
//! - [Install and Build](crate::guides::install_and_build)
//! - [Input Reference](crate::guides::input_reference)
//! - [Output Reference](crate::guides::output_reference)
//! - [Python Usage](crate::guides::python_usage)
//!
//! API pages that back those guides:
//!
//! - [Input schema types](crate::input::config)
//! - [Output writer](crate::output::write_json)
//! - [Python binding entrypoint](crate::python_bindings::run_simulation)


// declare crate modules at library root so both the CLI (main.rs) and the python
// extension see the same module tree
/// Bathymetry field construction and interpolation utilities.
mod bty;

/// Top-level simulation orchestration.
mod engine;

/// Gaussian beam influence accumulation on receivers.
mod influence;

/// Input schema and validation.
pub mod input;

/// User-facing guide pages rendered in `cargo doc`.
pub mod guides;

/// JSON output serialization.
pub mod output;

/// Ray state and numerical integration.
mod rays;

/// Surface and bottom reflection models.
mod reflect;

/// Sound-speed field interpolation and derivatives.
mod ssp;

/// Shared math helpers.
mod utils;

#[cfg(feature = "python")]
/// Python bindings exposed when the `python` feature is enabled.
pub mod python_bindings;

// library root intentionally minimal; python bindings are feature-gated in
// `python_bindings.rs` so normal cargo builds (without `--features python`)
// don't pull in pyo3/linker requirements.
