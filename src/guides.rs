//! User-facing documentation guides mirrored from `docs/` into `cargo doc`.
//!
//! This keeps workflow docs visible from Rust API docs and allows linking from
//! type/function pages directly to guide content.

/// Installation, build, and execution workflows.
#[doc = include_str!("../docs/01-install-and-build.md")]
pub mod install_and_build {}

/// JSON input schema and validation behavior.
#[doc = include_str!("../docs/02-input-reference.md")]
pub mod input_reference {}

/// JSON output structure and flattening conventions.
#[doc = include_str!("../docs/03-output-reference.md")]
pub mod output_reference {}

/// Python usage patterns for bindings and CLI workflows.
#[doc = include_str!("../docs/04-python-usage.md")]
pub mod python_usage {}
