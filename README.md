# rtrs

Underwater acoustic ray tracing in Rust.

## Current Features
- 3D ray / beam tracing
- Surface (vacuum) and bottom (rigid) reflections
- 3D sound speed profiles
- 2D bathymetry
- Efficient wideband beam tracing using Gaussian beams
- Python bindings via pyo3
- IO with json and netCDF files respectively
- Receiver representation for grids and arrays

## Planned Features
- Parallel procssing
- More boundary conditions (fluid-fluid and fluid-elastic)
- Absorption / losses (boundary and volume)
- Alternate ray time based time domain formulation

## Usage

Build with cargo:

```bash
cargo build --release
```

Build with python bindings:
```
maturin develop --release --features python
```

Run with json input file:
```bash
cargo run --release <path_to_input_file>.json
```

see `examples/` for more usage details, specifically `examples/munk_test_pyo3.py` for python bindings usage.

