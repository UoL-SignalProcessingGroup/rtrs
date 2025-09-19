# rtrs

Underwater acoustic ray tracing in Rust.

## Features (current)
- JSON input configuration input (see `examples/env_simple.json`)
- HDF5 output of ray paths
- 3D sound speed profiles with trilinear interpolation
- Flat surface pressure release boundary

## Planned features
- Acoustic Bottom reflections
- 2D batyhmetry grid
- Beam tracing and pressure calculation

## Usage

Build:
```
cargo build
```

Run with example:
```
cargo run --release -- examples/env_simple.json
```

