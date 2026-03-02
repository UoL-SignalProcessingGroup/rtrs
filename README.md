# rtrs

Underwater acoustic ray tracing in Rust.

## Current Features
- 3D ray / beam tracing
- Surface (vacuum) and bottom (rigid, lossy acoustic, lossy elastic) reflections
- 3D sound speed profiles
- 2D bathymetry
- Efficient wideband beam tracing using Gaussian beams ()
- Python bindings via pyo3
- IO with json files
- Receiver representation for grids and arrays
- Parallel processing over the number of rays / beams with rayon
- Optional no-ray-path mode to reduce memory and output size
- Euler and RK2 integration methods
- Input validation

## Planned Features
- Low frequency correction
- Performance improvements

## Usage

### Rust

Build with cargo:
```bash
cargo build --release
```

Run with a JSON input file:
```bash
cargo run --release <path_to_input_file>.json
```

### Python Bindings

Set up a virtual environment and install dependencies:
```bash
python -m venv --prompt rtrs .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

Build the Rust extension into the active environment:
```bash
maturin develop --release
```

Rebuild after any changes to the Rust code:
```bash
maturin develop --release
```

> **Note:** If both `VIRTUAL_ENV` and `CONDA_PREFIX` are set (e.g. conda base is active), maturin will error. Run `unset CONDA_PREFIX` before `maturin develop`, or `conda deactivate` first.

See `examples/` for usage, specifically `examples/munk_test_pyo3.py` for Python bindings.

### PyPI

When published to PyPI, the package can be installed with:
```bash
pip install rtrs
```

Pre-built wheels are distributed per platform and Python version via GitHub Actions using `maturin-action`, so no Rust toolchain is required.

## License
MIT License. See LICENSE file for details.

