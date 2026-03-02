# rtrs

3D Underwater acoustic ray tracing in Rust.

## Current Features
- 3D ray / beam tracing
- Surface (vacuum) and bottom (rigid, lossy acoustic, lossy elastic) reflections
- 3D sound speed profiles ($c(x,y,z)$)
- 2D bathymetry ($b(x,y)$)
- Efficient wideband beam tracing using Gaussian beams (ray geometry is traced as frequency independent, then the beam influence is calculated at receiver coordinates for each frequency)
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

Build with cargo (the resulting binary will be in `target/release/`):
```bash
cargo build --release
```

Run with a JSON input file (compiles and runs in release mode):
```bash
cargo run --release <path_to_input_file>.json
```

### Python Bindings Normal Build

Set up a virtual environment and install dependencies:
```bash
python -m venv --prompt rtrs .venv
source .venv/bin/activate
pip install -e .
```

or with conda:
```
conda activate <env_name>
pip install -e .
```

### Python Bindings with Development Mode

Useful for testing changes to the Rust code without needing to reinstall the package after every change. The `dev` extra also includes `maturin` as a dependency for building the Rust extension.

Set up a virtual environment and install dependencies:
```bash
python -m venv --prompt rtrs .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

or with conda:
```bash
conda activate <env_name>
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

See `examples/` for usage with and without python bindings. Note that in the examples which do not use the python bindings, the program is compiled and run in release mode with cargo from within a python script. 

### Future PyPI Release

When published to PyPI, the package can be installed with:
```bash
pip install rtrs
```

Hopefully, pre-built wheels are distributed per platform and Python version via GitHub Actions using `maturin-action`, so no Rust toolchain is required.

## License
MIT License. See LICENSE file for details.

