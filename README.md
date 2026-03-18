# rtrs

3D Underwater acoustic ray tracing with arbitrary Cartesian receiver coordinates, geometry independent broadband beam tracing, and Python bindings.

## Current Features
- Exclusively 3D ray / beam tracing
- Surface (vacuum) and bottom (rigid, lossy acoustic, lossy elastic) reflections
- 3D sound speed profiles ( $c(x,y,z)$ )
- 2D bathymetry ( $b(x,y)$ )
- Efficient wideband beam tracing using Gaussian beams (ray geometry is traced as frequency independent, then the beam influence is calculated at receiver coordinates for each frequency)
- Python bindings via pyo3
- IO with json files
- Receiver representation for grids and arrays
- Parallel processing over the number of rays / beams with rayon
- Optional no-ray-path mode to reduce memory and output size
- Euler and RK2 integration methods
- Input validation

## Installing, Building, and Running

There are 3 main ways to use rtrs: directly with Rust and Cargo, with Python bindings, or using CLI with the pre-built binary and JSON input files (Not covered here). The first step is to clone the repository and navigate to the project directory:
```bash
git clone https://github.com/fincb/rtrs.git
cd rtrs
```
Rust, Cargo, and are required for all methods, Python is also required for the Python bindings. The recommended way (for linux and macos) to install Rust and Cargo is with `rustup`:
```bash
curl https://sh.rustup.rs -sSf | sh
```

Then check:
```bash
rustc --version
cargo --version
rustup --version
```

### Building the Docs

The documentation is built with Cargo. To build the docs locally, run:
```bash
cargo doc --no-deps --features python
```

Or to open the docs in the browser after building:
```bash
cargo doc --no-deps --features python --open
```

The generated docs include user-facing guide pages under `rtrs::guides`:

- `rtrs::guides::install_and_build`
- `rtrs::guides::input_reference`
- `rtrs::guides::output_reference`
- `rtrs::guides::python_usage`


### Rust & Cargo Build and Run

Build with cargo (the resulting binary will be in `target/release/`):
```bash
cargo build --release
```

Run with a JSON input file (compiles and runs in release mode):
```bash
cargo run --release <path_to_input_file>.json
```

### Python Bindings Normal Build and Run

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

Then in a Python script the package can be imported and used:
```python
import rtrs

env = # <Dictionary with environment parameters>
result = rtrs.run_simulation(env)
```

### Python Bindings with plotting for examples

For running the examples with plotting, install the `viz` extra dependencies:
```bash
python -m venv --prompt rtrs .venv
source .venv/bin/activate
pip install -e ".[viz]"
```

or with conda:
```bash
conda activate <env_name>
pip install -e ".[viz]"
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

For plotting and development, install both extras:
```bash
pip install -e ".[dev,viz]"
```

Build the Rust extension into the active environment and to rebuild after any changes to the Rust code:
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

## Possible Future Features
- Low frequency correction
- Performance improvements
- Multiple sources and Source directivity patterns

## Documentation
- Primary docs via `cargo doc` (guides and API together)
- JSON schema: `docs/schema/simulation-config.schema.json`


## License
MIT License. See LICENSE file for details.

