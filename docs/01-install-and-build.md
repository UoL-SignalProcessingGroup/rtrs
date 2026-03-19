# Install and Build

This project is primarily a standalone Rust binary that can be accessed through a python interface or that reads a JSON config and writes a JSON output file.

## Prerequisites

- Rust toolchain (`cargo`, `rustc`)
- Python 3.8+ only if you want to run Python examples or bindings

## Build and Run (Rust Binary)

Build release binary:

```bash
cargo build --release
```

Run with a config file:

```bash
cargo run --release examples/testlinear.json
```

The output path is derived from the input path by replacing the extension with `.out.json`.

Example:

- Input: `examples/testlinear.json`
- Output: `examples/testlinear.out.json`

## Python Environment (Optional, but recommended)

If you run example scripts that call the binary and/or use Python plotting helpers:

```bash
python -m venv --prompt rtrs .venv
source .venv/bin/activate
pip install -e ".[viz]"
```

## Python Bindings (Optional, but recommended)

For iterative Rust + Python development:

```bash
python -m venv --prompt rtrs .venv
source .venv/bin/activate
pip install -e ".[dev]"
maturin develop --release
```

Re-run `maturin develop --release` after Rust changes that affect bindings.

## Common Commands

Run a specific example script:

```bash
.venv/bin/python examples/munk_test.py
```

Generate Rust docs for local code inspection:

```bash
cargo doc --no-deps --features python --open
```