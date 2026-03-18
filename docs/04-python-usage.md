# Python Usage

This page shows how to run `rtrs` from Python.

There are two common paths:

- PyO3 bindings (`import rtrs`) for direct in-memory calls.
- Subprocess CLI (`cargo run ...`) for file-based workflows.

## 1. Setup

Create and activate an environment:

```bash
python -m venv --prompt rtrs .venv
source .venv/bin/activate
```

Install dependencies and the package:

```bash
pip install -e ".[dev]"
maturin develop --release
```

If you only need plotting helpers from examples:

```bash
pip install -e ".[viz]"
```

## 2. Build Input Dictionary

Input structure is documented in `docs/02-input-reference.md`.

Important:

- `ssp.c_m_s` must be flattened in C-order with `z` fastest, then `y`, then `x`.
- `bathymetry.z_bty_m` must be flattened in C-order with `y` fastest, then `x`.

Minimal pattern:

```python
import numpy as np

z = np.linspace(0.0, 150.0, 2)
c_grid = np.array([
    [[1500.0, 1503.0], [1500.0, 1503.0]],
    [[1500.0, 1503.0], [1500.0, 1503.0]],
], dtype=float)  # shape (nx, ny, nz)

env = {
    "ssp": {
        "x_ssp_m": [0.0, 30000.0],
        "y_ssp_m": [0.0, 30000.0],
        "z_ssp_m": z.tolist(),
        "c_m_s": c_grid.flatten(order="C").tolist(),
    },
    "bathymetry": {
        "x_bty_m": [0.0, 30000.0],
        "y_bty_m": [0.0, 30000.0],
        "z_bty_m": np.array([[150.0, 150.0], [150.0, 150.0]], dtype=float).flatten(order="C").tolist(),
        "bottom_model": {"model": "rigid"},
    },
    "source": {
        "position": [0.0, 0.0, 50.0],
        "freq_hz": [100.0],
        "launch_elev_deg": np.linspace(-10.0, 10.0, 21).tolist(),
        "launch_azim_deg": [0.0],
    },
    "receivers": {
        "config_type": "grid",
        "x_rcvr_m": [0.0],
        "y_rcvr_m": np.linspace(0.0, 30000.0, 200).tolist(),
        "z_rcvr_m": np.linspace(10.0, 140.0, 130).tolist(),
    },
    "beam": {
        "step_m": 10.0,
        "max_steps": 40000,
        "max_range_m": 30000.0,
        "store_ray_paths": False,
        "integration_method": "euler",
    },
}
```

## 3. Run with PyO3 (`import rtrs`)

```python
import rtrs

result = rtrs.run_simulation(env)
```

Return object keys:

- `result["pressure_field"]` always present.
- `result["ray_paths"]` present only when `beam.store_ray_paths = true`.

### Reconstruct Complex Pressure

```python
import numpy as np

pf = result["pressure_field"]
shape = tuple(pf["shape"])  # (nfreq, nx, ny, nz)

re = np.asarray(pf["pressure_re"], dtype=np.float32).reshape(shape)
im = np.asarray(pf["pressure_im"], dtype=np.float32).reshape(shape)
pressure = re + 1j * im

# Optional outputs
amplitude = np.asarray(pf["amplitude"], dtype=np.float32).reshape(shape[1:])
delay_s = np.asarray(pf["delay_s"], dtype=np.float32).reshape(shape[1:])
```

### Grid vs Array Receivers

Grid mode (`receivers.config_type = "grid"`):

```python
x_m = np.asarray(pf["x_m"], dtype=float)
y_m = np.asarray(pf["y_m"], dtype=float)
z_m = np.asarray(pf["z_m"], dtype=float)
```

Array mode (`receivers.config_type = "array"`):

```python
receiver_positions_m = np.asarray(pf["receiver_positions_m"], dtype=float)
# shape is still (nfreq, nx, ny, nz), with nx*ny*nz equal to number of receivers
```

## 4. Run via Subprocess (CLI JSON Files)

This is useful when you want explicit input/output files.

```python
import json
import subprocess
from pathlib import Path

in_path = Path("examples/my_case.json")
out_path = in_path.with_suffix(".out.json")

in_path.write_text(json.dumps(env, indent=2))
subprocess.run(["cargo", "run", "--release", str(in_path)], check=True)

result_json = json.loads(out_path.read_text())
print(result_json.keys())
```

Output structure is documented in `docs/03-output-reference.md`.

## 5. Common Pitfalls

- Typoed keys are ignored at runtime with warnings.
  - Example: `"intergration_method"` is ignored; use `"integration_method"`.
- Flattening order mistakes can silently produce wrong physics.
- In array mode, receiver coordinates are returned under `pressure_field.receiver_positions_m` instead of `x_m/y_m/z_m`.
- After Rust changes affecting bindings, re-run:

```bash
maturin develop --release
```

## 6. Related Docs

- Install/build: `docs/01-install-and-build.md`
- Input dictionary: `docs/02-input-reference.md`
- Output format: `docs/03-output-reference.md`
