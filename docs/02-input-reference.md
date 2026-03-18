# Input Reference

`rtrs` reads a JSON configuration matching the `SimulationConfig` structure in `src/input.rs`.

For editor tooling and machine validation, see `docs/schema/simulation-config.schema.json`.

## Top-Level Structure

Required top-level keys:

- `ssp`
- `bathymetry`
- `source`
- `receivers`
- `beam`

Units use SI unless noted.
Depth uses positive-down convention (`z >= 0` expected).

## `ssp`

- `x_ssp_m`: number array, required
- `y_ssp_m`: number array, required
- `z_ssp_m`: number array, required
- `c_m_s`: number array, required

Meaning:

- `x_ssp_m`, `y_ssp_m`, `z_ssp_m` define axis coordinates for the 3D SSP grid.
- `c_m_s` contains flattened sound-speed values in m/s.

Validation behavior:

- `c_m_s` must be non-empty.
- All `c_m_s` values must be `> 0`.
- Negative `z_ssp_m` values are corrected to absolute value with a warning.

Flattening convention for `c_m_s` (important):

- Grid shape is `(nx, ny, nz)` where:
  - `nx = len(x_ssp_m)`
  - `ny = len(y_ssp_m)`
  - `nz = len(z_ssp_m)`
- `c_m_s` length must be `nx * ny * nz`.
- Values must be flattened with `z` fastest, then `y`, then `x` (equivalent to nested loops `for x` -> `for y` -> `for z`).

Example index mapping:

- `flat_index = ((ix * ny) + iy) * nz + iz`
- `c_m_s[flat_index] = c[ix, iy, iz]`

## `bathymetry`

- `x_bty_m`: number array, required
- `y_bty_m`: number array, required
- `z_bty_m`: number array, required
- `water_density_g_cm3`: number, optional
- `bottom_model`: object, optional (default: rigid)

Validation behavior:

- `z_bty_m` must be non-empty.
- `water_density_g_cm3` must be `> 0` when provided.

Flattening convention for `z_bty_m`:

- Grid shape is `(nx, ny)` where:
  - `nx = len(x_bty_m)`
  - `ny = len(y_bty_m)`
- `z_bty_m` length should be `nx * ny`.
- Values should be flattened with `y` fastest, then `x` (equivalent to nested loops `for x` -> `for y`).

Example index mapping:

- `flat_index = ix * ny + iy`
- `z_bty_m[flat_index] = z[ix, iy]`

### `bathymetry.bottom_model`

Tagged union with discriminator key `model`:

- `{"model": "rigid"}`
- `{"model": "acoustic", ...}`
- `{"model": "elastic", ...}`

`acoustic` fields:

- `compressional_speed_m_s`: number, required, `> 0`
- `density_g_cm3`: number, required, `> 0`
- `compressional_attenuation_db_per_wavelength`: number, optional, `>= 0`

`elastic` fields:

- `compressional_speed_m_s`: number, required, `> 0`
- `shear_speed_m_s`: number, required, `> 0`
- `density_g_cm3`: number, required, `> 0`
- `compressional_attenuation_db_per_wavelength`: number, optional, `>= 0`
- `shear_attenuation_db_per_wavelength`: number, optional, `>= 0`

## `source`

- `position`: number array of length 3 `[x, y, z]` in meters, required
- `freq_hz`: number array, required
- `launch_elev_deg`: number array, required
- `launch_azim_deg`: number array, required

Validation behavior:

- `launch_elev_deg` values must be within `[-90, 90]`.
- Negative source depth (`position[2]`) is corrected to absolute value with a warning.

## `receivers`

- `config_type`: string, required: `"grid"` or `"array"`
- `x_rcvr_m`: number array, required
- `y_rcvr_m`: number array, required
- `z_rcvr_m`: number array, required

Meaning:

- If `config_type = "grid"`, `x_rcvr_m`, `y_rcvr_m`, `z_rcvr_m` are axis vectors.
- If `config_type = "array"`, vectors represent explicit receiver coordinates.

Validation behavior:

- `config_type` must be `"grid"` or `"array"`.
- Negative receiver depths are corrected to absolute value with a warning.

## `beam`

- `step_m`: number, required, `> 0`
- `max_steps`: integer, required, `> 0`
- `max_range_m`: number, required, `> 0`
- `store_ray_paths`: boolean, optional, default `false`
- `show_progress`: boolean, optional, default `false`
- `atomic_progress_counter`: boolean, optional, default `false`
- `integration_method`: string, optional, `"euler"` (default) or `"rk2"`

## Unknown Keys

Runtime behavior currently ignores unknown keys and emits warnings.

Example warning form:

- `WARNING unknown input key is ignored: config.<path>`

Schema note:

- The schema in `docs/schema/simulation-config.schema.json` is intentionally strict (`additionalProperties: false`) to catch typos early.

## Minimal Example

```json
{
  "ssp": {
    "x_ssp_m": [0.0, 30000.0],
    "y_ssp_m": [0.0, 30000.0],
    "z_ssp_m": [0.0, 150.0],
    "c_m_s": [1500.0, 1503.0, 1500.0, 1503.0, 1500.0, 1503.0, 1500.0, 1503.0]
  },
  "bathymetry": {
    "x_bty_m": [0.0, 30000.0],
    "y_bty_m": [0.0, 30000.0],
    "z_bty_m": [150.0, 150.0, 150.0, 150.0],
    "bottom_model": {
      "model": "elastic",
      "compressional_speed_m_s": 1700.0,
      "shear_speed_m_s": 400.0,
      "density_g_cm3": 1.6,
      "compressional_attenuation_db_per_wavelength": 0.5,
      "shear_attenuation_db_per_wavelength": 0.35
    }
  },
  "source": {
    "position": [0.0, 0.0, 50.0],
    "freq_hz": [100.0],
    "launch_elev_deg": [-10.0, 0.0, 10.0],
    "launch_azim_deg": [0.0]
  },
  "receivers": {
    "config_type": "grid",
    "x_rcvr_m": [0.0, 1000.0],
    "y_rcvr_m": [0.0],
    "z_rcvr_m": [20.0, 40.0]
  },
  "beam": {
    "step_m": 5.0,
    "max_steps": 20000,
    "max_range_m": 30000.0,
    "store_ray_paths": false,
    "show_progress": false,
    "atomic_progress_counter": false,
    "integration_method": "euler"
  }
}
```
