# Output Reference

The standalone binary writes one JSON file per run.

Output path rule:

- Input `foo.json` -> Output `foo.out.json`

Root structure is written in `src/output.rs` and contains:

- `src` (useful for debug metadata and post-processing, not available in Python bindings)
- `ray_paths` (only when enabled)
- `bty`(useful for debug metadata and post-processing, not available in Python bindings)
- `pressure_field`

## `src`

Metadata copied from input source settings:

- `frequency_hz`: number array
- `source_position_m`: `[x, y, z]`
- `launch_elev_deg`: number array
- `launch_azim_deg`: number array

## `ray_paths` (optional)

Present only when `beam.store_ray_paths = true`.

Structure:

- Object keyed by `ray_<index>`
- Each value is an array of `[x, y, z]` points along that ray

Example:

```json
"ray_paths": {
  "ray_0": [[0.0, 0.0, 50.0], [5.0, 0.0, 49.8]],
  "ray_1": [[0.0, 0.0, 50.0], [5.0, 0.0, 50.2]]
}
```

## `bty`

Bathymetry vectors copied from input:

- `x_bty_m`: number array
- `y_bty_m`: number array
- `z_bty_m`: number array

## `pressure_field`

Contains per-receiver and per-frequency results.

Common keys:

- `frequency_hz`: number array
- `delay_s`: flattened 3D array container
- `amplitude`: flattened 3D array container
- `pressure_re`: flattened 4D real-part container
- `pressure_im`: flattened 4D imag-part container

Receiver coordinate keys depend on receiver mode:

- Grid mode (`receivers.config_type = "grid"`):
  - `x_m`, `y_m`, `z_m` axis vectors
- Array mode (`receivers.config_type = "array"`):
  - `receiver_positions_m` explicit receiver positions

### Flattened Array Containers

`delay_s` and `amplitude` format:

```json
{
  "shape": [nx, ny, nz],
  "data": [ ... nx*ny*nz values ... ]
}
```

`pressure_re` and `pressure_im` format:

```json
{
  "shape": [nfreq, nx, ny, nz],
  "data": [ ... nfreq*nx*ny*nz values ... ]
}
```

Flattening order is row-major nested loops:

- 3D: `ix` then `iy` then `iz`
- 4D: `ifreq` then `ix` then `iy` then `iz`

## Notes for Post-Processing

- Complex pressure is reconstructed from `pressure_re` and `pressure_im`.
- `shape` metadata is required to reshape `data` correctly.
- For broadband workflows in this repo, see the phasor-convention note in `/memories/repo/rtrs_conventions.md`.
