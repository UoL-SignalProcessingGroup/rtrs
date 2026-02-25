"""
Example using the PyO3 binding rtrs.run_simulation to run the Munk SSP test in-memory.
Build and install the extension first with:

    pip install maturin
    maturin develop --release --features python

Then run this script with the Python environment where the extension was installed.
"""

import numpy as np
import python_utils
import matplotlib.pyplot as plt
import rtrs

# build the same configuration as examples/munk_test.py
z = np.linspace(0.0, 5000.0, 50)
munk_ssp = python_utils.munk(z)
munk_ssp_3d = np.tile(munk_ssp, (2, 2, 1))
munk_ssp_3d_flat = munk_ssp_3d.flatten(order='C')

env_m = {
    "ssp": {
        "x_ssp_m": [0.0, 5000.0],
        "y_ssp_m": [0.0, 5000.0],
        "z_ssp_m": list(z),
        "c_m_s": list(munk_ssp_3d_flat)
    },
    "bathymetry": {
        "x_bty_m": [0.0, 50000.0],
        "y_bty_m": [0.0, 50000.0],
        "z_bty_m": np.array([[5000.0, 5000.0], [5000.0, 5000.0]]).flatten(order='C').tolist(),
        "density_g_cm3": 1.6,
        "c_bty_m_s": 1700.0,
        "attenuation_db": 0.0
    },
    "source": {
        "position": [0.0, 0.0, -1000.0],
        "freq_hz": [50.0],
        # "freq_hz": np.linspace(1.0, 1000.0, 1000).tolist(),
        "launch_elev_deg": np.linspace(-20.0, 20.0, 200).tolist(),
        # "launch_elev_deg": [0.0],
        # "launch_azim_deg": [0.0]
        "launch_azim_deg": np.linspace(-0.5, 0.5, 3).tolist()
    },
    "receivers": {
        "config_type": "array",
        "x_rcvr_m": np.full(500, 0.0).tolist(),
        "y_rcvr_m": np.linspace(0.0, 50000.0, 500).tolist(),
        # "z_rcvr_m": [1000.0] 
        "z_rcvr_m": np.full(500, 1000.0).tolist()
    },
    "beam": {
        "step_m": 10.0,
        "max_steps": 100_000,
        "max_range_m": 50_000.0
    }
}

# Call the Rust engine via PyO3 binding
result = rtrs.run_simulation(env_m)

# Extract ray paths and pressure_field
ray_paths = result["ray_paths"]
pf = result["pressure_field"]

freq = pf["frequency_hz"]
receiver_positions_m = np.array(pf['receiver_positions_m'][:])
shape = tuple(pf["shape"])  # For array mode: (nfreq, nreceivers, 1, 1)

# Reconstruct complex pressure array from flat re/im
re = np.array(pf["pressure_re"], dtype=np.float32).reshape(shape)
im = np.array(pf["pressure_im"], dtype=np.float32).reshape(shape)
pressure_raw = re + 1j * im

# Extract unique coordinates (matching python_utils.load_cmpx_pressure behavior)
x_m = np.unique(receiver_positions_m[:, 0])
y_m = np.unique(receiver_positions_m[:, 1])
z_m = np.unique(receiver_positions_m[:, 2])

# Reshape pressure to match expected format: (nfreq, nx, ny, nz)
pressure = np.reshape(pressure_raw, (len(freq), len(x_m), len(y_m), len(z_m)))
tl = - 20 * np.log10(np.abs(pressure))
# print(pressure)

# Convert returned ray_paths (lists) into numpy arrays expected by python_utils
ray_paths_np = [np.array(r) for r in ray_paths]

# Use the same plotting functions as munk_array.py
python_utils.plot_rays_xz(ray_paths_np)
python_utils.plot_rays_yz(ray_paths_np)
# python_utils.plot_rays_xy(ray_paths_np)
# python_utils.plot_rays_3d(ray_paths_np)

# Plot TL along y at z=1000m
z_val = 1000.0
z_idx = (np.abs(z_m - z_val)).argmin()
python_utils.plot_line_tl_y(tl[0,:,:,:], x_m, y_m, z_m, x_idx=len(x_m)//2, z_idx=z_idx)
# python_utils.plot_line_tl_z(tl[0,:,:,:], x_m, y_m, z_m, x_idx=len(x_m)//2, y_idx=len(y_m)//2)
# python_utils.plot_tl_yz(tl[0,:,:,:], x_m, y_m, z_m, x_idx=len(x_m)//2)
# python_utils.plot_pressure_yz(np.real(pressure[0,:,:,:]), x_m, y_m, z_m, x_idx=len(x_m)//2)
# python_utils.plot_pressure_freq(pressure, freq, x_m, y_m, z_m, x_idx=-1, y_idx=-1, z_idx=len(z_m)//2)
plt.show()
