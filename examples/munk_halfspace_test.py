import numpy as np
import os
import json
import python_utils
import matplotlib.pyplot as plt
import time


z = np.linspace(0.0, 5000.0, 50)
munk_ssp = python_utils.munk(z)
munk_ssp_3d = np.tile(munk_ssp, (2, 2, 1))
munk_ssp_3d_flat = munk_ssp_3d.flatten(order='C')

name     = "testm_halfspace"
jsonfile = f"examples/{name}.json"
outfile  = f"examples/{name}.out.json"

if os.path.exists(outfile):
    os.remove(outfile)


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
        "z_bty_m": np.array([[5000.0, 5000.0],
                              [5000.0, 5000.0]]).flatten(order='C').tolist(),
        # Halfspace parameters
        "bottom_p_wave_speed_m_s": 1600.0,      # bottom (m/s)
        "bottom_density_g_cm3": 1.8,             # g/cm3
        "water_density_g_cm3": 1.0,              # g/cm3
    },
    "source": {
        "position": [0.0, 0.0, 1000.0],
        "freq_hz": [200.0],
        "launch_elev_deg": np.linspace(-50.0, 50.0, 500).tolist(),
        "launch_azim_deg": np.linspace(-0.5, 0.5, 3).tolist()
    },
    "receivers": {
        "config_type": "grid",
        "x_rcvr_m": [0.0],
        "y_rcvr_m": np.linspace(0.0, 50000.0, 500).tolist(),
        "z_rcvr_m": np.linspace(0.0, 5000.0, 500).tolist()
    },
    "beam": {
        "step_m": 10.0,
        "max_steps": 100_000,
        "max_range_m": 50_000.0
    }
}

with open(jsonfile, "w") as f:
    json.dump(env_m, f, indent=2)


t0 = time.time()
os.system(f"cargo run --release {jsonfile}")
t1 = time.time()
print(f"completed in {t1 - t0:.2f} seconds.")


rays = python_utils.load_rays(outfile)
x_bty, y_bty, z_bty = python_utils.load_bty(outfile)
freq, x_m, y_m, z_m, pressure = python_utils.load_cmpx_pressure(outfile)
pressure = np.reshape(pressure, (len(freq), len(x_m), len(y_m), len(z_m)))
tl = -20 * np.log10(np.abs(pressure))

python_utils.plot_rays_xz(rays)
python_utils.plot_rays_yz(rays)

z_val = 1000.0
z_idx = (np.abs(z_m - z_val)).argmin()
python_utils.plot_line_tl_y(tl[0, :, :, :], x_m, y_m, z_m,
                            x_idx=len(x_m) // 2, z_idx=z_idx)
python_utils.plot_tl_yz(tl[0, :, :, :], x_m, y_m, z_m, x_idx=len(x_m) // 2)
python_utils.plot_pressure_yz(np.real(pressure[0, :, :, :]), x_m, y_m, z_m,
                               x_idx=len(x_m) // 2)
plt.show()
