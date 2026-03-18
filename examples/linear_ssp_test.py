import numpy as np
import os
import json
import python_utils
import matplotlib.pyplot as plt

# Linear SSP: c(z) = 1500 + 0.02 * z  [m/s]
water_depth = 150.0
ssp_gradient = 0.02  # s^-1

z_ssp = np.array([0.0, water_depth])
c_ssp = 1500.0 + ssp_gradient * z_ssp  # [1500.0, 1503.0]

ssp_3d = np.tile(c_ssp, (2, 2, 1))
ssp_3d_flat = ssp_3d.flatten(order='C')

name = "testlinear"
jsonfile = f"examples/{name}.json"
outfile = f"examples/{name}.out.json"

if os.path.exists(outfile):
    os.remove(outfile)

env = {
    "ssp": {
        "x_ssp_m": [0.0, 30000.0],
        "y_ssp_m": [0.0, 30000.0],
        "z_ssp_m": z_ssp.tolist(),
        "c_m_s": ssp_3d_flat.tolist()
    },
    "bathymetry": {
        "x_bty_m": [0.0, 30000.0],
        "y_bty_m": [0.0, 30000.0],
        "z_bty_m": np.array([[water_depth, water_depth],
                              [water_depth, water_depth]]).flatten(order='C').tolist(),
        "bottom_model": {
            "model": "acoustic",
            "compressional_speed_m_s": 1600.0,
            "density_g_cm3": 1.5,
            "compressional_attenuation_db_per_wavelength": 0.2
        },
        "water_density_g_cm3": 1.0,              # g/cm3
    },
    "source": {
        "position": [0.0, 0.0, 50.0],
        "freq_hz": [100.0],
        "launch_elev_deg": np.linspace(-15.0, 15.0, 15).tolist(),
        "launch_azim_deg": np.linspace(-0.1, 0.1, 3).tolist()
    },
    "receivers": {
        "config_type": "grid",
        "x_rcvr_m": [0.0],
        "y_rcvr_m": np.linspace(0.0, 30000.0, 500).tolist(),
        "z_rcvr_m": np.linspace(0.0, water_depth, 100).tolist()
    },
    "beam": {
        "step_m": 10.0,
        "max_steps": 10_000,
        "max_range_m": 30_000.0,
        "store_ray_paths": True,
        "integration_method": "rk2",
    }
}

with open(jsonfile, "w") as f:
    json.dump(env, f, indent=2)

os.system(f"cargo run --release  {jsonfile}")

rays = python_utils.load_rays(outfile)
x_bty, y_bty, z_bty = python_utils.load_bty(outfile)
freq, x_m, y_m, z_m, pressure = python_utils.load_cmpx_pressure(outfile)
pressure = np.reshape(pressure, (len(freq), len(x_m), len(y_m), len(z_m)))
tl = -20 * np.log10(np.abs(pressure))

python_utils.plot_rays_xz(rays)
python_utils.plot_rays_yz(rays)
python_utils.plot_line_tl_x(tl[0,:,:,:], x_m, y_m, z_m, y_idx=len(y_m)//2, z_idx=len(z_m)//2)
python_utils.plot_line_tl_y(tl[0,:,:,:], x_m, y_m, z_m, x_idx=len(x_m)//2, z_idx=len(z_m)//2)
python_utils.plot_line_tl_z(tl[0,:,:,:], x_m, y_m, z_m, x_idx=len(x_m)//2, y_idx=len(y_m)//2)
python_utils.plot_tl_yz(tl[0,:,:,:], x_m, y_m, z_m, x_idx=len(x_m)//2)
python_utils.plot_pressure_freq(pressure, freq, x_m, y_m, z_m, x_idx=-1, y_idx=-1, z_idx=len(z_m)//2)
plt.show()
