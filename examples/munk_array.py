import numpy as np
from time import time
import os 
import json
import python_utils
import matplotlib.pyplot as plt

z = np.linspace(0.0, 5000.0, 50)
munk_ssp = python_utils.munk(z)
munk_ssp_3d = np.tile(munk_ssp, (2, 2, 1))
munk_ssp_3d_flat = munk_ssp_3d.flatten(order='C')
name = "testm_array"
jsonfile = f"examples/{name}.json"
outfile = f"examples/{name}.out.json"
# remove outfile if present
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
        "z_bty_m": np.array([[5000.0, 5000.0], [5000.0, 5000.0]]).flatten(order='C').tolist(),
        "bottom_model": {
                "model": "elastic",
                "compressional_speed_m_s": 1700.0,
                "shear_speed_m_s": 400.0,
                "density_g_cm3": 1.6,
                "compressional_attenuation_db_per_wavelength": 0.5,
                "shear_attenuation_db_per_wavelength": 0.35,
            },
    },
    "source": {
        "position": [0.0, 0.0, 1000.0],
        "freq_hz": [500.0],
        # "freq_hz": np.linspace(-1000.0, 1000.0, 2000).tolist(),
        # "freq_hz": np.linspace(50.0, 1000.0, 1000).tolist(),
        "launch_elev_deg": np.linspace(-20.0, 20.0, 200).tolist(),
        # "launch_elev_deg": [0.0],
        # "launch_azim_deg": [0.0]
        "launch_azim_deg": np.linspace(-1.5, 1.5, 15).tolist()
    },
    "receivers": {
        "config_type": "array",
        # "x_rcvr_m": np.full(500, 0.0).tolist(),
        "x_rcvr_m": np.linspace(-200.0, 200.0, 1000).tolist(),
        # "y_rcvr_m": np.linspace(0.0, 50000.0, 500).tolist(),
        "y_rcvr_m": np.linspace(25000.0, 25200.0, 1000).tolist(),
        # "z_rcvr_m": [1000.0],
        "z_rcvr_m": np.full(1000, 1000.0).tolist()
    },
    "beam": {
        "step_m": 10.0,
        "max_steps": 100_000,
        "max_range_m": 30_000.0,
        "store_ray_paths": False,
        "show_progress": True,
        "integration_method": "rk2",
    }
}


with open(jsonfile, "w") as f:
    json.dump(env_m, f, indent=2)

t0 = time()
os.system(f"cargo run --release  {jsonfile}")
t1 = time()
print(f"Execution time: {t1 - t0:.2f} seconds")

x_bty, y_bty, z_bty = python_utils.load_bty(outfile)
freq, x_m, y_m, z_m, pressure = python_utils.load_cmpx_pressure(outfile)
try:
    pressure = np.reshape(pressure, (len(freq), len(x_m), len(y_m), len(z_m)))
    tl = - 20 * np.log10(np.abs(pressure))
    print("simulating a field in array mode, so reshaping pressure to 4D")
    z_val = 1000.0
    z_idx = (np.abs(z_m - z_val)).argmin()
    python_utils.plot_line_tl_y(tl[0,:,:,:], x_m, y_m, z_m, x_idx=len(x_m)//2, z_idx=z_idx)
    # python_utils.plot_line_tl_z(tl[0,:,:,:], x_m, y_m, z_m, x_idx=len(x_m)//2, y_idx=len(y_m)//2)
    python_utils.plot_tl_yz(tl[0,:,:,:], x_m, y_m, z_m, x_idx=len(x_m)//2)
except:
    print(f"simulating an array, not a field, so cannot reshape pressure to 4D")
    print(pressure.shape)
    print(f"array length: {np.sqrt((np.max(x_m) - np.min(x_m))**2 + (np.max(y_m) - np.min(y_m))**2)}")
    print(f"array element spacing: {np.sqrt( (x_m[1] - x_m[0])**2 + (y_m[1] - y_m[0])**2 )} m")
    tl = - 20 * np.log10(np.abs(pressure))
    python_utils.plot_tl_along_array(tl[0, :, 0, 0])
    python_utils.plot_pressure_along_array(pressure[0,:,0,0])
    python_utils.plot_array_geometry(x_m, y_m, z_m)
    




try:
    rays = python_utils.load_rays(outfile)

    # python_utils.plot_rays_xz(rays)
    # python_utils.plot_rays_yz(rays)
    python_utils.plot_rays_xy(rays)
    # python_utils.plot_rays_3d(rays)
    # python_utils.plot_rays_bty_3d(rays, x_bty, y_bty, z_bty)
except:
    print(f"No rays")




plt.show()  

