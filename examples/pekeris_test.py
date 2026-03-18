import numpy as np
import os 
import json
import python_utils
import matplotlib.pyplot as plt

z_pekeris = np.array([0.0, 100.0])
ssp_pekeris = np.array([1500.0, 1500.0])
ssp_pekeris_3d = np.tile(ssp_pekeris, (2, 2, 1))
ssp_pekeris_3d_flat = ssp_pekeris_3d.flatten(order='C')
name = "testp"
jsonfile = f"examples/{name}.json"
outfile = f"examples/{name}.out.json"
# remove outfile if present
if os.path.exists(outfile):
    os.remove(outfile)

env_p = {
    "ssp": {
        "x_ssp_m": [0.0, 30000.0],
        "y_ssp_m": [0.0, 30000.0],
        "z_ssp_m": [0.0, 100.0],
        "c_m_s": list(ssp_pekeris_3d_flat)
    },
    "bathymetry": {
        "x_bty_m": [0.0, 30000.0],
        "y_bty_m": [0.0, 30000.0],
        "z_bty_m": np.array([[100.0, 100.0], [100.0, 100.0]]).flatten(order='C').tolist(),
        "bottom_model": {
            "model": "acoustic",
            "compressional_speed_m_s": 1600.0,
            # "shear_speed_m_s": 400.0,
            "density_g_cm3": 1.5,
            "compressional_attenuation_db_per_wavelength": 0.2,
            # "shear_attenuation_db_per_wavelength": 0.35
        }
    },
    "source": {
        "position": [0.0, 0.0, 25.0],
        "freq_hz": [50.0],
        "launch_elev_deg": np.linspace(-15.0, 15.0, 500).tolist(),
        "launch_azim_deg": np.linspace(-0.1, 0.1, 3).tolist()
    },
    "receivers": {
        "config_type": "grid",
        "x_rcvr_m": [0.0],
        "y_rcvr_m": np.linspace(0.0, 30000.0, 500).tolist(),
        "z_rcvr_m": np.linspace(0.0, 100.0, 100).tolist()
    },
    "beam": {
        "step_m": 20.0,
        "max_steps": 500_000,
        "max_range_m": 30_000.0,
        "store_ray_paths": False,
    }
}


with open(jsonfile, "w") as f:
    json.dump(env_p, f, indent=2)

os.system(f"cargo run --release  {jsonfile}")

# rays = python_utils.load_rays(outfile)
x_bty, y_bty, z_bty = python_utils.load_bty(outfile)
freq, x_m, y_m, z_m, pressure = python_utils.load_cmpx_pressure(outfile)
pressure = np.reshape(pressure, (len(freq), len(x_m), len(y_m), len(z_m)))
tl = - 20 * np.log10(np.abs(pressure))


# python_utils.plot_rays_xz(rays)
# python_utils.plot_rays_yz(rays)
# python_utils.plot_rays_xy(rays)
# python_utils.plot_rays_3d(rays)
# python_utils.plot_rays_bty_3d(rays, x_bty, y_bty, z_bty)
# python_utils.plot_line_tl_x(tl[0,:,:,:], x_m, y_m, z_m, y_idx=len(y_m)//2, z_idx=len(z_m)//2)
python_utils.plot_line_tl_y(tl[0,:,:,:], x_m, y_m, z_m, x_idx=len(x_m)//2, z_idx=len(z_m)//2)
# python_utils.plot_line_tl_z(tl[0,:,:,:], x_m, y_m, z_m, x_idx=len(x_m)//2, y_idx=len(y_m)//2)
python_utils.plot_tl_yz(tl[0,:,:,:], x_m, y_m, z_m, x_idx=len(x_m)//2, vmax=120)
# python_utils.plot_pressure_freq(pressure, freq, x_m, y_m, z_m, x_idx=-1, y_idx=-1, z_idx=len(z_m)//2)
plt.show()  

