import numpy as np
import os 
import json
import python_utils
import matplotlib.pyplot as plt

z_pekeris = np.array([0.0, 100.0])
ssp_pekeris = np.array([1500.0, 1500.0])
ssp_pekeris_3d = np.tile(ssp_pekeris, (2, 2, 1))
ssp_pekeris_3d_flat = ssp_pekeris_3d.flatten(order='C')
name = "testbty"
jsonfile = f"examples/{name}.json"
outfile = f"examples/{name}.out.json"
# remove outfile if present
if os.path.exists(outfile):
    os.remove(outfile)

# Create bathymetry (bty) array: at x=0, z=1000; at x=5000, z=0
x_bty = np.array([0.0, 50000.0])
y_bty = np.array([0.0, 50000.0])
z_bty = np.array([[1000.0, -100.0], [1000.0, -100.0]])  # shape (2,2): constant along y

# Flatten z_bty in row-major order
z_bty_flat = z_bty.flatten(order='C')

env_bty = {
    "ssp": {
        "x_ssp_m": [0.0, 5000.0],
        "y_ssp_m": [0.0, 5000.0],
        "z_ssp_m": [0.0, 100.0],
        "c_m_s": list(ssp_pekeris_3d_flat)
    },
    "bathymetry": {
        "kind": "tabulated",
        "flat_depth_m": None,
        "x_bty_m": list(x_bty),
        "y_bty_m": list(y_bty),
        "z_bty_m": list(z_bty_flat),
    },
    "source": {
        "position": [0.0, 25000.0, 50.0],
        "freq_hz": [1000.0],
        "launch_elev_deg": np.linspace(-25.0, 25.0, 1).tolist(),
        "launch_azim_deg": np.linspace(30.0, 150.0, 5).tolist()
    },
    "receivers": {
        "config_type": "grid",
        "x_rcvr_m": [1000.0],
        "y_rcvr_m": [1000.0],
        "z_rcvr_m": [1000.0]
    },
    "beam": {
        "step_m": 25.0,
        "max_steps": 1_000_000,
        "max_range_m": 30_000.0
    }
}

with open(jsonfile, "w") as f:
    json.dump(env_bty, f, indent=2)

os.system(f"cargo run --release  {jsonfile}")

rays = python_utils.load_rays(outfile)
x_bty, y_bty, z_bty = python_utils.load_bty(outfile)
freq, x_m, y_m, z_m, pressure = python_utils.load_cmpx_pressure(outfile)
pressure = np.reshape(pressure, (len(freq), len(x_m), len(y_m), len(z_m)))
tl = - 20 * np.log10(np.abs(pressure))


python_utils.plot_rays_xz(rays)
python_utils.plot_rays_yz(rays)
python_utils.plot_rays_xy(rays)
# python_utils.plot_rays_3d(rays)
python_utils.plot_rays_bty_3d(rays, x_bty, y_bty, z_bty)
# python_utils.plot_line_tl_x(tl[0,:,:,:], x_m, y_m, z_m, y_idx=len(y_m)//2, z_idx=len(z_m)//2)
# python_utils.plot_line_tl_y(tl[0,:,:,:], x_m, y_m, z_m, x_idx=len(x_m)//2, z_idx=len(z_m)//2)
# python_utils.plot_line_tl_z(tl[0,:,:,:], x_m, y_m, z_m, x_idx=len(x_m)//2, y_idx=len(y_m)//2)
# python_utils.plot_tl_yz(tl[0,:,:,:], x_m, y_m, z_m, x_idx=len(x_m)//2)
# python_utils.plot_pressure_freq(pressure, freq, x_m, y_m, z_m, x_idx=-1, y_idx=-1, z_idx=len(z_m)//2)
plt.show()  

