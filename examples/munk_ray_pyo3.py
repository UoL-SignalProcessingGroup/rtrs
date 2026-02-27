import time
import numpy as np
import matplotlib.pyplot as plt
import rtrs
import python_utils

# build ssp for rtrs
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
        "z_bty_m": np.array([[5000.0, 5000.0],
                              [5000.0, 5000.0]]).flatten(order='C').tolist(),
    },
    "source": {
        "position": [0.0, 0.0, 1000.0],
        "freq_hz": [1000.0],
        "launch_elev_deg": np.linspace(-20.0, 20.0, 10).tolist(),
        "launch_azim_deg": [0.0]
    },
    "receivers": {
        "config_type": "grid",
        "x_rcvr_m": [0.0],
        "y_rcvr_m": np.linspace(0.0, 50000.0, 5).tolist(),
        "z_rcvr_m": np.linspace(0.0, 5000.0, 5).tolist()
    },
    "beam": {
        "step_m": 10.0,
        "max_steps": 100_000,
        "max_range_m": 55_000.0,
        "store_ray_paths": True,
    }
}

# Call the Rust engine via PyO3 binding
t0 = time.time()
result = rtrs.run_simulation(env_m)
t1 = time.time()
time_rt = t1 - t0
print(f"RTRS simulation took {time_rt:.2f} seconds")
ray_paths = result["ray_paths"]
ray_paths_np = [np.array(r) for r in ray_paths]


plt.figure()
for ray in ray_paths_np:
    y = ray[:,1]
    z = ray[:,2]
    plt.plot(y/1000, z, color="red")
plt.plot([], [], color="black", label="Bellhop")
plt.plot([], [], color="red", label="RTRS")
plt.gca().invert_yaxis()
plt.ylim([5000, 0.0])
plt.ylabel("Depth (m)")
plt.xlabel("Range (km)")
plt.title(f"ray trace")
plt.grid()
plt.legend()

plt.show()