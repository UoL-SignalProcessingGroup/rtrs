import json
import numpy as np

def munk(Z, min_c=1500.0, epsilon=0.00737, min_z1=1300.0, min_z2=1300.0):
    """
    Make munk ssp
    :param Z: (numpy.ndarray)
    :param min_c: (float)
    :param epsilon: (float)
    :param min_z1: (float)
    :param min_z2: (float)
    :return: ssp, c (list, numpy.ndarray)
    """

    # Calculate Munk profile
    zbar = (2 * (Z - min_z1)) / min_z2
    c = min_c * (1.0 + epsilon * (zbar - 1 + np.exp(-zbar)))

    return c


# munk stuff
z = np.linspace(0.0, 5000.0, 50)
munk_ssp = munk(z)
# print("Munk SSP (1D):", munk_ssp)
# Tile munk_ssp into a 3D array of shape (2, 2, len(munk_ssp))
munk_ssp_3d = np.tile(munk_ssp, (2, 2, 1))
# print("Munk SSP (3D):", munk_ssp_3d)
# Flatten munk_ssp_3d to a 1D array in row-major order (C order)
munk_ssp_3d_flat = munk_ssp_3d.flatten(order='C')
# print("Munk SSP (flattened 3D):", munk_ssp_3d_flat)

# perkeris
z_pekeris = np.array([0.0, 100.0])
ssp_pekeris = np.array([1500.0, 1500.0])
ssp_pekeris_3d = np.tile(ssp_pekeris, (2, 2, 1))
ssp_pekeris_3d_flat = ssp_pekeris_3d.flatten(order='C')

# print(" SSP (1D):", ssp_pekeris)
# print(" SSP (3D):", ssp_pekeris_3d)
# print(" SSP (flattened 3D):", ssp_pekeris_3d_flat)

x_rcvr = np.linspace(0.0, 5000.0, 10)
alpha = np.linspace(-10.0, 10.0, 50)



env_p = {
    "ssp": {
        "interp_type": "tabulated",
        "x_ssp_m": [0.0, 5000.0],
        "y_ssp_m": [0.0, 5000.0],
        "z_ssp_m": [0.0, 100.0],
        "c_m_s": list(ssp_pekeris_3d_flat)
    },
    "bathymetry": {
        "x_bty_m": [0.0, 5000.0],
        "y_bty_m": [0.0, 5000.0],
        "z_bty_m": np.array([[100.0, 100.0], [100.0, 100.0]]).flatten(order='C').tolist(),
        "density_g_cm3": [1.6],
        "c_bty_m_s": [1700.0],
        # "attenuation_db_per_wavelength": [0.2]
    },
    "source": {
        "position": [0.0, 0.0, 20.0],
        "freq_hz": 1000.0,
        "launch_elev_deg": np.linspace(-10.0, 10.0, 200).tolist(),
        "launch_azim_deg": [-0.1, 0.0, 0.1]
    },
    "receivers": {
        "x_rcvr_m": [0.0],
        "y_rcvr_m": np.linspace(0.0, 5000.0, 500).tolist(),
        "z_rcvr_m": np.linspace(0.0, 100.0, 50).tolist()
    },
    "beam": {
        "step_m": 10.0,
        "max_steps": 10_000,
        "max_range_m": 10_000.0
    }
}

env_m = {
    "ssp": {
        "interp_type": "tabulated",
        "x_ssp_m": [0.0, 5000.0],
        "y_ssp_m": [0.0, 5000.0],
        "z_ssp_m": list(z),
        "c_m_s": list(munk_ssp_3d_flat)
    },
    "bathymetry": {
        "x_bty_m": [0.0, 50000.0],
        "y_bty_m": [0.0, 50000.0],
        "z_bty_m": np.array([[5000.0, 5000.0], [5000.0, 5000.0]]).flatten(order='C').tolist(),
        "density_g_cm3": [1.6],
        "c_bty_m_s": [1700.0],
        # "attenuation_db_per_wavelength": [0.2]
    },
    "source": {
        "position": [0.0, 0.0, 1000.0],
        "freq_hz": 100.0,
        "launch_elev_deg": np.linspace(-10.0, 10.0, 4000).tolist(),
        "launch_azim_deg": [-0.001, 0.0, 0.001]
    },
    "receivers": {
        "x_rcvr_m": [0.0],
        "y_rcvr_m": np.linspace(0.0, 50000.0, 300).tolist(),
        "z_rcvr_m": [1000.0]
    },
    "beam": {
        "step_m": 100.0,
        "max_steps": 100_000,
        "max_range_m": 50_000.0
    }
}

# Create bathymetry (bty) array: at x=0, z=1000; at x=5000, z=0
x_bty = np.array([0.0, 50000.0])
y_bty = np.array([0.0, 50000.0])
z_bty = np.array([[1000.0, 0.0], [1000.0, 0.0]])  # shape (2,2): constant along y

# Flatten z_bty in row-major order
z_bty_flat = z_bty.flatten(order='C')

env_bty = {
    "ssp": {
        "interp_type": "tabulated",
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
        "density_g_cm3": [2.0],
        "c_bty_m_s": [1600.0],
        # "attenuation_db_per_wavelength": [0.5]
    },
    "source": {
        "position": [0.0, 25000.0, 50.0],
        "freq_hz": 1000.0,
        "launch_elev_deg": np.linspace(-50.0, 50.0, 1).tolist(),
        "launch_azim_deg": np.linspace(30.0, 150.0, 5).tolist()
    },
    "receivers": {
        "x_rcvr_m": [1000.0],
        "y_rcvr_m": [1000.0],
        "z_rcvr_m": [1000.0]
    },
    "beam": {
        "step_m": 1.0,
        "max_steps": 1_000_000,
        "max_range_m": 30_000.0
    }
}


with open("examples/testp.json", "w") as f:
    json.dump(env_p, f, indent=2)

with open("examples/testm.json", "w") as f:
    json.dump(env_m, f, indent=2)

with open("examples/testb.json", "w") as f:
    json.dump(env_bty, f, indent=2)




