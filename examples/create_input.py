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
z = np.linspace(0.0, 5000.0, 15)
munk_ssp = munk(z)
print("Munk SSP (1D):", munk_ssp)
# Tile munk_ssp into a 3D array of shape (2, 2, len(munk_ssp))
munk_ssp_3d = np.tile(munk_ssp, (2, 2, 1))
print("Munk SSP (3D):", munk_ssp_3d)
# Flatten munk_ssp_3d to a 1D array in row-major order (C order)
munk_ssp_3d_flat = munk_ssp_3d.flatten(order='C')
print("Munk SSP (flattened 3D):", munk_ssp_3d_flat)

# perkeris
z_pekeris = np.array([0.0, 100.0])
ssp_pekeris = np.array([1500.0, 1500.0])
ssp_pekeris_3d = np.tile(ssp_pekeris, (2, 2, 1))
ssp_pekeris_3d_flat = ssp_pekeris_3d.flatten(order='C')

print(" SSP (1D):", ssp_pekeris)
print(" SSP (3D):", ssp_pekeris_3d)
print(" SSP (flattened 3D):", ssp_pekeris_3d_flat)

rcvr_ranges = np.linspace(0.0, 5000.0, 500)
alpha = np.linspace(-20.0, 10.0, 8)



env_p = {
    "ssp": {
        "interp_type": "tabulated",
        "x_ssp_m": [0.0, 5000.0],
        "y_ssp_m": [0.0, 5000.0],
        "z_ssp_m": [0.0, 100.0],
        "c_m_s": list(ssp_pekeris_3d_flat)
    },
    "bathymetry": {
        "kind": "flat",
        "flat_depth_m": 100.0,
        "x_bty_m": None,
        "y_bty_m": None,
        "z_bty_m": None,
        "density_g_cm3": [2.0],
        "c_bty_m_s": [1600.0],
        "attenuation_db_per_wavelength": [0.5]
    },
    "source": {
        "position": [0.0, 0.0, 10.0],
        "freq_hz": 1000.0,
        "launch_elev_deg": list(alpha),
        "launch_azim_deg": [0.0]
    },
    "receivers": {
        "kind": "cylindrical",
        "ranges_m": list(rcvr_ranges),
        "bearings_deg": [0.0],
        "depths_m": [50.0],
        "x_rcvr_m": None,
        "y_rcvr_m": None,
        "z_rcvr_m": None
    },
    "beam": {
        "step_m": 1.0,
        "beam_type": "ray"
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
        "kind": "flat",
        "flat_depth_m": 5000.0,
        "x_bty_m": None,
        "y_bty_m": None,
        "z_bty_m": None,
        "density_g_cm3": [1.6],
        "c_bty_m_s": [1700.0],
        "attenuation_db_per_wavelength": [0.2]
    },
    "source": {
        "position": [0.0, 0.0, 1000.0],
        "freq_hz": 1000.0,
        "launch_elev_deg": list(alpha),
        "launch_azim_deg": [0.0]
    },
    "receivers": {
        "kind": "cylindrical",
        "ranges_m": list(rcvr_ranges),
        "bearings_deg": [0.0],
        "depths_m": [50.0],
        "x_rcvr_m": None,
        "y_rcvr_m": None,
        "z_rcvr_m": None
    },
    "beam": {
        "step_m": 1.0,
        "beam_type": "ray"
    }
}


with open("examples/testp.json", "w") as f:
    json.dump(env_p, f, indent=2)

with open("examples/testm.json", "w") as f:
    json.dump(env_m, f, indent=2)




