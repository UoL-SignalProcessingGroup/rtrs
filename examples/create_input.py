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


def linear_towed_array(center, heading_deg, length_m=None, n_elements=None, spacing_m=None, depth_m=50.0):
    """
    Generate receiver coordinates for a linear horizontal towed array.

    Parameters
    - center: (x, y, z) center position of the array (meters). z is ignored for horizontal placement
    - heading_deg: float, array heading in degrees clockwise from +x (East) (same convention as typical navigation)
    - length_m: total physical length of the array in meters (optional if spacing_m and n_elements provided)
    - n_elements: integer number of hydrophones/elements along the array (required unless length_m+spacing or spacing+makes sense)
    - spacing_m: spacing between adjacent elements in meters (optional)
    - depth_m: depth of array (positive value, meters)

    Returns
    - x_list, y_list, z_list: plain Python lists suitable for JSON serialization (meters)

    Notes
    - If n_elements and spacing_m are provided, length_m is derived as (n_elements-1)*spacing_m.
    - If length_m and n_elements are provided, spacing_m = length_m/(n_elements-1).
    - Heading 0 deg => array aligned along +x; 90 deg => along +y.
    """
    import math

    # Validate inputs and infer missing values
    if n_elements is None and spacing_m is None:
        raise ValueError("Either n_elements or spacing_m (or both with length_m) must be provided")

    if n_elements is not None and n_elements < 1:
        raise ValueError("n_elements must be >= 1")

    if n_elements is None:
        # derive n_elements from length and spacing
        if length_m is None:
            raise ValueError("When n_elements not given, length_m must be provided along with spacing_m")
        n_elements = int(math.floor(length_m / spacing_m)) + 1

    if spacing_m is None:
        if length_m is None:
            raise ValueError("When spacing_m not given, length_m and n_elements must be provided")
        if n_elements == 1:
            spacing_m = 0.0
        else:
            spacing_m = length_m / (n_elements - 1)

    # Compute element positions along a line centered on the center point
    # Parametric positions along the array from -half_len to +half_len
    half_len = spacing_m * (n_elements - 1) / 2.0
    offsets = [(-half_len + i * spacing_m) for i in range(n_elements)]

    # Convert heading to radians; heading is clockwise from +x so convert to math angle (counter-clockwise)
    # math_angle = -heading_deg
    theta = math.radians(heading_deg)

    cx, cy, cz = center
    x_list = []
    y_list = []
    z_list = []
    for s in offsets:
        # Position offset in global coords: dx = s*cos(theta), dy = s*sin(theta)
        dx = s * math.cos(theta)
        dy = s * math.sin(theta)
        x_list.append(cx + dx)
        y_list.append(cy + dy)
        # z: place at specified depth relative to sea surface (positive down)
        z_list.append(depth_m)

    return list(map(float, x_list)), list(map(float, y_list)), list(map(float, z_list))

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

env_p = {
    "ssp": {
        
        "x_ssp_m": [0.0, 10000.0],
        "y_ssp_m": [0.0, 10000.0],
        "z_ssp_m": [0.0, 100.0],
        "c_m_s": list(ssp_pekeris_3d_flat)
    },
    "bathymetry": {
        "x_bty_m": [0.0, 10000.0],
        "y_bty_m": [0.0, 10000.0],
        "z_bty_m": np.array([[100.0, 100.0], [100.0, 100.0]]).flatten(order='C').tolist(),
        "density_g_cm3": 1.6,
        "c_bty_m_s": 1700.0,
        "attenuation_db": 0.5
    },
    "source": {
        "position": [0.0, 0.0, 25.0],
        "freq_hz": [1000.0],
        "launch_elev_deg": np.linspace(-10.0, 10.0, 200).tolist(),
        "launch_azim_deg": [-0.1, 0.0, 0.1]
    },
    "receivers": {
        "x_rcvr_m": [0.0],
        "y_rcvr_m": np.linspace(0.0, 10000.0, 500).tolist(),
        "z_rcvr_m": np.linspace(0.0, 100.0, 100).tolist()
    },
    "beam": {
        "step_m": 15.0,
        "max_steps": 10_000,
        "max_range_m": 10_000.0
    }
}

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
        "attenuation_db": 1.0
    },
    "source": {
        "position": [0.0, 0.0, 1000.0],
        "freq_hz": [1000.0],
        # "freq_hz": np.linspace(1.0, 1000.0, 1000).tolist(),
        "launch_elev_deg": np.linspace(-10.0, 10.0, 4000).tolist(),
        "launch_azim_deg": [-0.001, 0.0, 0.001]
    },
    "receivers": {
        "x_rcvr_m": [0.0],
        "y_rcvr_m": np.linspace(0.0, 50000.0, 2000).tolist(),
        "z_rcvr_m": [1000.0] # np.linspace(0.0, 5000.0, 500).tolist()
    },
    "beam": {
        "step_m": 100.0,
        "max_steps": 100_000,
        "max_range_m": 50_000.0
    }
}

env_m2d = {
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
        "attenuation_db": 1.0
    },
    "source": {
        "position": [0.0, 0.0, 1000.0],
        "freq_hz": [1000.0],
        # "freq_hz": np.linspace(1.0, 1000.0, 1000).tolist(),
        "launch_elev_deg": np.linspace(-10.0, 10.0, 4000).tolist(),
        "launch_azim_deg": [-0.001, 0.0, 0.001]
    },
    "receivers": {
        "x_rcvr_m": [0.0],
        "y_rcvr_m": np.linspace(0.0, 50000.0, 2000).tolist(),
        "z_rcvr_m": [1000.0] # np.linspace(0.0, 5000.0, 500).tolist()
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
        "density_g_cm3": 2.0,
        "c_bty_m_s": 1600.0,
        "attenuation_db": 0.5
    },
    "source": {
        "position": [0.0, 25000.0, 50.0],
        "freq_hz": [1000.0],
        "launch_elev_deg": np.linspace(-25.0, 25.0, 1).tolist(),
        "launch_azim_deg": np.linspace(30.0, 150.0, 5).tolist()
    },
    "receivers": {
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

x,y,z = linear_towed_array(center=(0.0, 5000.0, 50.0), heading_deg=45.0, length_m=100.0, n_elements=11, depth_m=50.0)
print()
print(x)
print(y)
print(z)

env_lha = {
    "ssp": {
        "x_ssp_m": [0.0, 10000.0],
        "y_ssp_m": [0.0, 10000.0],
        "z_ssp_m": [0.0, 100.0],
        "c_m_s": list(ssp_pekeris_3d_flat)
    },
    "bathymetry": {
        "x_bty_m": [0.0, 10000.0],
        "y_bty_m": [0.0, 10000.0],
        "z_bty_m": np.array([[100.0, 100.0], [100.0, 100.0]]).flatten(order='C').tolist(),
        "density_g_cm3": 1.6,
        "c_bty_m_s": 1700.0,
        "attenuation_db": 0.5
    },
    "source": {
        "position": [0.0, 0.0, 25.0],
        "freq_hz": [1000.0],
        "launch_elev_deg": np.linspace(-10.0, 10.0, 200).tolist(),
        "launch_azim_deg": [-0.1, 0.0, 0.1]
    },
    "receivers": {
        "x_rcvr_m": x,
        "y_rcvr_m": y,
        "z_rcvr_m": z
    },
    "beam": {
        "step_m": 15.0,
        "max_steps": 10_000,
        "max_range_m": 10_000.0
    }
}

env_bb = {
    "ssp": {
        "x_ssp_m": [0.0, 10000.0],
        "y_ssp_m": [0.0, 10000.0],
        "z_ssp_m": [0.0, 100.0],
        "c_m_s": list(ssp_pekeris_3d_flat)
    },
    "bathymetry": {
        "x_bty_m": [0.0, 10000.0],
        "y_bty_m": [0.0, 10000.0],
        "z_bty_m": np.array([[100.0, 100.0], [100.0, 100.0]]).flatten(order='C').tolist(),
        "density_g_cm3": 1.6,
        "c_bty_m_s": 1700.0,
        "attenuation_db": 0.5
    },
    "source": {
        "position": [0.0, 0.0, 25.0],
        "freq_hz": np.linspace(1.0, 1000.0, 1000).tolist(),
        "launch_elev_deg": np.linspace(-10.0, 10.0, 200).tolist(),
        "launch_azim_deg": [-0.1, 0.0, 0.1]
    },
    "receivers": {
        "x_rcvr_m": [0.0],
        "y_rcvr_m": [30000.0],
        "z_rcvr_m": [20.0]
    },
    "beam": {
        "step_m": 15.0,
        "max_steps": 40_000,
        "max_range_m": 40_000.0
    }
}


with open("examples/testp.json", "w") as f:
    json.dump(env_p, f, indent=2)

with open("examples/testm.json", "w") as f:
    json.dump(env_m, f, indent=2)

with open("examples/testb.json", "w") as f:
    json.dump(env_bty, f, indent=2)

with open("examples/testlha.json", "w") as f:
    json.dump(env_lha, f, indent=2)

with open("examples/testbb.json", "w") as f:
    json.dump(env_bb, f, indent=2)

with open("examples/testm2d.json", "w") as f:
    json.dump(env_m2d, f, indent=2)




