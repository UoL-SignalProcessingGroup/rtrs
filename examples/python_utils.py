import json
import numpy as np
import matplotlib.pyplot as plt

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

def load_rays(jpath):
    with open(jpath, 'r') as f:
        data = json.load(f)
    ray_group = data['ray_paths']
    rays_path = []
    for key in sorted(ray_group.keys(), key=lambda k: int(k.split('_')[1])):
        rays_path.append(np.array(ray_group[key], dtype=np.float32))
    return rays_path

def load_src(jpath):
    with open(jpath, 'r') as f:
        data = json.load(f)
    src = data['src']
    launch_elev_deg = np.array(src['launch_elev_deg'])
    launch_azim_deg = np.array(src['launch_azim_deg'])
    source_pos = np.array(src['source_position_m'])
    freq_hz = np.array(src['frequency_hz'])
    return launch_elev_deg, launch_azim_deg, source_pos, freq_hz


def load_bty(jpath):
    with open(jpath, 'r') as f:
        data = json.load(f)
    bty = data['bty']
    x_bty_m = np.array(bty['x_bty_m'])
    y_bty_m = np.array(bty['y_bty_m'])
    z_bty_m = np.array(bty['z_bty_m'])
    z_bty_m = np.reshape(z_bty_m, (len(x_bty_m), len(y_bty_m)))
    return x_bty_m, y_bty_m, z_bty_m.T


def load_cmpx_pressure(jpath):
    with open(jpath, 'r') as f:
        data = json.load(f)
    pf = data['pressure_field']
    frequency_hz = np.array(pf['frequency_hz'])
    if 'receiver_positions_m' in pf:
        recs = np.array(pf['receiver_positions_m'])
        x_m = np.unique(recs[:, 0])
        y_m = np.unique(recs[:, 1])
        z_m = np.unique(recs[:, 2])
    else:
        x_m = np.array(pf['x_m'])
        y_m = np.array(pf['y_m'])
        z_m = np.array(pf['z_m'])
    re_entry = pf['pressure_re']
    im_entry = pf['pressure_im']
    pressure_re = np.array(re_entry['data'], dtype=np.float32).reshape(re_entry['shape'])
    pressure_im = np.array(im_entry['data'], dtype=np.float32).reshape(im_entry['shape'])
    return frequency_hz, x_m, y_m, z_m, pressure_re + 1j * pressure_im



def plot_rays_3d(rays):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for ray in rays:
        x = ray[:,0]
        y = ray[:,1]
        z = ray[:,2]
        ax.plot(x, y, z, lw=0.8)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('depth (m)')
    ax.set_title('Ray paths (3D)')
    ax.invert_zaxis()
    plt.tight_layout()

def plot_rays_bty_3d(rays, x_bty, y_bty, z_bty):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # Plot rays
    for ray in rays:
        x = ray[:,0]
        y = ray[:,1]
        z = ray[:,2]
        ax.plot(x, y, z, lw=0.8)
    # Plot bathymetry surface
    X, Y = np.meshgrid(x_bty, y_bty)
    ax.plot_surface(X, Y, z_bty, cmap='terrain', alpha=0.5, linewidth=0, antialiased=False)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('depth (m)')
    ax.set_title('Ray paths (3D) with Bathymetry')
    ax.set_xlim(np.min(x_bty), np.max(x_bty))
    ax.set_ylim(np.min(y_bty), np.max(y_bty))
    ax.invert_zaxis()
    plt.tight_layout()

def plot_rays_yz(rays):
    plt.figure()
    for ray in rays:
        x = ray[:,1]
        z = ray[:,2]
        plt.plot(x, z, lw=0.8)
    plt.xlabel('y (m)')
    plt.ylabel('depth (m)')
    plt.title('Ray paths (y-z plane)')
    plt.gca().invert_yaxis()
    plt.grid()
    plt.tight_layout()

def plot_rays_yz_bty(rays, x_bty, y_bty, z_bty, x_idx=None):
    plt.figure()
    for ray in rays:
        y = ray[:,1]
        z = ray[:,2]
        plt.plot(y, z, lw=0.8)

    if x_idx is None:
        z_profile = np.min(z_bty, axis=1)
        bty_label = 'bathymetry (projected min depth over x)'
    else:
        z_profile = z_bty[:, x_idx]
        bty_label = f'bathymetry at x={x_bty[x_idx]:.1f} m'

    plt.plot(y_bty, z_profile, 'k-', lw=1.8, label=bty_label)
    plt.xlabel('y (m)')
    plt.ylabel('depth (m)')
    plt.title('Ray paths (y-z plane) with Bathymetry')
    plt.gca().invert_yaxis()
    plt.grid()
    plt.legend()
    plt.tight_layout()

def plot_rays_xz(rays):
    plt.figure()
    for ray in rays:
        y = ray[:,0]
        z = ray[:,2]
        plt.plot(y, z, lw=0.8)
    plt.xlabel('x (m)')
    plt.ylabel('depth (m)')
    plt.title('Ray paths (x-z plane)')
    plt.gca().invert_yaxis()
    plt.grid()
    plt.tight_layout()

def plot_rays_xy(rays):
    plt.figure()
    for ray in rays:
        x = ray[:,0]
        y = ray[:,1]
        plt.plot(x, y, lw=0.8)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Ray paths (x-y plane)')
    plt.axis('equal')
    plt.grid()
    plt.tight_layout()

def plot_line_tl_x(tl, x_m, y_m, z_m, y_idx, z_idx):
    plt.figure()
    plt.plot(x_m, tl[:, y_idx, z_idx])
    plt.xlabel('x (m)')
    plt.ylabel('Transmission Loss (dB)')
    plt.title(f'Transmission Loss along x at y={y_m[y_idx]:.1f} m, z={z_m[z_idx]:.1f} m')
    plt.gca().invert_yaxis()
    plt.grid()
    plt.tight_layout()

def plot_line_tl_y(tl, x_m, y_m, z_m, x_idx, z_idx):
    plt.figure()
    plt.plot(y_m, tl[x_idx, :, z_idx])
    plt.xlabel('y (m)')
    plt.ylabel('Transmission Loss (dB)')
    plt.title(f'Transmission Loss along y at x={x_m[x_idx]:.1f} m, z={z_m[z_idx]:.1f} m')
    plt.grid()
    plt.gca().invert_yaxis()
    plt.tight_layout()

def plot_line_tl_z(tl, x_m, y_m, z_m, x_idx, y_idx):
    plt.figure()
    plt.plot(z_m, tl[x_idx, y_idx, :])
    plt.xlabel('depth (m)')
    plt.ylabel('Transmission Loss (dB)')
    plt.title(f'Transmission Loss along z at x={x_m[x_idx]:.1f} m, y={y_m[y_idx]:.1f} m')
    plt.gca().invert_yaxis()
    plt.grid()
    plt.tight_layout()

def plot_tl_yz(tl, x_m, y_m, z_m, x_idx, vmin=50, vmax=150):
    plt.figure()
    Y, Z = np.meshgrid(y_m, z_m)
    plt.pcolormesh(Y, Z, tl[x_idx, :, :].T, shading='auto', cmap='jet_r', vmin=vmin, vmax=vmax)
    plt.colorbar(label='Transmission Loss (dB)')
    plt.xlabel('y (m)')
    plt.ylabel('depth (m)')
    plt.title(f'Transmission Loss (y-z plane) at x={x_m[x_idx]:.1f} m')
    plt.gca().invert_yaxis()
    plt.tight_layout()

def plot_pressure_yz(pressure, x_m, y_m, z_m, x_idx):
    plt.figure()
    Y, Z = np.meshgrid(y_m, z_m)
    plt.pcolormesh(Y, Z, pressure[x_idx, :, :].T, shading='auto', cmap='bwr', vmin=-0.0001, vmax=0.0001)
    plt.colorbar(label='Pressure')
    plt.xlabel('y (m)')
    plt.ylabel('depth (m)')
    plt.title(f'Pressure (y-z plane) at x={x_m[x_idx]:.1f} m')
    plt.gca().invert_yaxis()
    plt.tight_layout()


def plot_pressure_freq(pressure, freq, x_m, y_m, z_m, x_idx, y_idx, z_idx):
    plt.figure()
    plt.plot(freq, pressure[:, x_idx, y_idx, z_idx])
    plt.xlabel('Frequency Index')
    plt.ylabel('Pressure Amplitude')
    plt.title(f'Pressure Amplitude at x={x_m[x_idx]:.1f} m, y={y_m[y_idx]:.1f} m, z={z_m[z_idx]:.1f} m')
    plt.grid()
    plt.tight_layout()