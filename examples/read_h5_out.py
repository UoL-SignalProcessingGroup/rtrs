import h5py
import numpy as np
import matplotlib.pyplot as plt


def load_rays(h5path):
    with h5py.File(h5path, 'r') as f:
        ray_group = f['ray_paths']
        rays_path = []
        for key in ray_group.keys():
            rays_path.append(ray_group[key][:])
    return rays_path

def load_src(h5path):
    with h5py.File(h5path, 'r') as f:
        src = f['source']
        launch_elev_deg = src['launch_elev_deg'][:]
        launch_azim_deg = src['launch_azim_deg'][:]
        source_pos = src['position'][:]
        freq_hz = src['frequency_hz'][()]
    return launch_elev_deg, launch_azim_deg, source_pos, freq_hz


def load_bty(h5path):
    with h5py.File(h5path, 'r') as f:
        bty = f['bty']
        x_bty_m = bty['x_bty_m'][:]
        y_bty_m = bty['y_bty_m'][:]
        z_bty_m = bty['z_bty_m'][:]
        z_bty_m = np.reshape(z_bty_m, (len(x_bty_m), len(y_bty_m)))
    return x_bty_m, y_bty_m, z_bty_m.T


def load_cmpx_pressure(h5path):
    with h5py.File(h5path, 'r') as f:
        pressure_field = f['pressure_field']
        x_m = pressure_field['x_m'][:]
        y_m = pressure_field['y_m'][:]
        z_m = pressure_field['z_m'][:]
        pressure_im = pressure_field['pressure_im'][:]
        pressure_re = pressure_field['pressure_re'][:]
    return x_m, y_m, z_m, pressure_re + 1j * pressure_im


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

def plot_rays_xz(rays):
    plt.figure()
    for ray in rays:
        x = ray[:,1]
        z = ray[:,2]
        plt.plot(x, z, lw=0.8)
    plt.xlabel('x (m)')
    plt.ylabel('depth (m)')
    plt.title('Ray paths (x-z plane)')
    plt.gca().invert_yaxis()
    plt.grid()
    plt.tight_layout()

def plot_rays_yz(rays):
    plt.figure()
    for ray in rays:
        y = ray[:,0]
        z = ray[:,2]
        plt.plot(y, z, lw=0.8)
    plt.xlabel('y (m)')
    plt.ylabel('depth (m)')
    plt.title('Ray paths (y-z plane)')
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

def plot_tl_yz(tl, x_m, y_m, z_m, x_idx):
    plt.figure()
    Y, Z = np.meshgrid(y_m, z_m)
    plt.pcolormesh(Y, Z, tl[x_idx, :, :].T, shading='auto', cmap='jet_r')
    plt.colorbar(label='Transmission Loss (dB)')
    plt.xlabel('y (m)')
    plt.ylabel('depth (m)')
    plt.title(f'Transmission Loss (y-z plane) at x={x_m[x_idx]:.1f} m')
    plt.gca().invert_yaxis()
    plt.tight_layout()

def main():
    h5file = "examples/testm.h5"
    rays = load_rays(h5file)
    x_bty, y_bty, z_bty = load_bty(h5file)
    x_m, y_m, z_m, pressure = load_cmpx_pressure(h5file)
    pressure = np.reshape(pressure, (len(x_m), len(y_m), len(z_m)))
    print(pressure.shape)
    tl = - 20 * np.log10(np.abs(pressure))
    print(tl)
    

    plot_rays_xz(rays)
    # plot_rays_yz(rays)
    # plot_rays_xy(rays)
    # plot_rays_3d(rays)
    plot_rays_bty_3d(rays, x_bty, y_bty, z_bty)
    # plot_line_tl_x(tl, x_m, y_m, z_m, y_idx=len(y_m)//2, z_idx=len(z_m)//2)
    plot_line_tl_y(tl, x_m, y_m, z_m, x_idx=len(x_m)//2, z_idx=len(z_m)//2)
    # plot_line_tl_z(tl, x_m, y_m, z_m, x_idx=len(x_m)//2, y_idx=len(y_m)//2)
    plot_tl_yz(tl, x_m, y_m, z_m, x_idx=len(x_m)//2)
    plt.show()  
    

if __name__ == '__main__':
    main()
    
