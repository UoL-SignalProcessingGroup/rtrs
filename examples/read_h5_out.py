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

def plot_tl_line(tl, points):
    rcvr_range = []
    for point in points:
        rcvr_range.append(point[0])
    plt.figure()
    plt.plot(rcvr_range, tl)
    plt.xlabel('Range (m)')
    plt.ylabel('Transmission Loss (dB)')
    plt.title(f'Transmission Loss at Receiver Depth {points[0][2]} m')
    plt.gca().invert_yaxis()
    plt.grid()

def main():
    h5file = "examples/testm.h5"
    rays = load_rays(h5file)
    plot_rays_3d(rays)
    plot_rays_xz(rays)
    plt.show()  
    

if __name__ == '__main__':
    main()
    
