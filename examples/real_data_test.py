import json
import os
import time

import cmocean
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

import python_utils


GEBCO_FILE = "bins/env_data/GEBCO_11_Apr_2025_cd1b685d47c9/gebco_2024_n61.25_s59.0_w-8.0_e-5.0.nc"
COP_TEMP_FILE = "bins/env_data/cmems_mod_glo_phy-thetao_anfc_0.083deg_PT6H-i_thetao_8.00W-5.00W_59.00N-61.25N_0.49-5274.78m_2025-03-19.nc"
COP_SAL_FILE = "bins/env_data/cmems_mod_glo_phy-so_anfc_0.083deg_PT6H-i_so_8.00W-5.00W_59.00N-61.25N_0.49-5274.78m_2025-03-19.nc"
LAUNCH_SSP_DASHBOARD = False
SSP_DASHBOARD_PORT = 8050


def latlon_to_xy_m(lat_deg, lon_deg, lat0_deg, lon0_deg):
    lat0_rad = np.deg2rad(lat0_deg)
    dlon_rad = np.deg2rad(lon_deg - lon0_deg)
    dlat_rad = np.deg2rad(lat_deg - lat0_deg)

    a = 6_378_137.0
    f = 1.0 / 298.257223563
    e2 = f * (2.0 - f)

    sin_lat0 = np.sin(lat0_rad)
    w = np.sqrt(1.0 - e2 * sin_lat0**2)
    n = a / w
    m = a * (1.0 - e2) / (w**3)

    x = dlon_rad * n * np.cos(lat0_rad)
    y = dlat_rad * m
    return x, y


def leroy_sound_speed(z, temp, sal, lat_deg):
    z4 = z[None, :, None, None]
    lat4 = lat_deg[None, None, :, None]
    c = (
        1402.5 + 5.0 * temp - 5.44e-2 * temp**2 + 2.1e-4 * temp**3
        + 1.33 * sal - 1.23e-2 * sal * temp + 8.7e-5 * sal * temp**2
        + 1.56e-2 * z4 + 2.55e-7 * z4**2 - 7.3e-12 * z4**3
        + 1.2e-6 * z4 * (lat4 - 45.0) - 9.5e-13 * temp * z4**3
        + 3e-7 * temp**2 * z4 + 1.43e-5 * sal * z4
    )
    return c


def fill_nans_2d(arr):
    out = np.array(arr, dtype=float, copy=True)
    ny, nx = out.shape
    for j in range(ny):
        row = out[j, :]
        valid = np.isfinite(row)
        if valid.any():
            out[j, :] = np.interp(np.arange(nx), np.where(valid)[0], row[valid])
    for i in range(nx):
        col = out[:, i]
        valid = np.isfinite(col)
        if valid.any():
            out[:, i] = np.interp(np.arange(ny), np.where(valid)[0], col[valid])
    return out


def interp_2d_regular(z_old, x_old, y_old, x_new, y_new):
    z_old = np.array(z_old, dtype=float)
    z_x = np.vstack([np.interp(x_new, x_old, row) for row in z_old])
    z_xy = np.vstack([np.interp(y_new, y_old, z_x[:, i]) for i in range(z_x.shape[1])]).T
    return z_xy


def interp_2d_nearest_regular(z_old, x_old, y_old, x_new, y_new):
    z_old = np.asarray(z_old, dtype=float)
    x_old = np.asarray(x_old, dtype=float)
    y_old = np.asarray(y_old, dtype=float)
    x_new = np.asarray(x_new, dtype=float)
    y_new = np.asarray(y_new, dtype=float)

    ix = np.searchsorted(x_old, x_new)
    ix = np.clip(ix, 1, len(x_old) - 1)
    ix = np.where(np.abs(x_new - x_old[ix - 1]) <= np.abs(x_old[ix] - x_new), ix - 1, ix)

    iy = np.searchsorted(y_old, y_new)
    iy = np.clip(iy, 1, len(y_old) - 1)
    iy = np.where(np.abs(y_new - y_old[iy - 1]) <= np.abs(y_old[iy] - y_new), iy - 1, iy)

    return z_old[np.ix_(iy, ix)]


def interp_3d_horizontal(c_zyx, x_old, y_old, x_new, y_new):
    nz = c_zyx.shape[0]
    out = np.empty((nz, len(y_new), len(x_new)), dtype=float)
    for k in range(nz):
        layer = fill_nans_2d(c_zyx[k, :, :])
        out[k, :, :] = interp_2d_regular(layer, x_old, y_old, x_new, y_new)
    return out


def extrapolate_columns_to_depth(c_zyx, z_in, z_out, c_fill=1480.0):
    ny, nx = c_zyx.shape[1], c_zyx.shape[2]
    out = np.full((len(z_out), ny, nx), c_fill, dtype=float)
    for j in range(ny):
        for i in range(nx):
            col = c_zyx[:, j, i]
            valid = np.isfinite(col)
            if valid.sum() == 0:
                continue
            if valid.sum() == 1:
                out[:, j, i] = col[valid][0]
                continue
            zv = z_in[valid]
            cv = col[valid]
            out[:, j, i] = np.interp(z_out, zv, cv)
            deep = z_out > zv[-1]
            slope = (cv[-1] - cv[-2]) / (zv[-1] - zv[-2])
            out[deep, j, i] = cv[-1] + slope * (z_out[deep] - zv[-1])
    return out


def read_gebco_bathymetry(path):
    with nc.Dataset(path, "r") as ds:
        lat = np.asarray(ds.variables["lat"][:], dtype=float)
        lon = np.asarray(ds.variables["lon"][:], dtype=float)
        elev = np.asarray(ds.variables["elevation"][:], dtype=float)
    bty_m = -elev
    return lat, lon, bty_m


def read_copernicus_temp_sal(temp_path, sal_path):
    def _to_float_with_nan(var_data):
        arr = np.ma.array(var_data)
        arr = np.ma.filled(arr, np.nan)
        arr = np.asarray(arr, dtype=float)
        arr[~np.isfinite(arr)] = np.nan
        arr[np.abs(arr) > 1.0e4] = np.nan
        return arr

    with nc.Dataset(temp_path, "r") as ds_t:
        t = _to_float_with_nan(ds_t.variables["thetao"][0, :, :, :])
        z = np.asarray(ds_t.variables["depth"][:], dtype=float)
        lat = np.asarray(ds_t.variables["latitude"][:], dtype=float)
        lon = np.asarray(ds_t.variables["longitude"][:], dtype=float)
    with nc.Dataset(sal_path, "r") as ds_s:
        s = _to_float_with_nan(ds_s.variables["so"][0, :, :, :])
    return z, lat, lon, t, s


def check_depth_conventions(ssp_z, bty_reg):
    if ssp_z.ndim != 1 or ssp_z.size < 2:
        raise ValueError("Copernicus depth axis must be a 1D array with at least 2 values.")
    if ssp_z[0] < 0.0 or np.any(np.diff(ssp_z) <= 0.0):
        raise ValueError("Copernicus depth axis must be +z downward (non-negative and strictly increasing).")
    if np.nanmax(bty_reg) <= 0.0:
        raise ValueError("Bathymetry has no +z underwater values after conversion from GEBCO elevation.")


def make_rtrs_grids(nx=120, ny=120, nz=120):
    bty_lat, bty_lon, bty = read_gebco_bathymetry(GEBCO_FILE)
    ssp_z, ssp_lat, ssp_lon, temp, sal = read_copernicus_temp_sal(COP_TEMP_FILE, COP_SAL_FILE)

    lat0 = 0.5 * (bty_lat.min() + bty_lat.max())
    lon0 = 0.5 * (bty_lon.min() + bty_lon.max())

    x_bty_raw, _ = latlon_to_xy_m(np.full_like(bty_lon, lat0), bty_lon, lat0, lon0)
    _, y_bty_raw = latlon_to_xy_m(bty_lat, np.full_like(bty_lat, lon0), lat0, lon0)

    x_reg = np.linspace(x_bty_raw.min(), x_bty_raw.max(), nx)
    y_reg = np.linspace(y_bty_raw.min(), y_bty_raw.max(), ny)
    bty_reg = interp_2d_nearest_regular(bty, x_bty_raw, y_bty_raw, x_reg, y_reg)
    if np.nanmin(bty) < 0.0 and np.nanmin(bty_reg) >= 0.0:
        print(
            "warning: source GEBCO contains altimetry (land, -z in RTRS convention), "
            "but none remains on the current regular grid; increase nx/ny to retain it."
        )
    check_depth_conventions(ssp_z, bty_reg)

    c = leroy_sound_speed(ssp_z, temp[None, :, :, :], sal[None, :, :, :], ssp_lat)[0]
    c[np.isnan(c)] = np.nan

    x_ssp_raw, _ = latlon_to_xy_m(np.full_like(ssp_lon, lat0), ssp_lon, lat0, lon0)
    _, y_ssp_raw = latlon_to_xy_m(ssp_lat, np.full_like(ssp_lat, lon0), lat0, lon0)

    c_reg_h = interp_3d_horizontal(c, x_ssp_raw, y_ssp_raw, x_reg, y_reg)
    z_reg = np.linspace(0.0, float(np.max(bty_reg)), nz)
    c_reg = extrapolate_columns_to_depth(c_reg_h, ssp_z, z_reg, c_fill=1480.0)

    return x_reg, y_reg, bty_reg, z_reg, c_reg


def make_env_json(x_reg, y_reg, bty_reg, z_reg, c_reg):
    c_flat = np.transpose(c_reg, (2, 1, 0)).flatten(order="C")
    bty_flat = np.transpose(np.asarray(bty_reg, dtype=float), (1, 0)).flatten(order="C")
    # print(bty_flat)

    # ray paths
    # return {
    #     "ssp": {
    #         "x_ssp_m": x_reg.tolist(),
    #         "y_ssp_m": y_reg.tolist(),
    #         "z_ssp_m": z_reg.tolist(),
    #         "c_m_s": c_flat.tolist(),
    #     },
    #     "bathymetry": {
    #         "x_bty_m": x_reg.tolist(),
    #         "y_bty_m": y_reg.tolist(),
    #         "z_bty_m": bty_flat.tolist(),
    #         "water_density_g_cm3": 1.0,
    #         "bottom_model": {
    #             "model": "elastic",
    #             "compressional_speed_m_s": 1700.0,
    #             "shear_speed_m_s": 400.0,
    #             "density_g_cm3": 1.6,
    #             "compressional_attenuation_db_per_wavelength": 0.5,
    #             "shear_attenuation_db_per_wavelength": 0.35,
    #         },
    #     },
    #     "source": {
    #         "position": [0.0, 0.0, 100.0],
    #         "freq_hz": [200.0],
    #         "launch_elev_deg": np.linspace(-15.0, 15.0, 5).tolist(),
    #         "launch_azim_deg": np.linspace(0.0, 359.0, 9).tolist(),
    #     },
    #     "receivers": {
    #         "config_type": "grid",
    #         # "x_rcvr_m": np.linspace(x_reg.min(), x_reg.max(), 260).tolist(),
    #         # "x_rcvr_m": np.linspace(-5000.0, 5000.0, 260).tolist(),
    #         "x_rcvr_m": [0.0],
    #         # "y_rcvr_m": np.linspace(y_reg.min(), y_reg.max(), 260).tolist(),
    #         "y_rcvr_m": np.linspace(0.0, 50000, 260).tolist(),
    #         # "z_rcvr_m": [100.0],
    #         # "z_rcvr_m": np.linspace(z_reg.min(), z_reg.max(), 260).tolist(),
    #         "z_rcvr_m": np.linspace(0.0, 2000.0, 260).tolist(),
    #     },
    #     "beam": {
    #         "step_m": 25.0,
    #         "max_steps": 100_000,
    #         "max_range_m": 75_000.0,
    #         "store_ray_paths": True,
    #         "integration_method": "rk2",
    #     },
    # }

    # x-y tl
    return {
        "ssp": {
            "x_ssp_m": x_reg.tolist(),
            "y_ssp_m": y_reg.tolist(),
            "z_ssp_m": z_reg.tolist(),
            "c_m_s": c_flat.tolist(),
        },
        "bathymetry": {
            "x_bty_m": x_reg.tolist(),
            "y_bty_m": y_reg.tolist(),
            "z_bty_m": bty_flat.tolist(),
            "water_density_g_cm3": 1.0,
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
            "position": [0.0, 0.0, 50.0],
            "freq_hz": [100.0],
            "launch_elev_deg": np.linspace(-15.0, 15.0, 15).tolist(),
            "launch_azim_deg": np.linspace(0.0, 359.0, 360).tolist(),
        },
        "receivers": {
            "config_type": "grid",
            # "x_rcvr_m": np.linspace(x_reg.min(), x_reg.max(), 260).tolist(),
            "x_rcvr_m": np.linspace(-25000.0, 25000.0, 400).tolist(),
            # "x_rcvr_m": [0.0],
            # "y_rcvr_m": np.linspace(y_reg.min(), y_reg.max(), 260).tolist(),
            "y_rcvr_m": np.linspace(-25000.0, 25000, 400).tolist(),
            "z_rcvr_m": [50.0],
            # "z_rcvr_m": np.linspace(z_reg.min(), z_reg.max(), 260).tolist(),
            # "z_rcvr_m": np.linspace(0.0, 2000.0, 260).tolist(),
        },
        "beam": {
            "step_m": 15.0,
            "max_steps": 50_000,
            "max_range_m": 25_000.0,
            "store_ray_paths": False,
            "show_progress": True,
            "integration_method": "rk2",
        },
    }


def plot_bathymetry(x_m, y_m, bty_m):
    x_km = x_m / 1000.0
    y_km = y_m / 1000.0
    xg, yg = np.meshgrid(x_km, y_km)
    fig, ax = plt.subplots(figsize=(7, 6))
    cmap = cmocean.cm.topo_r
    z = np.asarray(bty_m, dtype=float)
    zmin = float(np.nanmin(z))
    zmax = float(np.nanmax(z))
    eps = max(1e-9, 1e-6 * max(abs(zmin), abs(zmax), 1.0))
    vmin = zmin if zmin < 0.0 else -eps
    vmax = zmax if zmax > 0.0 else eps
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    cont = ax.contourf(xg, yg, z, levels=60, cmap=cmap, norm=norm)
    ax.scatter(0.0, 0.0, color="b", marker="o", s=40, label="Source")
    ax.set_xlabel("X (km)")
    ax.set_ylabel("Y (km)")
    ax.set_title("GEBCO Bathymetry")
    ax.legend(loc="upper right")
    cbar = fig.colorbar(cont, ax=ax)
    cbar.set_label("Bathymetry z (m, +down)")
    cbar.ax.invert_yaxis()
    fig.tight_layout()


def main():
    name = "test_real_data"
    jsonfile = f"examples/{name}.json"
    outfile = f"examples/{name}.out.json"

    x_reg, y_reg, bty_reg, z_reg, c_reg = make_rtrs_grids(nx=120, ny=120, nz=120)

    if LAUNCH_SSP_DASHBOARD:
        python_utils.launch_sound_speed_dashboard(
            x_reg,
            y_reg,
            z_reg,
            c_reg,
            bty_xy_m=bty_reg,
            port=SSP_DASHBOARD_PORT,
        )
        return

    env = make_env_json(x_reg, y_reg, bty_reg, z_reg, c_reg)

    if os.path.exists(outfile):
        os.remove(outfile)

    with open(jsonfile, "w") as f:
        json.dump(env, f, indent=2)

    t0 = time.time()
    os.system(f"cargo run --release {jsonfile}")
    print(f"completed in {time.time() - t0:.2f} seconds.")

    x_bty, y_bty, z_bty = python_utils.load_bty(outfile)
    freq, x_m, y_m, z_m, pressure = python_utils.load_cmpx_pressure(outfile)
    pressure = np.reshape(pressure, (len(freq), len(x_m), len(y_m), len(z_m)))
    tl = -20 * np.log10(np.maximum(np.abs(pressure), 1e-30))

    try: 
        rays = python_utils.load_rays(outfile)
        # python_utils.plot_rays_yz_bty(rays, x_bty, y_bty, z_bty)
        # python_utils.plot_rays_xz(rays)
        # python_utils.plot_rays_yz(rays)
        # python_utils.plot_rays_bty_3d(rays, x_bty, y_bty, z_bty)
        python_utils.plot_rays_bty_3d_plotly(rays, x_bty, y_bty, z_bty, show=True, z_scale=0.15, force_equal_xy=True)
    except:
        print("no ray paths stored, skipping ray plots")


    z_val = 10.0
    z_idx = (np.abs(z_m - z_val)).argmin()
    python_utils.plot_line_tl_y(tl[0, :, :, :], x_m, y_m, z_m, x_idx=len(x_m) // 2, z_idx=z_idx)
    python_utils.plot_line_tl_x(tl[0, :, :, :], x_m, y_m, z_m, y_idx=len(y_m) // 2, z_idx=z_idx)
    python_utils.plot_tl_yz(tl[0, :, :, :], x_m, y_m, z_m, x_idx=len(x_m) // 2)
    python_utils.plot_tl_xy(tl[0, :, :, :], x_m, y_m, z_m, z_idx=z_idx)
    plot_bathymetry(x_bty, y_bty, z_bty)
    plt.show()


if __name__ == "__main__":
    main()
