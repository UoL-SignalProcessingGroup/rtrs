import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmocean
import plotly.graph_objects as go


def launch_sound_speed_dashboard(
    x_ssp_m,
    y_ssp_m,
    z_ssp_m,
    c_zyx,
    bty_xy_m=None,
    host='127.0.0.1',
    port=8050,
    debug=False,
):
    from dash import Dash, dcc, html, Input, Output, State

    x_ssp_m = np.asarray(x_ssp_m, dtype=float)
    y_ssp_m = np.asarray(y_ssp_m, dtype=float)
    z_ssp_m = np.asarray(z_ssp_m, dtype=float)
    c_zyx = np.asarray(c_zyx, dtype=float)

    if c_zyx.ndim != 3:
        raise ValueError('c_zyx must have shape (nz, ny, nx).')

    nz, ny, nx = c_zyx.shape
    if len(x_ssp_m) != nx or len(y_ssp_m) != ny or len(z_ssp_m) != nz:
        raise ValueError('Grid axes lengths must match c_zyx dimensions (nz, ny, nx).')

    if bty_xy_m is not None:
        bty_xy_m = np.asarray(bty_xy_m, dtype=float)
        if bty_xy_m.shape != (ny, nx):
            raise ValueError('bty_xy_m must have shape (ny, nx), matching c_zyx[0].')

    dark_layout = dict(
        paper_bgcolor='#1a1a2e',
        plot_bgcolor='#16213e',
        font=dict(color='#eaeaea', size=11),
        margin=dict(l=50, r=20, t=42, b=42),
    )

    def _mpl_cmap_to_plotly(cmap, n=256):
        vals = np.linspace(0.0, 1.0, n)
        rgba = cmap(vals)
        return [
            [float(v), f'rgb({int(r * 255)},{int(g * 255)},{int(b * 255)})']
            for v, (r, g, b, _) in zip(vals, rgba)
        ]

    def _two_slope_colorscale(cmap, zmin, zmax, vcenter=0.0, n=256):
        if not np.isfinite(zmin) or not np.isfinite(zmax) or zmax <= zmin:
            return _mpl_cmap_to_plotly(cmap, n=n)

        t0 = float(np.clip((vcenter - zmin) / (zmax - zmin), 0.001, 0.999))
        n_lo = max(2, int(n * t0))
        n_hi = max(2, n - n_lo)
        cs = []

        for i, c in enumerate(np.linspace(0.0, 0.5, n_lo)):
            pos = t0 * i / (n_lo - 1)
            r, g, b, _ = cmap(c)
            cs.append([round(float(pos), 6), f'rgb({int(r * 255)},{int(g * 255)},{int(b * 255)})'])

        for i, c in enumerate(np.linspace(0.5, 1.0, n_hi)):
            if i == 0:
                continue
            pos = t0 + (1.0 - t0) * i / (n_hi - 1)
            r, g, b, _ = cmap(c)
            cs.append([round(float(pos), 6), f'rgb({int(r * 255)},{int(g * 255)},{int(b * 255)})'])

        return cs

    colorscales = {
        'viridis': 'Viridis',
        'thermal': _mpl_cmap_to_plotly(cmocean.cm.thermal),
        'bathymetry': _mpl_cmap_to_plotly(cmocean.cm.topo_r),
    }

    def _nearest_index(values, target):
        return int(np.argmin(np.abs(values - float(target))))

    ix0 = nx // 2
    iy0 = ny // 2
    cmin = float(np.nanmin(c_zyx))
    cmax = float(np.nanmax(c_zyx))

    nmarks = min(8, nz)
    mark_indices = np.linspace(0, nz - 1, nmarks, dtype=int)
    marks = {int(i): f'{z_ssp_m[i]:.0f} m' for i in mark_indices}

    def _map_figure(depth_idx, ix, iy, cmap_key, selector_mode):
        depth_idx = int(np.clip(depth_idx, 0, nz - 1))
        if selector_mode == 'bathymetry' and bty_xy_m is not None:
            z_map = bty_xy_m
            zmin = float(np.nanmin(z_map))
            zmax = float(np.nanmax(z_map))
            colorscale = _two_slope_colorscale(cmocean.cm.topo_r, zmin, zmax, vcenter=0.0)
            title = 'Bathymetry Selector'
            hovertemplate = 'x=%{x:.0f} m<br>y=%{y:.0f} m<br>z_bty=%{z:.2f} m<extra></extra>'
            colorbar_title = 'z_bty (m)'
        else:
            z_map = c_zyx[depth_idx, :, :]
            zmin = cmin
            zmax = cmax
            colorscale = colorscales[cmap_key]
            title = f'SSP X-Y Slice Selector | depth = {z_ssp_m[depth_idx]:.1f} m'
            hovertemplate = 'x=%{x:.0f} m<br>y=%{y:.0f} m<br>c=%{z:.2f} m/s<extra></extra>'
            colorbar_title = 'c (m/s)'

        fig = go.Figure()
        fig.add_trace(
            go.Heatmap(
                x=x_ssp_m,
                y=y_ssp_m,
                z=z_map,
                colorscale=colorscale,
                zmin=zmin,
                zmax=zmax,
                colorbar=dict(title=dict(text=colorbar_title, side='right'), thickness=12, tickfont=dict(size=9)),
                hovertemplate=hovertemplate,
            )
        )
        fig.add_trace(
            go.Scatter(
                x=[x_ssp_m[ix]],
                y=[y_ssp_m[iy]],
                mode='markers',
                marker=dict(symbol='star', color='yellow', size=13, line=dict(color='black', width=1)),
                name='Selected',
                hovertemplate='Selected<br>x=%{x:.1f} m<br>y=%{y:.1f} m<extra></extra>',
            )
        )
        fig.update_layout(
            title=dict(text=title, x=0.5),
            xaxis=dict(title='x (m)', gridcolor='#2a2a4a'),
            yaxis=dict(title='y (m)', gridcolor='#2a2a4a', scaleanchor='x', scaleratio=1),
            legend=dict(x=0.01, y=0.99, bgcolor='rgba(0,0,0,0.5)', font=dict(size=10)),
            **dark_layout,
        )
        return fig

    def _profile_figure(ix, iy):
        profile = c_zyx[:, iy, ix]
        finite = profile[np.isfinite(profile)]
        if finite.size:
            x_min = float(np.nanmin(finite)) - 5.0
            x_max = float(np.nanmax(finite)) + 5.0
        else:
            x_min = 1450.0
            x_max = 1550.0

        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                x=profile,
                y=z_ssp_m,
                mode='lines',
                line=dict(width=2, color='#4fc3f7'),
                hovertemplate='c=%{x:.2f} m/s<br>z=%{y:.1f} m<extra></extra>',
                name='c(z)',
            )
        )

        if bty_xy_m is not None:
            bty_depth = float(bty_xy_m[iy, ix])
            fig.add_trace(
                go.Scatter(
                    x=[x_min, x_max],
                    y=[bty_depth, bty_depth],
                    mode='lines',
                    line=dict(color='saddlebrown', width=1.5, dash='dot'),
                    showlegend=False,
                    hovertemplate=f'Bottom: {bty_depth:.1f} m<extra></extra>',
                    name='Bathymetry',
                )
            )

            z_max = max(float(np.nanmax(z_ssp_m)), bty_depth)
            yaxis = dict(title='depth (m)', autorange='reversed', gridcolor='#2a2a4a', range=[z_max, float(np.nanmin(z_ssp_m))])
        else:
            yaxis = dict(title='depth (m)', autorange='reversed', gridcolor='#2a2a4a')

        fig.update_layout(
            title=dict(text='Sound Speed Profile', x=0.5),
            xaxis=dict(title='c (m/s)', gridcolor='#2a2a4a', range=[x_min, x_max]),
            yaxis=yaxis,
            showlegend=False,
            **dark_layout,
        )
        return fig

    app = Dash(__name__)
    app.title = 'Sound Speed Dashboard'

    app.layout = html.Div(
        style={'fontFamily': 'Arial, sans-serif', 'backgroundColor': '#1a1a2e', 'color': '#eaeaea', 'padding': '12px'},
        children=[
            html.H2('Sound Speed Dashboard', style={'textAlign': 'center', 'color': '#e0e0ff', 'marginBottom': '6px'}),
            html.Div(
                id='info-bar',
                style={'textAlign': 'center', 'marginBottom': '10px', 'fontSize': '13px', 'color': '#aaa'},
                children='Click any point on the SSP selector to inspect c(z) at that location.',
            ),
            dcc.Store(id='selected-indices', data={'ix': ix0, 'iy': iy0}),
            html.Div(
                style={'padding': '6px 10px', 'backgroundColor': '#16213e', 'borderRadius': '4px', 'marginBottom': '10px'},
                children=[
                    html.Div(
                        style={'display': 'flex', 'alignItems': 'center', 'gap': '10px', 'marginBottom': '6px'},
                        children=[
                            html.Span('Selector data:', style={'color': '#aaa', 'fontSize': '12px'}),
                            dcc.Tabs(
                                id='selector-mode',
                                value='ssp_xy',
                                children=[
                                    dcc.Tab(
                                        label='SSP slices in x-y',
                                        value='ssp_xy',
                                        style={'backgroundColor': '#16213e', 'color': '#aaa', 'borderColor': '#2a2a4a'},
                                        selected_style={'backgroundColor': '#2a2a4a', 'color': '#eaeaea', 'borderColor': '#4a4a8a'},
                                    ),
                                    dcc.Tab(
                                        label='Bathymetry map',
                                        value='bathymetry',
                                        style={'backgroundColor': '#16213e', 'color': '#aaa', 'borderColor': '#2a2a4a'},
                                        selected_style={'backgroundColor': '#2a2a4a', 'color': '#eaeaea', 'borderColor': '#4a4a8a'},
                                    ),
                                ],
                                style={'borderBottom': '1px solid #2a2a4a', 'width': '360px'},
                            ),
                            html.Span('Selector colormap:', style={'color': '#aaa', 'fontSize': '12px'}),
                            dcc.Dropdown(
                                id='cmap-selector',
                                options=[
                                    {'label': 'Viridis', 'value': 'viridis'},
                                    {'label': 'Thermal', 'value': 'thermal'},
                                    {'label': 'Bathymetry (topo_r)', 'value': 'bathymetry'},
                                ],
                                value='bathymetry',
                                clearable=False,
                                style={'width': '260px', 'color': '#111'},
                            ),
                        ],
                    ),
                    html.Span('Depth (m):', style={'color': '#aaa', 'fontSize': '11px', 'display': 'block', 'marginBottom': '2px'}),
                    dcc.Slider(
                        id='depth-slider',
                        min=0,
                        max=nz - 1,
                        step=1,
                        value=0,
                        marks=marks,
                        tooltip={'placement': 'bottom', 'always_visible': False},
                        updatemode='drag',
                    ),
                ],
            ),
            html.Div(
                style={'display': 'flex', 'gap': '10px'},
                children=[
                    dcc.Graph(id='ssp-map', style={'flex': '2', 'minWidth': '0', 'height': '62vh'}, config={'scrollZoom': True}),
                    dcc.Graph(id='ssp-profile', style={'flex': '1', 'minWidth': '0', 'height': '62vh'}),
                ],
            ),
            html.Div(id='selected-point-label', style={'textAlign': 'center', 'marginTop': '10px', 'fontSize': '13px', 'color': '#aaa'}),
        ],
    )

    @app.callback(
        Output('selected-indices', 'data'),
        Input('ssp-map', 'clickData'),
        State('selected-indices', 'data'),
    )
    def _update_selected_indices(click_data, selected_data):
        if not click_data or 'points' not in click_data or len(click_data['points']) == 0:
            return selected_data
        point = click_data['points'][0]
        if 'x' not in point or 'y' not in point:
            return selected_data
        ix = _nearest_index(x_ssp_m, point['x'])
        iy = _nearest_index(y_ssp_m, point['y'])
        return {'ix': ix, 'iy': iy}

    @app.callback(
        Output('ssp-map', 'figure'),
        Input('depth-slider', 'value'),
        Input('selected-indices', 'data'),
        Input('cmap-selector', 'value'),
        Input('selector-mode', 'value'),
    )
    def _update_map(depth_idx, selected_data, cmap_key, selector_mode):
        ix = int(selected_data['ix'])
        iy = int(selected_data['iy'])
        cmap_key = cmap_key if cmap_key in colorscales else 'bathymetry'
        selector_mode = selector_mode if selector_mode in {'ssp_xy', 'bathymetry'} else 'ssp_xy'
        return _map_figure(depth_idx, ix, iy, cmap_key, selector_mode)

    @app.callback(
        Output('depth-slider', 'disabled'),
        Output('cmap-selector', 'disabled'),
        Input('selector-mode', 'value'),
    )
    def _toggle_controls(selector_mode):
        is_bathymetry = selector_mode == 'bathymetry'
        return is_bathymetry, is_bathymetry

    @app.callback(
        Output('ssp-profile', 'figure'),
        Output('selected-point-label', 'children'),
        Input('selected-indices', 'data'),
    )
    def _update_profile(selected_data):
        ix = int(selected_data['ix'])
        iy = int(selected_data['iy'])
        label = f'Selected point: x={x_ssp_m[ix]:.1f} m, y={y_ssp_m[iy]:.1f} m'
        if bty_xy_m is not None:
            label = f'{label}, z_bty={bty_xy_m[iy, ix]:.1f} m'
        return _profile_figure(ix, iy), label

    print(f'\n  Sound speed dashboard running → http://{host}:{port}/\n')
    app.run(host=host, port=int(port), debug=bool(debug))


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

def plot_rays_bty_3d(rays, x_bty, y_bty, z_bty, figsize=(8, 6), z_scale=1.0):
    if z_scale <= 0:
        raise ValueError("z_scale must be > 0")

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')
    # Plot rays
    for ray in rays:
        x = ray[:,0]
        y = ray[:,1]
        z = ray[:,2]
        ax.plot(x, y, z, lw=0.8)
    # Plot bathymetry surface
    X, Y = np.meshgrid(x_bty, y_bty)
    z = np.asarray(z_bty, dtype=float)
    zmin = float(np.nanmin(z))
    zmax = float(np.nanmax(z))
    eps = max(1e-9, 1e-6 * max(abs(zmin), abs(zmax), 1.0))
    vmin = zmin if zmin < 0.0 else -eps
    vmax = zmax if zmax > 0.0 else eps
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    cmap = cmocean.cm.topo_r
    ax.plot_surface(X, Y, z, cmap=cmap, norm=norm, alpha=0.5, linewidth=0, antialiased=False)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('depth (m)')
    ax.set_title('Ray paths (3D) with Bathymetry')
    ax.set_xlim(np.min(x_bty), np.max(x_bty))
    ax.set_ylim(np.min(y_bty), np.max(y_bty))
    ax.set_box_aspect((1.0, 1.0, float(z_scale)))
    ax.invert_zaxis()
    plt.tight_layout()

def plot_rays_bty_3d_plotly(
    rays,
    x_bty,
    y_bty,
    z_bty,
    show=True,
    figsize=(1000, 700),
    z_scale=1.0,
    force_equal_xy=False,
):
    if z_scale <= 0:
        raise ValueError("z_scale must be > 0")

    X, Y = np.meshgrid(x_bty, y_bty)

    z = np.asarray(z_bty, dtype=float)
    zmin = float(np.nanmin(z))
    zmax = float(np.nanmax(z))
    eps = max(1e-9, 1e-6 * max(abs(zmin), abs(zmax), 1.0))
    vmin = zmin if zmin < 0.0 else -eps
    vmax = zmax if zmax > 0.0 else eps
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    z_norm = norm(z)

    cmap = cmocean.cm.topo_r
    colorscale = [
        [i / 255.0, f"rgb({int(r * 255)}, {int(g * 255)}, {int(b * 255)})"]
        for i, (r, g, b, _) in enumerate(cmap(np.linspace(0.0, 1.0, 256)))
    ]

    x_min = float(np.min(x_bty))
    x_max = float(np.max(x_bty))
    y_min = float(np.min(y_bty))
    y_max = float(np.max(y_bty))
    if force_equal_xy:
        x_mid = 0.5 * (x_min + x_max)
        y_mid = 0.5 * (y_min + y_max)
        xy_span = max(x_max - x_min, y_max - y_min)
        half_span = 0.5 * xy_span
        x_min, x_max = x_mid - half_span, x_mid + half_span
        y_min, y_max = y_mid - half_span, y_mid + half_span

    fig = go.Figure()

    fig.add_trace(
        go.Surface(
            x=X,
            y=Y,
            z=z,
            surfacecolor=z_norm,
            colorscale=colorscale,
            cmin=0.0,
            cmax=1.0,
            opacity=1.0,
            showscale=False,
            name='Bathymetry',
        )
    )

    for idx, ray in enumerate(rays):
        x = ray[:, 0]
        y = ray[:, 1]
        z = ray[:, 2]
        fig.add_trace(
            go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode='lines',
                line=dict(width=3),
                name=f'Ray {idx + 1}',
                showlegend=False,
            )
        )

    fig.update_layout(
        title='Ray paths (3D) with Bathymetry',
        width=figsize[0],
        height=figsize[1],
        scene=dict(
            aspectmode='manual',
            aspectratio=dict(x=1.0, y=1.0, z=float(z_scale)),
            xaxis=dict(title='x (m)', range=[x_min, x_max]),
            yaxis=dict(title='y (m)', range=[y_min, y_max]),
            zaxis=dict(title='depth (m)', autorange='reversed'),
        ),
        margin=dict(l=0, r=0, b=0, t=40),
    )

    if show:
        fig.show()

    return fig

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

def plot_tl_along_array(tl):
    plt.figure()
    plt.plot(tl)
    plt.xlabel('Array index')
    plt.ylabel('Transmission Loss (dB)')
    plt.title('Transmission Loss along array')
    plt.gca().invert_yaxis()
    plt.grid()
    plt.tight_layout()

def plot_pressure_along_array(pressure):
    plt.figure()
    plt.plot(np.real(pressure))
    plt.xlabel('Array index')
    plt.ylabel('Pressure (real)')
    plt.title('Pressure along array')
    plt.grid()
    plt.tight_layout()

def plot_phase_along_array(pressure):
    plt.figure()
    plt.plot(np.angle(pressure))
    plt.xlabel('Array index')
    plt.ylabel('Phase (radians)')
    plt.title('Phase along array')
    plt.grid()
    plt.tight_layout()

def plot_array_geometry(x_m, y_m, z_m=None, source_pos=None, receiver_positions=None):
    plt.figure()

    # plt.plot(x_m, z_m, 'k-', lw=1.5, label='Array depth profile')
    plt.plot(x_m, y_m, 'k-*', lw=1.5, label='Array horizontal profile')

    if source_pos is not None:
        plt.plot(source_pos[0], source_pos[1], 'r*', markersize=12, label='Source')
    if receiver_positions is not None:
        plt.scatter(receiver_positions[:, 0], receiver_positions[:, 1], c='b', marker='o', label='Receivers')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Array Geometry (top-down view)')
    plt.axis('equal')
    plt.grid()
    plt.legend()
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

def plot_tl_xy(tl, x_m, y_m, z_m, z_idx, vmin=50, vmax=150):
    plt.figure()
    X, Y = np.meshgrid(x_m, y_m)
    plt.pcolormesh(X, Y, tl[:, :, z_idx].T, shading='auto', cmap='jet_r', vmin=vmin, vmax=vmax)
    plt.colorbar(label='Transmission Loss (dB)')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title(f'Transmission Loss (x-y plane) at z={z_m[z_idx]:.1f} m')
    # plt.gca().invert_yaxis()
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