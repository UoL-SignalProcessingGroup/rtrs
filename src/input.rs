pub mod config {
    use ndarray::Array3;
    use serde::Deserialize;

    #[derive(Debug, Deserialize)]
    pub struct SoundSpeed {
        pub interp_type: String, // "tabulated" (for now)
        pub x_ssp_m: Vec<f64>, // x (m)
        pub y_ssp_m: Vec<f64>, // y (m)
        pub z_ssp_m: Vec<f64>, // z (m, positive down)
        pub c_m_s: Vec<f64>, // sound speed (m/s)
    }

    #[derive(Debug, Deserialize)]
    pub struct Bathymetry {
        pub kind: String, // "flat" or "3D"
        pub flat_depth_m: Option<f64>, // for "flat" kind
        pub x_bty_m: Option<Vec<f64>>, // x (m)
        pub y_bty_m: Option<Vec<f64>>, // y (m)
        pub z_bty_m: Option<Vec<f64>>, // z (m, positive down)
        pub density_g_cm3: Vec<f64>, // bottom density (g/cm^3)
        pub c_bty_m_s: Vec<f64>, // bottom sound speed (m/s)
        pub attenuation_db_per_wavelength: Vec<f64>, // bottom attenuation (dB/wavelength)
    }

    #[derive(Debug, Deserialize)]
    pub struct Source {
        pub position: [f64; 3], // x, y, z (m)
        pub freq_hz: f64,
        pub launch_elev_deg: Vec<f64>, // "alpha" angles
        pub launch_azim_deg: Vec<f64>,  // "beta" angles
    }

    #[derive(Debug, Deserialize)]
    pub struct Receivers {
        pub kind: String, // "clyindrical" or "cartesian"
        // Cylindrical (r, theta, z) coordinates
        pub ranges_m: Option<Vec<f64>>,
        pub bearings_deg: Option<Vec<f64>>, // "theta" angles
        pub depths_m: Option<Vec<f64>>,
        // Cartesian (x, y, z) coordinates
        pub x_rcvr_m: Option<Vec<f64>>,
        pub y_rcvr_m: Option<Vec<f64>>,
        pub z_rcvr_m: Option<Vec<f64>>,
    }

    #[derive(Debug, Deserialize)]
    pub struct BeamSettings {
        pub step_m: f64,
        pub beam_type: String, 
    }

    #[derive(Debug, Deserialize)]
    pub struct SimulationConfig {
        pub ssp: SoundSpeed,
        pub bathymetry: Bathymetry,
        pub source: Source,
        pub receivers: Receivers,
        pub beam: BeamSettings,
    }
}




/*
Example JSON template for SimulationConfig:

{
    "ssp": {
        "interp_type": "tabulated",
        "x_ssp_m": [0.0, 5000.0],
        "y_ssp_m": [0.0, 5000.0],
        "z_ssp_m": [0.0, 100.0],
        // flattened 3D array: c_m_s[x][y][z] as a flat Vec<f64>
        "c_m_s": [1500.0, 1520.0, 1510.0, 1530.0, 1490.0, 1510.0, 1500.0, 1520.0]
    }
    },
    "bathymetry": {
        "kind": "flat",
        "flat_depth_m": 100.0,
        "x_bty_m": null,
        "y_bty_m": null,
        "z_bty_m": null,
        "density_g_cm3": [2.0],
        "c_bty_m_s": [1600.0],
        "attenuation_db_per_wavelength": [0.5]
    },
    "source": {
        "position": [0.0, 0.0, 10.0],
        "freq_hz": 1000.0,
        "launch_elev_deg": [10.0, 20.0, 30.0],
        "launch_azim_deg": [0.0, 90.0, 180.0]
    },
    "receivers": {
        "kind": "cartesian",
        "ranges_m": null,
        "bearings_deg": null,
        "depths_m": null,
        "x_rcvr_m": [100.0, 200.0],
        "y_rcvr_m": [0.0, 0.0],
        "z_rcvr_m": [20.0, 30.0]
    },
    "beam": {
        "step_m": 1.0,
        "beam_type": "ray"
    }
}
*/