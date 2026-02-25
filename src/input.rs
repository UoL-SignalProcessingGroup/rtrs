pub mod config {
    use serde::Deserialize;

    #[derive(Debug, Deserialize)]
    pub struct SoundSpeed {
        pub x_ssp_m: Vec<f32>, // x (m)
        pub y_ssp_m: Vec<f32>, // y (m)
        pub z_ssp_m: Vec<f32>, // z (m, positive down)
        pub c_m_s: Vec<f32>, // sound speed (m/s)
    }

    #[derive(Debug, Deserialize)]
    pub struct Bathymetry {
        pub x_bty_m: Vec<f32>, // x (m)
        pub y_bty_m: Vec<f32>, // y (m)
        pub z_bty_m: Vec<f32>, // z (m, positive down)
        // pub density_g_cm3: f32, // bottom density (g/cm^3)
        // pub c_bty_m_s: f32, // bottom sound speed (m/s)
        // pub attenuation_db: f32, // bottom attenuation (dB/meter)
    }

    #[derive(Debug, Deserialize)]
    pub struct Source {
        pub position: [f32; 3], // x, y, z (m)
        pub freq_hz: Vec<f32>,
        pub launch_elev_deg: Vec<f32>, // "alpha" angles
        pub launch_azim_deg: Vec<f32>,  // "beta" angles
    }

    #[derive(Debug, Deserialize)]
    pub struct Receivers {
        pub config_type: String, // "grid" or "array"
        // if the config_type is "grid" then the receivers are bellow are the axis vecotors
        // if "array" then these are the explicit receiver coordinate positions
        pub x_rcvr_m: Vec<f32>,
        pub y_rcvr_m: Vec<f32>,
        pub z_rcvr_m: Vec<f32>,
    }

    #[derive(Debug, Deserialize)]
    pub struct BeamSettings {
        pub step_m: f32,
        pub max_steps: usize,
        pub max_range_m: f32,
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