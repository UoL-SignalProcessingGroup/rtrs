pub mod config {
    use anyhow::{bail, Result};
    use serde::Deserialize;

    #[derive(Debug, Deserialize)]
    pub struct SoundSpeed {
        pub x_ssp_m: Vec<f32>, // x (m)
        pub y_ssp_m: Vec<f32>, // y (m)
        pub z_ssp_m: Vec<f32>, // z (m, positive down)
        pub c_m_s: Vec<f32>,   // sound speed (m/s)
    }

    impl SoundSpeed {
        fn validate(&mut self, errors: &mut Vec<String>) {
            if self.c_m_s.is_empty() {
                errors.push("ssp: c_m_s must not be empty".into());
            }
            if self.c_m_s.iter().any(|&c| c <= 0.0) {
                errors.push("ssp: all sound speeds must be positive".into());
            }
            for z in self.z_ssp_m.iter_mut().filter(|z| **z < 0.0) {
                eprintln!("WARNING ssp: z_ssp_m value {} is negative (positive down convention); correcting to {}", *z, z.abs());
                *z = z.abs();
            }
        }
    }

    #[derive(Debug, Deserialize)]
    pub struct Bathymetry {
        pub x_bty_m: Vec<f32>, // x (m)
        pub y_bty_m: Vec<f32>, // y (m)
        pub z_bty_m: Vec<f32>, // z (m, positive down)
        // not implemented:
        // pub density_g_cm3: f32, // bottom density (g/cm^3)
        // pub c_bty_m_s: f32, // bottom sound speed (m/s)
        // pub attenuation_db: f32, // bottom attenuation (dB/meter)
    }

    impl Bathymetry {
        fn validate(&self, errors: &mut Vec<String>) {
            if self.z_bty_m.is_empty() {
                errors.push("bathymetry: z_bty_m must not be empty".into());
            }
        }
    }

    #[derive(Debug, Deserialize)]
    pub struct Source {
        pub position: [f32; 3], // x, y, z (m)
        pub freq_hz: Vec<f32>,
        pub launch_elev_deg: Vec<f32>, // "alpha" angles
        pub launch_azim_deg: Vec<f32>, // "beta" angles
    }

    impl Source {
        fn validate(&mut self, errors: &mut Vec<String>) {
            if self.position[2] < 0.0 {
                eprintln!(
                    "WARNING source: z position {} is negative (positive down convention); correcting to {}",
                    self.position[2],
                    self.position[2].abs()
                );
                self.position[2] = self.position[2].abs();
            }
            if self.launch_elev_deg.iter().any(|&a| a < -90.0 || a > 90.0) {
                errors.push("source: launch elevation angles must be in [-90, 90] deg".into());
            }
        }
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

    impl Receivers {
        fn validate(&mut self, errors: &mut Vec<String>) {
            if self.config_type != "grid" && self.config_type != "array" {
                errors.push(format!(
                    "receivers: config_type must be \"grid\" or \"array\", got \"{}\"",
                    self.config_type
                ));
            }
            for z in self.z_rcvr_m.iter_mut().filter(|z| **z < 0.0) {
                eprintln!("WARNING receivers: z_rcvr_m value {} is negative (positive down convention); correcting to {}", *z, z.abs());
                *z = z.abs();
            }
        }
    }

    #[derive(Debug, Deserialize)]
    pub struct BeamSettings {
        pub step_m: f32,
        pub max_steps: usize,
        pub max_range_m: f32,
    }

    impl BeamSettings {
        fn validate(&self, errors: &mut Vec<String>) {
            if self.step_m <= 0.0 {
                errors.push("beam: step_m must be positive".into());
            }
            if self.max_steps == 0 {
                errors.push("beam: max_steps must be > 0".into());
            }
            if self.max_range_m <= 0.0 {
                errors.push("beam: max_range_m must be positive".into());
            }
        }
    }

    #[derive(Debug, Deserialize)]
    pub struct SimulationConfig {
        pub ssp: SoundSpeed,
        pub bathymetry: Bathymetry,
        pub source: Source,
        pub receivers: Receivers,
        pub beam: BeamSettings,
    }

    impl SimulationConfig {
        pub fn validate(&mut self) -> Result<()> {
            let mut errors: Vec<String> = Vec::new();
            self.ssp.validate(&mut errors);
            self.bathymetry.validate(&mut errors);
            self.source.validate(&mut errors);
            self.receivers.validate(&mut errors);
            self.beam.validate(&mut errors);
            if errors.is_empty() {
                Ok(())
            } else {
                bail!("Invalid simulation config:\n  - {}", errors.join("\n  - "))
            }
        }
    }
}