pub mod config {
    use anyhow::{Result, bail};
    use serde::Deserialize;

    fn default_none_f32() -> Option<f32> {
        None
    }
    fn default_false() -> bool {
        false
    }
    fn default_integration_method() -> IntegrationMethod {
        IntegrationMethod::Euler
    }

    #[derive(Debug, Deserialize)]
    pub struct SoundSpeed {
        pub x_ssp_m: Vec<f32>, // x (m)
        pub y_ssp_m: Vec<f32>, // y (m)
        pub z_ssp_m: Vec<f32>, // z (m, positive down)
        pub c_m_s: Vec<f32>,   // sound speed (m/s)
    }

    impl SoundSpeed {
        fn validate(&mut self, errors: &mut Vec<String>, warnings: &mut Vec<String>) {
            if self.c_m_s.is_empty() {
                errors.push("ssp: c_m_s must not be empty".into());
            }
            if self.c_m_s.iter().any(|&c| c <= 0.0) {
                errors.push("ssp: all sound speeds must be positive".into());
            }
            for z in self.z_ssp_m.iter_mut().filter(|z| **z < 0.0) {
                warnings.push(format!("ssp: z_ssp_m value {} is negative (positive down convention); correcting to {}", *z, z.abs()));
                *z = z.abs();
            }
        }
    }

    #[derive(Debug, Deserialize)]
    pub struct Bathymetry {
        pub x_bty_m: Vec<f32>, // x (m)
        pub y_bty_m: Vec<f32>, // y (m)
        pub z_bty_m: Vec<f32>, // z (m, positive down)
        #[serde(default = "default_none_f32")]
        pub water_density_g_cm3: Option<f32>,
        #[serde(default = "default_bottom_boundary_model")]
        pub bottom_model: BottomBoundaryModel,
    }

    fn default_bottom_boundary_model() -> BottomBoundaryModel {
        BottomBoundaryModel::Rigid
    }

    #[derive(Debug, Deserialize, Clone)]
    #[serde(rename_all = "snake_case", tag = "model")]
    pub enum BottomBoundaryModel {
        Rigid,
        Acoustic {
            compressional_speed_m_s: f32,
            density_g_cm3: f32,
            #[serde(default = "default_none_f32")]
            compressional_attenuation_db_per_wavelength: Option<f32>,
        },
        Elastic {
            compressional_speed_m_s: f32,
            shear_speed_m_s: f32,
            density_g_cm3: f32,
            #[serde(default = "default_none_f32")]
            compressional_attenuation_db_per_wavelength: Option<f32>,
            #[serde(default = "default_none_f32")]
            shear_attenuation_db_per_wavelength: Option<f32>,
        },
    }

    impl BottomBoundaryModel {
        fn validate(&self, errors: &mut Vec<String>) {
            match self {
                BottomBoundaryModel::Rigid => {}
                BottomBoundaryModel::Acoustic {
                    compressional_speed_m_s,
                    density_g_cm3,
                    compressional_attenuation_db_per_wavelength,
                } => {
                    if *compressional_speed_m_s <= 0.0 {
                        errors.push("bathymetry.bottom_model: acoustic compressional_speed_m_s must be positive".into());
                    }
                    if *density_g_cm3 <= 0.0 {
                        errors.push(
                            "bathymetry.bottom_model: acoustic density_g_cm3 must be positive"
                                .into(),
                        );
                    }
                    if let Some(value) = compressional_attenuation_db_per_wavelength {
                        if *value < 0.0 {
                            errors.push("bathymetry.bottom_model: acoustic compressional_attenuation_db_per_wavelength must be >= 0".into());
                        }
                    }
                }
                BottomBoundaryModel::Elastic {
                    compressional_speed_m_s,
                    shear_speed_m_s,
                    density_g_cm3,
                    compressional_attenuation_db_per_wavelength,
                    shear_attenuation_db_per_wavelength,
                } => {
                    if *compressional_speed_m_s <= 0.0 {
                        errors.push("bathymetry.bottom_model: elastic compressional_speed_m_s must be positive".into());
                    }
                    if *shear_speed_m_s <= 0.0 {
                        errors.push(
                            "bathymetry.bottom_model: elastic shear_speed_m_s must be positive"
                                .into(),
                        );
                    }
                    if *density_g_cm3 <= 0.0 {
                        errors.push(
                            "bathymetry.bottom_model: elastic density_g_cm3 must be positive"
                                .into(),
                        );
                    }
                    if let Some(value) = compressional_attenuation_db_per_wavelength {
                        if *value < 0.0 {
                            errors.push("bathymetry.bottom_model: elastic compressional_attenuation_db_per_wavelength must be >= 0".into());
                        }
                    }
                    if let Some(value) = shear_attenuation_db_per_wavelength {
                        if *value < 0.0 {
                            errors.push("bathymetry.bottom_model: elastic shear_attenuation_db_per_wavelength must be >= 0".into());
                        }
                    }
                }
            }
        }
    }

    impl Bathymetry {
        fn validate(&self, errors: &mut Vec<String>, _warnings: &mut Vec<String>) {
            if self.z_bty_m.is_empty() {
                errors.push("bathymetry: z_bty_m must not be empty".into());
            }
            if let Some(value) = self.water_density_g_cm3 {
                if value <= 0.0 {
                    errors.push(
                        "bathymetry: water_density_g_cm3 must be positive when provided".into(),
                    );
                }
            }
            self.bottom_model.validate(errors);
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
        fn validate(&mut self, errors: &mut Vec<String>, warnings: &mut Vec<String>) {
            if self.position[2] < 0.0 {
                warnings.push(format!(
                    "source: z position {} is negative (positive down convention); correcting to {}",
                    self.position[2],
                    self.position[2].abs()
                ));
                self.position[2] = self.position[2].abs();
            }
            // if self.freq_hz.iter().any(|&f| f < 0.0 ) {
            //     errors.push("source: frequencies in freq_hz must be non-negative".into());
            // }
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
        fn validate(&mut self, errors: &mut Vec<String>, warnings: &mut Vec<String>) {
            if self.config_type != "grid" && self.config_type != "array" {
                errors.push(format!(
                    "receivers: config_type must be \"grid\" or \"array\", got \"{}\"",
                    self.config_type
                ));
            }
            for z in self.z_rcvr_m.iter_mut().filter(|z| **z < 0.0) {
                warnings.push(format!("receivers: z_rcvr_m value {} is negative (positive down convention); correcting to {}", *z, z.abs()));
                *z = z.abs();
            }
        }
    }

    #[derive(Debug, Deserialize)]
    #[serde(rename_all = "snake_case")]
    pub enum IntegrationMethod {
        Euler,
        Rk2,
    }

    #[derive(Debug, Deserialize)]
    pub struct BeamSettings {
        pub step_m: f32,
        pub max_steps: usize,
        pub max_range_m: f32,
        #[serde(default = "default_false")]
        pub store_ray_paths: bool,
        #[serde(default = "default_false")]
        pub show_progress: bool,
        #[serde(default = "default_integration_method")]
        pub integration_method: IntegrationMethod,
    }

    impl BeamSettings {
        fn validate(&self, errors: &mut Vec<String>, _warnings: &mut Vec<String>) {
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
        pub fn validate(&mut self) -> Result<Vec<String>> {
            let mut errors: Vec<String> = Vec::new();
            let mut warnings: Vec<String> = Vec::new();
            self.ssp.validate(&mut errors, &mut warnings);
            self.bathymetry.validate(&mut errors, &mut warnings);
            self.source.validate(&mut errors, &mut warnings);
            self.receivers.validate(&mut errors, &mut warnings);
            self.beam.validate(&mut errors, &mut warnings);

            if errors.is_empty() {
                Ok(warnings)
            } else {
                bail!("Invalid simulation config:\n  - {}", errors.join("\n  - "))
            }
        }
    }
}
