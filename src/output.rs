use anyhow::Result;
use ndarray::{Array3, Array4};
use num_complex::Complex32;
use serde::ser::{SerializeMap, SerializeSeq};
use serde::{Serialize, Serializer};
use std::fs::File;
use std::io::BufWriter;

use crate::influence::PressureField;
use crate::input::config::SimulationConfig;

struct RayPathsOut<'a> {
    ray_paths: &'a [Vec<[f32; 3]>],
}

impl Serialize for RayPathsOut<'_> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut map = serializer.serialize_map(Some(self.ray_paths.len()))?;
        for (i, path) in self.ray_paths.iter().enumerate() {
            map.serialize_entry(&format!("ray_{}", i), path)?;
        }
        SerializeMap::end(map)
    }
}

struct FlatArray3<'a> {
    array: &'a Array3<f32>,
}

impl Serialize for FlatArray3<'_> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let (nx, ny, nz) = self.array.dim();
        let mut map = serializer.serialize_map(Some(2))?;
        map.serialize_entry("shape", &[nx, ny, nz])?;

        struct Data<'a> {
            array: &'a Array3<f32>,
        }
        impl Serialize for Data<'_> {
            fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
            where
                S: Serializer,
            {
                let (nx, ny, nz) = self.array.dim();
                let mut seq = serializer.serialize_seq(Some(nx * ny * nz))?;
                for ix in 0..nx {
                    for iy in 0..ny {
                        for iz in 0..nz {
                            seq.serialize_element(&self.array[(ix, iy, iz)])?;
                        }
                    }
                }
                SerializeSeq::end(seq)
            }
        }

        map.serialize_entry("data", &Data { array: self.array })?;
        SerializeMap::end(map)
    }
}

struct FlatArray4Complex<'a> {
    array: &'a Array4<Complex32>,
    imag: bool,
}

impl Serialize for FlatArray4Complex<'_> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let (nfreq, nx, ny, nz) = self.array.dim();
        let mut map = serializer.serialize_map(Some(2))?;
        map.serialize_entry("shape", &[nfreq, nx, ny, nz])?;

        struct Data<'a> {
            array: &'a Array4<Complex32>,
            imag: bool,
        }
        impl Serialize for Data<'_> {
            fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
            where
                S: Serializer,
            {
                let (nfreq, nx, ny, nz) = self.array.dim();
                let mut seq = serializer.serialize_seq(Some(nfreq * nx * ny * nz))?;
                for ifreq in 0..nfreq {
                    for ix in 0..nx {
                        for iy in 0..ny {
                            for iz in 0..nz {
                                let v = self.array[(ifreq, ix, iy, iz)];
                                let x = if self.imag { v.im } else { v.re };
                                seq.serialize_element(&x)?;
                            }
                        }
                    }
                }
                SerializeSeq::end(seq)
            }
        }

        map.serialize_entry(
            "data",
            &Data {
                array: self.array,
                imag: self.imag,
            },
        )?;
        SerializeMap::end(map)
    }
}

struct PressureFieldOut<'a> {
    simulation_config: &'a SimulationConfig,
    pressure_field: &'a PressureField,
}

#[derive(Serialize)]
struct SrcOut<'a> {
    frequency_hz: &'a Vec<f32>,
    source_position_m: [f32; 3],
    launch_elev_deg: &'a Vec<f32>,
    launch_azim_deg: &'a Vec<f32>,
}

#[derive(Serialize)]
struct BtyOut<'a> {
    x_bty_m: &'a Vec<f32>,
    y_bty_m: &'a Vec<f32>,
    z_bty_m: &'a Vec<f32>,
}

impl Serialize for PressureFieldOut<'_> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut map = serializer.serialize_map(None)?;
        map.serialize_entry("frequency_hz", &self.simulation_config.source.freq_hz)?;

        if self.pressure_field.is_array {
            let recs = self
                .pressure_field
                .receiver_positions_m
                .as_ref()
                .expect("receiver_positions_m present for array mode");
            map.serialize_entry("receiver_positions_m", recs)?;
        } else {
            map.serialize_entry("x_m", &self.pressure_field.x_m)?;
            map.serialize_entry("y_m", &self.pressure_field.y_m)?;
            map.serialize_entry("z_m", &self.pressure_field.z_m)?;
        }

        map.serialize_entry(
            "delay_s",
            &FlatArray3 {
                array: &self.pressure_field.delay_s,
            },
        )?;
        map.serialize_entry(
            "amplitude",
            &FlatArray3 {
                array: &self.pressure_field.amplitude,
            },
        )?;
        map.serialize_entry(
            "pressure_re",
            &FlatArray4Complex {
                array: &self.pressure_field.pressure,
                imag: false,
            },
        )?;
        map.serialize_entry(
            "pressure_im",
            &FlatArray4Complex {
                array: &self.pressure_field.pressure,
                imag: true,
            },
        )?;
        SerializeMap::end(map)
    }
}

pub fn write_json(
    file_path: &str,
    simulation_config: &SimulationConfig,
    ray_paths: Option<&Vec<Vec<[f32; 3]>>>,
    pressure_field: &PressureField,
) -> Result<()> {
    let file = File::create(file_path)?;
    let writer = BufWriter::new(file);
    let mut serializer = serde_json::Serializer::new(writer);

    let mut root = serializer.serialize_map(None)?;
    root.serialize_entry(
        "src",
        &SrcOut {
            frequency_hz: &simulation_config.source.freq_hz,
            source_position_m: simulation_config.source.position,
            launch_elev_deg: &simulation_config.source.launch_elev_deg,
            launch_azim_deg: &simulation_config.source.launch_azim_deg,
        },
    )?;

    if let Some(ray_paths) = ray_paths {
        root.serialize_entry("ray_paths", &RayPathsOut { ray_paths })?;
    }

    root.serialize_entry(
        "bty",
        &BtyOut {
            x_bty_m: &simulation_config.bathymetry.x_bty_m,
            y_bty_m: &simulation_config.bathymetry.y_bty_m,
            z_bty_m: &simulation_config.bathymetry.z_bty_m,
        },
    )?;
    root.serialize_entry(
        "pressure_field",
        &PressureFieldOut {
            simulation_config,
            pressure_field,
        },
    )?;
    SerializeMap::end(root)?;

    Ok(())
}
