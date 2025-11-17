use anyhow::Result;
use hdf5_metno as hdf5;
use ndarray::Array4;
// use num_complex

use crate::{input::config::SimulationConfig};
use crate::influence::PressureField;

// pressure_field: Array3<[f32; 2]> // alt option for (re, im)

pub fn write_hdf5(file_path: &str, simulation_config: &SimulationConfig, ray_paths: Vec<Vec<[f32; 3]>>, pressure_field: PressureField) -> Result<()> {

    // Create or overwrite HDF5 file
    let file = hdf5::File::create(file_path)?;

    // create meta group and write source info
    let src = file.create_group("src")?;
    src.new_dataset_builder().with_data(&simulation_config.source.freq_hz).create("frequency_hz")?;
    src.new_dataset_builder().with_data(&simulation_config.source.position).create("source_position_m")?;
    src.new_dataset_builder().with_data(&simulation_config.source.launch_elev_deg).create("launch_elev_deg")?;
    src.new_dataset_builder().with_data(&simulation_config.source.launch_azim_deg).create("launch_azim_deg")?;

    // ray path histories
    let rays = file.create_group("ray_paths")?;
    for (i, path) in ray_paths.iter().enumerate() {
        let ds_name = format!("ray_{}", i);
        rays.new_dataset_builder().with_data(path).create(ds_name.as_str())?;
    }

    let bty = file.create_group("bty")?;
    bty.new_dataset_builder().with_data(&simulation_config.bathymetry.x_bty_m).create("x_bty_m")?;
    bty.new_dataset_builder().with_data(&simulation_config.bathymetry.y_bty_m).create("y_bty_m")?;
    bty.new_dataset_builder().with_data(&simulation_config.bathymetry.z_bty_m).create("z_bty_m")?;

    let pressure = file.create_group("pressure_field")?;
    pressure.new_dataset_builder().with_data(&simulation_config.source.freq_hz).create("frequency_hz")?;
    if pressure_field.is_array {
        // write explicit receiver positions as an (N,3) dataset
        let recs = pressure_field.receiver_positions.as_ref().expect("receiver_positions present for array mode");
        let nrec = recs.len();
        let mut flat: Vec<f32> = Vec::with_capacity(nrec * 3);
        for r in recs.iter() { flat.push(r[0]); flat.push(r[1]); flat.push(r[2]); }
        let arr = ndarray::Array2::from_shape_vec((nrec, 3), flat)
            .expect("failed to shape receiver positions");
        pressure.new_dataset_builder().with_data(&arr).create("receiver_positions_m")?;
    } else {
        pressure.new_dataset_builder().with_data(&pressure_field.x_m).create("x_m")?;
        pressure.new_dataset_builder().with_data(&pressure_field.y_m).create("y_m")?;
            pressure.new_dataset_builder().with_data(&pressure_field.z_m).create("z_m")?;
    }
    // write per-receiver earliest arrival delay and amplitude
    pressure.new_dataset_builder().with_data(&pressure_field.delay_s).create("delay_s")?;
    pressure.new_dataset_builder().with_data(&pressure_field.amplitude).create("amplitude")?;
    
    // write complex pressure field as separate real and imaginary 3D datasets
    // pressure_field.pressure is an ndarray::Array4<num_complex::Complex32> with shape (nfreq, nx, ny, nz)
    let shape = pressure_field.pressure.dim(); // (nfreq, nx, ny, nz)
    let (nfreq, nx, ny, nz) = (shape.0, shape.1, shape.2, shape.3);

    // allocate flattened vectors for real and imag parts
    let mut re_flat: Vec<f32> = Vec::with_capacity(nfreq * nx * ny * nz);
    let mut im_flat: Vec<f32> = Vec::with_capacity(nfreq * nx * ny * nz);
    for ifreq in 0..nfreq {
        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    let v = pressure_field.pressure[(ifreq, ix, iy, iz)];
                    re_flat.push(v.re);
                    im_flat.push(v.im);
                }
            }
        }
    }

    // convert flattened data into 4D arrays with shape (nfreq, nx, ny, nz)
    let re_arr: Array4<f32> = Array4::from_shape_vec((nfreq, nx, ny, nz), re_flat)
        .expect("failed to reshape real pressure into 4D array");
    let im_arr: Array4<f32> = Array4::from_shape_vec((nfreq, nx, ny, nz), im_flat)
        .expect("failed to reshape imag pressure into 4D array");

    pressure.new_dataset_builder().with_data(&re_arr).create("pressure_re")?;
    pressure.new_dataset_builder().with_data(&im_arr).create("pressure_im")?;
    

    

    Ok(())
}
