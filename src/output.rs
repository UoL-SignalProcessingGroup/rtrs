use anyhow::Result;
use hdf5_metno as hdf5;

use crate::input::config::SimulationConfig;


pub fn write_hdf5(file_path: &str, simulation_config: &SimulationConfig, ray_paths: Vec<Vec<[f64; 3]>>) -> Result<()> {

    // Create or overwrite HDF5 file
    let file = hdf5::File::create(file_path)?;

    // create meta group and write source info
    let src = file.create_group("src")?;
    src.new_dataset_builder().with_data(&[simulation_config.source.freq_hz]).create("frequency_hz")?;
    src.new_dataset_builder().with_data(&simulation_config.source.position).create("source_position_m")?;
    src.new_dataset_builder().with_data(&simulation_config.source.launch_elev_deg).create("launch_elev_deg")?;
    src.new_dataset_builder().with_data(&simulation_config.source.launch_azim_deg).create("launch_azim_deg")?;

    // ray path histories
    let rays = file.create_group("ray_paths")?;
    for (i, path) in ray_paths.iter().enumerate() {
        let ds_name = format!("ray_{}", i);
        rays.new_dataset_builder().with_data(path).create(ds_name.as_str())?;
    }


    Ok(())
}
