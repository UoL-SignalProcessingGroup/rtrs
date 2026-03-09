mod bty;
mod engine;
mod influence;
mod input;
mod output;
mod rays;
mod reflect;
mod ssp;
mod utils;

use anyhow::Result;
use std::fs;

use input::config::SimulationConfig;
use output::write_json;

fn load_config(path: &str) -> Result<SimulationConfig> {
    let text = fs::read_to_string(path)?;
    let mut unknown_fields = Vec::new();
    let mut de = serde_json::Deserializer::from_str(&text);
    let mut cfg: SimulationConfig = serde_ignored::deserialize(&mut de, |path| {
        unknown_fields.push(path.to_string());
    })?;
    let mut warnings = cfg.validate()?;
    warnings.extend(
        unknown_fields
            .into_iter()
            .map(|field| format!("unknown input key is ignored: config.{}", field)),
    );

    for w in warnings {
        eprintln!("WARNING {}", w);
    }
    Ok(cfg)
}

fn main() -> Result<()> {
    // Load config from JSON
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!(
            "Usage: {} <config.json>",
            args.get(0).unwrap_or(&"rtrs".into())
        );
        return Ok(());
    }
    let in_path = &args[1];
    let config = load_config(in_path)?;

    let (ray_paths, pressure_field) = engine::core(&config);

    // Derive output path: replace extension with .json (append if none)
    let out_path = {
        let mut p = std::path::PathBuf::from(in_path);
        p.set_extension("out.json");
        p
    };

    // Write JSON output
    let out_path_str = out_path.to_str().expect("Invalid output path");
    write_json(out_path_str, &config, ray_paths.as_ref(), &pressure_field)?;
    Ok(())
}
