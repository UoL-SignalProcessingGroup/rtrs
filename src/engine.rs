use crate::bty;
use crate::influence::{PressureField, gaussian_beam_influence, init_pressure_field};
use crate::input::config::SimulationConfig;
use crate::rays::trace_ray;
use crate::ssp::init_ssp;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::Zip;
use rayon::prelude::*;
use std::f32::consts::PI;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::Instant;

fn merge_pressure_fields(dst: &mut PressureField, src: &PressureField) {
    Zip::from(&mut dst.pressure)
        .and(&src.pressure)
        .for_each(|d, s| *d += *s);

    Zip::from(&mut dst.delay_s)
        .and(&mut dst.amplitude)
        .and(&src.delay_s)
        .and(&src.amplitude)
        .for_each(|d_delay, d_amp, &s_delay, &s_amp| {
            if s_delay < *d_delay {
                *d_delay = s_delay;
                *d_amp = s_amp;
            }
        });
}

pub fn core(cfg: &SimulationConfig) -> (Option<Vec<Vec<[f32; 3]>>>, PressureField) {
    // convert angles to radians
    let launch_elev_rad: Vec<f32> = cfg
        .source
        .launch_elev_deg
        .iter()
        .map(|d| d.to_radians())
        .collect(); // "alpha" in Fortran
    let launch_azim_rad: Vec<f32> = cfg
        .source
        .launch_azim_deg
        .iter()
        .map(|d| d.to_radians())
        .collect(); // "beta" in Fortran
    let d_elev = if launch_elev_rad.len() >= 2 {
        (launch_elev_rad[1] - launch_elev_rad[0]).abs()
    } else {
        1.0_f32.to_radians()
    };
    let d_azim = if launch_azim_rad.len() >= 2 {
        (launch_azim_rad[1] - launch_azim_rad[0]).abs()
    } else {
        1.0_f32.to_radians()
    };

    // initialize environmental fields
    let ssp_field = init_ssp(cfg);
    let bty_field = bty::init_bty(cfg);

    // angular frequency
    let omega: Vec<f32> = cfg.source.freq_hz.iter().map(|f| 2.0 * PI * f).collect();

    // Create all angle pairs upfront for parallel processing
    let angle_pairs: Vec<(f32, f32)> = launch_azim_rad
        .iter()
        .flat_map(|&azim| launch_elev_rad.iter().map(move |&elev| (azim, elev)))
        .collect();

    // Process rays in parallel with worker-local accumulation and reduce
    let store_ray_paths = cfg.beam.store_ray_paths;
    let show_progress = cfg.beam.show_progress;
    let use_atomic_progress_counter = cfg.beam.atomic_progress_counter;
    let n_rays = angle_pairs.len();
    let target_chunks = (rayon::current_num_threads() * 2).max(1);
    let chunk_size = if n_rays == 0 {
        1
    } else {
        (n_rays + target_chunks - 1) / target_chunks
    };
    let progress = if show_progress && !use_atomic_progress_counter && n_rays > 0 {
        let pb = ProgressBar::new(n_rays as u64);
        let style = ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {percent:>3}% eta {eta_precise}",
        )
        .unwrap_or_else(|_| ProgressStyle::default_bar())
        .progress_chars("=>-");
        pb.set_style(style);
        Some(pb)
    } else {
        None
    };
    let progress_for_workers = progress.clone();
    let completed = if show_progress && use_atomic_progress_counter && n_rays > 0 {
        Some(Arc::new(AtomicUsize::new(0)))
    } else {
        None
    };
    let last_print_ms = if show_progress && use_atomic_progress_counter && n_rays > 0 {
        Some(Arc::new(AtomicU64::new(0)))
    } else {
        None
    };
    let start = if show_progress && use_atomic_progress_counter && n_rays > 0 {
        Some(Instant::now())
    } else {
        None
    };

    let (mut indexed_paths, pressure_field) = angle_pairs
        .par_chunks(chunk_size)
        .enumerate()
        .map(|(chunk_idx, chunk)| {
            let mut chunk_paths = if store_ray_paths {
                Vec::<(usize, Vec<[f32; 3]>)>::with_capacity(chunk.len())
            } else {
                Vec::new()
            };
            let mut chunk_pressure = init_pressure_field(cfg);
            let base_idx = chunk_idx * chunk_size;

            for (offset, &(azim, elev)) in chunk.iter().enumerate() {
                let idx = base_idx + offset;
                let mut ray_history = trace_ray(azim, elev, cfg, &ssp_field, &bty_field);

                gaussian_beam_influence(
                    &mut ray_history,
                    &mut chunk_pressure,
                    &bty_field,
                    elev,
                    d_azim,
                    d_elev,
                    &omega,
                );

                if store_ray_paths {
                    let path = ray_history
                        .iter()
                        .map(|r| r.position)
                        .collect::<Vec<[f32; 3]>>();
                    chunk_paths.push((idx, path));
                }
                if let Some(pb) = &progress_for_workers {
                    pb.inc(1);
                }
            }

            if let (Some(completed), Some(last_print_ms), Some(start)) =
                (&completed, &last_print_ms, &start)
            {
                let done = completed.fetch_add(chunk.len(), Ordering::Relaxed) + chunk.len();
                let elapsed = start.elapsed();
                let elapsed_ms = elapsed.as_millis() as u64;
                let prev_ms = last_print_ms.load(Ordering::Relaxed);

                // Limit prints to roughly every 5 seconds, always print at completion.
                if elapsed_ms.saturating_sub(prev_ms) >= 5000 || done == n_rays {
                    if last_print_ms
                        .compare_exchange(prev_ms, elapsed_ms, Ordering::SeqCst, Ordering::Relaxed)
                        .is_ok()
                    {
                        let pct = 100.0 * done as f64 / n_rays as f64;
                        let secs = elapsed.as_secs_f64();
                        let rate = done as f64 / secs.max(1e-9);
                        let remaining = n_rays.saturating_sub(done);
                        let eta = remaining as f64 / rate.max(1e-9);

                        println!(
                            "Progress: {done}/{n_rays} ({pct:.1}%) | elapsed {secs:.1}s | {rate:.0} rays/s | ETA {eta:.1}s"
                        );
                    }
                }
            }

            (chunk_paths, chunk_pressure)
        })
        .reduce(
            || {
                (
                    Vec::<(usize, Vec<[f32; 3]>)>::new(),
                    init_pressure_field(cfg),
                )
            },
            |mut left, right| {
                left.0.extend(right.0);
                merge_pressure_fields(&mut left.1, &right.1);
                left
            },
        );
    if let Some(pb) = progress {
        pb.finish();
    }

    let ray_paths = if store_ray_paths {
        indexed_paths.sort_by_key(|(idx, _)| *idx);
        Some(
            indexed_paths
                .into_iter()
                .map(|(_, path)| path)
                .collect::<Vec<_>>(),
        )
    } else {
        None
    };

    return (ray_paths, pressure_field);
}
