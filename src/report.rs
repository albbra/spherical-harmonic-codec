//! Coefficient dump writers: human-readable `.log` and machine-readable `.csv`.

use std::fs::File;
use std::io::{BufWriter, Write};

use crate::codec::RoundTripStats;
use crate::error::Error;
use crate::vec3::Vec3;

#[allow(clippy::too_many_arguments)]
pub fn write_log(
    prefix: &str,
    input: &str,
    w: usize,
    h: usize,
    lmax: usize,
    hit_l: usize,
    rmse_tol: f64,
    a_lm: &[Vec3],
    aw: usize,
    stats: &RoundTripStats,
) -> Result<(), Error> {
    let path = format!("{}_L{:03}.log", prefix, hit_l);
    let fp = File::create(&path).map_err(Error::Io)?;
    let mut out = BufWriter::new(fp);

    writeln!(out, "# Spherical harmonic decomposition")?;
    writeln!(out, "# input      : {}", input)?;
    writeln!(out, "# image size : {} x {}", w, h)?;
    writeln!(out, "# channels   : RGB (3)")?;
    writeln!(out, "# maxL       : {}", lmax)?;
    writeln!(out, "# rmseTol    : {:.3e}", rmse_tol)?;
    writeln!(out, "# chosen L   : {}", hit_l)?;
    writeln!(
        out,
        "# coefficients: {} per channel, {} total",
        (hit_l + 1) * (hit_l + 1),
        3 * (hit_l + 1) * (hit_l + 1)
    )?;
    writeln!(out, "# Parseval RMSE  : {:.6}", stats.parseval_rmse)?;
    writeln!(out, "#")?;
    writeln!(out, "# Reconstruction quality:")?;
    writeln!(out, "#   RMSE (rgb avg)   : {:.6}", stats.rmse_rgb)?;
    writeln!(out, "#   RMSE (weighted)  : {:.6}", stats.rmse_weighted)?;
    writeln!(out, "#   MAE              : {:.6}", stats.mae_rgb)?;
    writeln!(out, "#   MAX error        : {:.6}", stats.max_err)?;
    writeln!(out, "#   PSNR             : {:.2} dB", stats.psnr_db)?;
    writeln!(out, "#   RMSE top 5%      : {:.6}", stats.band_rmse[0])?;
    writeln!(out, "#   RMSE middle 90%  : {:.6}", stats.band_rmse[1])?;
    writeln!(out, "#   RMSE bottom 5%   : {:.6}", stats.band_rmse[2])?;
    writeln!(out, "#")?;
    writeln!(out, "# Convention:")?;
    writeln!(out, "#   theta in [0,pi] (polar, 0=N pole), phi in [0,2pi)")?;
    writeln!(out, "#   Y_lm real, orthonormal, no Condon-Shortley phase.")?;
    writeln!(out, "#   For scipy.sph_harm, multiply m>0 by (-1)^m.")?;
    writeln!(out, "#")?;
    writeln!(out, "# Columns: l  m  a_lm_R  a_lm_G  a_lm_B")?;
    writeln!(out, "#")?;

    for l in 0..=hit_l {
        for m in -(l as i32)..=(l as i32) {
            let c = a_lm[l * aw + (m + lmax as i32) as usize];
            writeln!(
                out,
                "{:4} {:4}  {:.10e}  {:.10e}  {:.10e}",
                l, m, c.r, c.g, c.b
            )?;
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn write_csv(
    prefix: &str,
    input: &str,
    w: usize,
    h: usize,
    lmax: usize,
    hit_l: usize,
    a_lm: &[Vec3],
    aw: usize,
    stats: &RoundTripStats,
) -> Result<(), Error> {
    let path = format!("{}_L{:03}.csv", prefix, hit_l);
    let fp = File::create(&path).map_err(Error::Io)?;
    let mut out = BufWriter::new(fp);

    writeln!(
        out,
        "# input={} W={} H={} L={} rmse={:.6} psnr={:.2}",
        input, w, h, hit_l, stats.rmse_rgb, stats.psnr_db
    )?;
    writeln!(out, "l,m,a_r,a_g,a_b")?;

    for l in 0..=hit_l {
        for m in -(l as i32)..=(l as i32) {
            let c = a_lm[l * aw + (m + lmax as i32) as usize];
            writeln!(out, "{},{},{:.10e},{:.10e},{:.10e}", l, m, c.r, c.g, c.b)?;
        }
    }
    Ok(())
}
