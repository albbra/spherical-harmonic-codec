//! Round-trip quality metrics: RMSE, MAE, PSNR, band RMSE, and per-channel
//! max error, all computed on reconstructed vs. original pixels.

use crate::codec::RoundTripStats;
use crate::vec3::Vec3;

#[allow(clippy::too_many_arguments)]
pub fn round_trip(
    original: &[Vec3],
    recon: &[Vec3],
    w: usize,
    h: usize,
    sintheta: &[f64],
    dtheta: f64,
    dphi: f64,
    total_w: f64,
) -> RoundTripStats {
    let n = (w * h) as f64;
    let mut sse = [0.0f64; 3];
    let mut sae = [0.0f64; 3];
    let mut max_err = [0.0f64; 3];
    let mut sse_w = 0.0f64;

    for (y, &st) in sintheta.iter().enumerate().take(h) {
        let wgt = st * dtheta * dphi;
        for x in 0..w {
            let i = y * w + x;
            let o = original[i];
            let r = recon[i];
            let e = [r.r - o.r, r.g - o.g, r.b - o.b];
            for c in 0..3 {
                sse[c] += e[c] * e[c];
                sae[c] += e[c].abs();
                if e[c].abs() > max_err[c] {
                    max_err[c] = e[c].abs();
                }
            }
            sse_w += wgt * (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
        }
    }

    let mut total_sse = 0.0;
    let mut total_sae = 0.0;
    let mut total_max = 0.0;
    for c in 0..3 {
        total_sse += sse[c];
        total_sae += sae[c];
        if max_err[c] > total_max {
            total_max = max_err[c];
        }
    }
    let rmse_rgb = (total_sse / (3.0 * n)).sqrt();
    let mae_rgb = total_sae / (3.0 * n);
    let psnr_db = -20.0 * rmse_rgb.max(1e-12).log10();
    let rmse_weighted = (sse_w / (3.0 * total_w)).sqrt();

    let band = |y0: usize, y1: usize| -> f64 {
        let mut s = 0.0;
        let mut cnt = 0usize;
        for y in y0..y1 {
            for x in 0..w {
                let o = original[y * w + x];
                let r = recon[y * w + x];
                s += (r.r - o.r).powi(2) + (r.g - o.g).powi(2) + (r.b - o.b).powi(2);
                cnt += 3;
            }
        }
        (s / cnt.max(1) as f64).sqrt()
    };
    let top_n = (h / 20).max(1);
    let bot_n = (h / 20).max(1);
    let band_rmse = [band(0, top_n), band(top_n, h - bot_n), band(h - bot_n, h)];

    RoundTripStats {
        rmse_weighted,
        rmse_rgb,
        mae_rgb,
        psnr_db,
        band_rmse,
        max_err: total_max,
        parseval_rmse: 0.0,
    }
}
