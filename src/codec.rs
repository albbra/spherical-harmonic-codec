use crate::accumulate;
use crate::config::Config;
use crate::error::Error;
use crate::fourier;
use crate::metrics;
use crate::report;
use crate::synthesis;
use crate::vec3::Vec3;

/// Result of a complete codec run.
///
/// Populated by [`Codec::run`], then consumed by
/// [`print_summary`](Codec::print_summary),
/// [`write_images`](Codec::write_images), and
/// [`write_reports`](Codec::write_reports).
pub struct Codec {
    /// Image width in pixels.
    pub width: usize,
    /// Image height in pixels.
    pub height: usize,
    /// Nyquist cap on band `L` (from config or `min(W, H) / 2`).
    pub lmax: usize,
    /// Band actually used — smallest `L` whose Parseval residual met
    /// `rmse_tol`, or `lmax` if the tolerance couldn't be met.
    pub hit_l: usize,
    /// Tolerance actually used (copied from [`Config`]).
    pub rmse_tol: f64,
    /// Flat coefficient array, indexed `l * aw + (m + lmax)`.
    /// Length `(lmax + 1) * aw`, where `aw = 2 * lmax + 1`.
    pub a_lm: Vec<Vec3>,
    /// `2 * lmax + 1` — width of the `m` axis in [`a_lm`](Self::a_lm).
    pub aw: usize,
    /// Original pixel data (linearised RGB in `[0, 1]`, row-major).
    pub original: Vec<Vec3>,
    /// Reconstructed pixel data at band `hit_l`.
    pub recon: Vec<Vec3>,
    /// Round-trip quality metrics.
    pub stats: RoundTripStats,
    /// Input path as a display string (used in written reports).
    pub input: String,
}

/// Round-trip quality metrics computed after synthesis.
///
/// All values are in `[0, 1]` pixel range unless noted.
pub struct RoundTripStats {
    /// Solid-angle-weighted RMSE — what Parseval was minimising.
    pub rmse_weighted: f64,
    /// Unweighted RGB RMSE (every pixel equal weight).
    pub rmse_rgb: f64,
    /// Unweighted RGB mean absolute error.
    pub mae_rgb: f64,
    /// Peak signal-to-noise ratio in dB. `> 40 dB` is visually lossless.
    pub psnr_db: f64,
    /// RMSE over three horizontal bands: top 5% of rows, middle 90%,
    /// bottom 5%. Higher top/bottom than middle indicates polar ringing.
    pub band_rmse: [f64; 3],
    /// Largest per-channel absolute error over any single pixel.
    pub max_err: f64,
    /// Parseval residual-RMSE estimate at the chosen `L`.
    ///
    /// May be `0.0` when the estimate clamps (numerical overshoot on
    /// discrete grids). See [`Codec::run`] for details.
    pub parseval_rmse: f64,
}

impl Codec {
    /// Run the full forward → select-L → inverse pipeline.
    ///
    /// Steps:
    ///
    /// 1. Load the input as linearised `[0, 1]` RGB.
    /// 2. Forward SH transform ([`fourier::run`] then [`accumulate::run`]).
    /// 3. Parseval-based band selection — smallest `L` whose residual
    ///    RMSE meets [`Config::rmse_tol`].
    /// 4. Inverse transform ([`synthesis::run`]) at the chosen `L`.
    /// 5. Round-trip quality metrics ([`metrics::round_trip`]).
    ///
    /// # Errors
    ///
    /// Returns [`Error::Decode`] if the input can't be read as an image,
    /// or [`Error::Io`] for any intermediate I/O failure.
    ///
    /// # Numerical note — Parseval clamping
    ///
    /// The Parseval residual `‖f‖² - Σ|a_lm|²` can go slightly negative
    /// on a discrete grid because the discrete SH basis is not exactly
    /// orthonormal. The estimator clamps to `0.0`, so `parseval_rmse ==
    /// 0.0` should be read as "below numerical resolution", not
    /// "exact". This does not affect L selection (which only ever tests
    /// `rmse < tol` for `tol > 0`), but tests must not assume
    /// proportionality near zero.
    pub fn run(cfg: Config) -> Result<Self, Error> {
        let (f, w, h) = crate::load_rgb(&cfg.input)?;
        let lmax = cfg.max_l.unwrap_or_else(|| (w.min(h) / 2).max(1));
        // NB: lmax = min(W,H)/2 is Nyquist. That's the right default for
        // arbitrary images, but full-Lmax reconstruction at Nyquist aliases
        // (the discrete SH basis is non-orthogonal there). Tests that need
        // near-perfect reconstruction should cap max_l well below this.

        let dphi = 2.0 * std::f64::consts::PI / w as f64;
        let dtheta = std::f64::consts::PI / h as f64;

        let costheta: Vec<f64> = (0..h).map(|y| ((y as f64 + 0.5) * dtheta).cos()).collect();
        let sintheta: Vec<f64> = (0..h).map(|y| ((y as f64 + 0.5) * dtheta).sin()).collect();

        let mut total_w = 0.0f64;
        let mut f2 = 0.0f64;
        for y in 0..h {
            let wgt = sintheta[y] * dtheta * dphi;
            total_w += wgt;
            for x in 0..w {
                f2 += wgt * f[y * w + x].norm_sq();
            }
        }

        let stride = lmax + 1;
        let mut cos_t = vec![0.0f64; w * stride];
        let mut sin_t = vec![0.0f64; w * stride];
        for x in 0..w {
            let phi = (x as f64 + 0.5) * dphi;
            for m in 0..=lmax {
                cos_t[x * stride + m] = (m as f64 * phi).cos();
                sin_t[x * stride + m] = (m as f64 * phi).sin();
            }
        }

        println!("step 1/3: Fourier over phi ...");
        let (a_buf, b_buf) = fourier::run(&f, w, h, lmax, &cos_t, &sin_t, stride, dphi);

        println!("step 2/3: Legendre accumulation ...");
        let (a_lm, aw) = accumulate::run(&a_buf, &b_buf, h, lmax, &costheta, &sintheta, dtheta);
        drop(a_buf);
        drop(b_buf);

        // Parseval-based L selection
        println!("energy analysis ...");
        let mut running = f2;
        let mut hit_l = lmax;
        let mut parseval_rmse = 0.0;
        for l in 0..=lmax {
            for m in -(l as i32)..=(l as i32) {
                let idx = l * aw + (m + lmax as i32) as usize;
                running -= a_lm[idx].norm_sq();
            }
            let rmse = (running.max(0.0) / (3.0 * total_w)).sqrt();
            if l < 20 || l % 25 == 0 || rmse < cfg.rmse_tol {
                println!("  L={:3}  est.RMSE={:.6}", l, rmse);
            }
            if rmse < cfg.rmse_tol {
                hit_l = l;
                parseval_rmse = rmse;
                println!("chosen L = {}", l);
                break;
            }
            parseval_rmse = rmse;
        }
        if parseval_rmse > cfg.rmse_tol {
            println!("  tolerance not met; using L={}", lmax);
        }

        println!("step 3/3: synthesis at L={} ...", hit_l);
        let recon = synthesis::run(
            &a_lm, w, h, aw, &costheta, &sintheta, &cos_t, &sin_t, stride, hit_l,
        );

        let mut stats = metrics::round_trip(&f, &recon, w, h, &sintheta, dtheta, dphi, total_w);
        stats.parseval_rmse = parseval_rmse;

        Ok(Self {
            width: w,
            height: h,
            lmax,
            hit_l,
            rmse_tol: cfg.rmse_tol,
            a_lm,
            aw,
            original: f,
            recon,
            stats,
            input: cfg.input.display().to_string(),
        })
    }

    /// Print a human-readable reconstruction summary to stdout.
    pub fn print_summary(&self) {
        let s = &self.stats;
        println!("\n=== Reconstruction quality (L={}) ===", self.hit_l);
        println!(
            "  overall RMSE       {:.6}  ({:.4}% of range)",
            s.rmse_rgb,
            s.rmse_rgb * 100.0
        );
        println!("  overall MAE        {:.6}", s.mae_rgb);
        println!("  overall MAX err    {:.6}", s.max_err);
        println!("  weighted RMSE      {:.6}", s.rmse_weighted);
        println!("  PSNR               {:.2} dB", s.psnr_db);
        println!("  Parseval RMSE      {:.6}", s.parseval_rmse);
        println!("  RMSE top 5%%        {:.6}", s.band_rmse[0]);
        println!("  RMSE middle 90%%     {:.6}", s.band_rmse[1]);
        println!("  RMSE bottom 5%%     {:.6}", s.band_rmse[2]);
        println!("=====================================\n");
    }

    /// Write the reconstructed PNG and the 4× error heatmap.
    ///
    /// Files:
    /// - `<prefix>_L<###>.png` — reconstruction from SH coefficients only
    /// - `<prefix>_L<###>_diff4x.png` — 4× amplified `|error|` heatmap
    ///
    /// # Errors
    ///
    /// Returns [`Error::Encode`] if either PNG can't be written.
    pub fn write_images(&self, prefix: &str) -> Result<(), Error> {
        let clamp = |v: f64| (v.clamp(0.0, 1.0) * 255.0).round() as u8;

        let mut buf = Vec::with_capacity(self.width * self.height * 3);
        for p in &self.recon {
            buf.push(clamp(p.r));
            buf.push(clamp(p.g));
            buf.push(clamp(p.b));
        }
        let img = image::RgbImage::from_raw(self.width as u32, self.height as u32, buf)
            .ok_or_else(|| Error::Config("failed to build output image".into()))?;
        let path = format!("{}_L{:03}.png", prefix, self.hit_l);
        img.save(&path).map_err(|e| Error::Encode {
            path: path.clone().into(),
            source: e,
        })?;
        println!("wrote {}", path);

        let mut diff = Vec::with_capacity(self.width * self.height * 3);
        for i in 0..self.original.len() {
            let o = self.original[i];
            let r = self.recon[i];
            diff.push((((r.r - o.r).abs() * 4.0 * 255.0).min(255.0)) as u8);
            diff.push((((r.g - o.g).abs() * 4.0 * 255.0).min(255.0)) as u8);
            diff.push((((r.b - o.b).abs() * 4.0 * 255.0).min(255.0)) as u8);
        }
        let dimg = image::RgbImage::from_raw(self.width as u32, self.height as u32, diff)
            .ok_or_else(|| Error::Config("failed to build diff image".into()))?;
        let dpath = format!("{}_L{:03}_diff4x.png", prefix, self.hit_l);
        dimg.save(&dpath).map_err(|e| Error::Encode {
            path: dpath.clone().into(),
            source: e,
        })?;
        println!("wrote {}", dpath);

        Ok(())
    }

    /// Write the coefficient `.log` and `.csv` files.
    ///
    /// Files:
    /// - `<prefix>_L<###>.log` — human-readable coefficient dump + metrics
    /// - `<prefix>_L<###>.csv` — machine-readable `l, m, a_r, a_g, a_b`
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] if either file can't be created or written.
    pub fn write_reports(&self, prefix: &str) -> Result<(), Error> {
        report::write_log(
            prefix,
            &self.input,
            self.width,
            self.height,
            self.lmax,
            self.hit_l,
            self.rmse_tol,
            &self.a_lm,
            self.aw,
            &self.stats,
        )?;
        report::write_csv(
            prefix,
            &self.input,
            self.width,
            self.height,
            self.lmax,
            self.hit_l,
            &self.a_lm,
            self.aw,
            &self.stats,
        )?;
        println!("wrote {}_L{:03}.log", prefix, self.hit_l);
        println!("wrote {}_L{:03}.csv", prefix, self.hit_l);
        Ok(())
    }
}
