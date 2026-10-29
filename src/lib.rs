//! # Spherical Harmonic Codec
//!
//! Forward and inverse real spherical-harmonic transform for equirectangular
//! RGB images, with Parseval-based rate control and full round-trip quality
//! metrics (RMSE, MAE, PSNR, band RMSE).
//!
//! # Pipeline
//!
//! The codec runs in three stages, plus a rate-selection step in between:
//!
//! 1. **Fourier over φ** ([`fourier`]) — per row θ, take a real Fourier
//!    transform over the azimuth φ. Parallel over `m`.
//! 2. **Legendre accumulation** ([`accumulate`]) — combine per-`m` Fourier
//!    coefficients with fully-normalised associated Legendre functions to
//!    produce SH coefficients `a_lm`. Parallel with per-thread accumulators,
//!    merged via a reduce tree.
//! 3. **Parseval selection** — pick the smallest band `L` whose residual
//!    energy falls below the caller's RMSE tolerance. No synthesis required.
//! 4. **Synthesis** ([`synthesis`]) — inverse transform back to pixels at
//!    the chosen band `L`. Parallel over rows.
//!
//! # Example
//!
//! ```no_run
//! use spherical_harmonic_codec::{Codec, Config};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let cfg = Config::from_path("sky.png")?
//!     .with_rmse_tol(1e-3)
//!     .with_max_l(64);
//!
//! let codec = Codec::run(cfg)?;
//!
//! println!(
//!     "chosen L = {} ({} coefficients/channel), PSNR = {:.2} dB",
//!     codec.hit_l,
//!     (codec.hit_l + 1) * (codec.hit_l + 1),
//!     codec.stats.psnr_db,
//! );
//!
//! codec.write_images("sky_sh")?;   // reconstructed + diff heatmap PNGs
//! codec.write_reports("sky_sh")?;  // .log + .csv coefficient dumps
//! # Ok(())
//! # }
//! ```
//!
//! # Conventions
//!
//! - `theta ∈ [0, π]` — polar angle, `0 = north pole`
//! - `phi ∈ [0, 2π)` — azimuth
//! - Real spherical harmonics, orthonormal on the unit sphere, **no**
//!   Condon–Shortley phase. For `scipy.sph_harm` compatibility, multiply
//!   `m > 0` coefficients by `(-1)^m`.
//!
//! # Numerical notes
//!
//! - The default band cap `lmax = min(W, H) / 2` is the Nyquist limit. On a
//!   discrete equirectangular grid, reconstructing *at* Nyquist aliases —
//!   the discrete SH basis is not exactly orthonormal there. Tests that
//!   want near-perfect reconstruction should cap `max_l` well below this.
//! - Parseval's residual `‖f‖² - Σ|a_lm|²` can go slightly negative from
//!   floating-point cancellation. [`Codec::run`] clamps to `0.0`, so
//!   `parseval_rmse == 0.0` means "below numerical resolution", not
//!   "exact". See [`Codec::run`] for details.

pub mod accumulate;
pub mod config;
pub mod error;
pub mod fourier;
pub mod legendre;
pub mod metrics;
pub mod report;
pub mod synthesis;
pub mod vec3;

mod codec;

pub use codec::{Codec, RoundTripStats};
pub use config::Config;
pub use error::Error;

use image::RgbImage;
use std::path::Path;

/// Load an image as RGB and return `(pixels, width, height)`.
///
/// Pixels are normalised to `[0.0, 1.0]` and stored row-major.
pub(crate) fn load_rgb(path: &Path) -> Result<(Vec<vec3::Vec3>, usize, usize), Error> {
    let img = image::open(path).map_err(|e| Error::Decode {
        path: path.to_owned(),
        source: e,
    })?;
    let rgb: RgbImage = img.to_rgb8();
    let (w, h) = (rgb.width() as usize, rgb.height() as usize);
    let raw = rgb.as_raw();
    let mut out = Vec::with_capacity(w * h);
    for i in (0..raw.len()).step_by(3) {
        out.push(vec3::Vec3 {
            r: raw[i] as f64 / 255.0,
            g: raw[i + 1] as f64 / 255.0,
            b: raw[i + 2] as f64 / 255.0,
        });
    }
    Ok((out, w, h))
}
