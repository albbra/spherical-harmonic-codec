//! Codec configuration: input path, band cap, and RMSE tolerance.

use crate::error::Error;
use std::path::PathBuf;

/// Configuration for a single codec run.
///
/// Build with [`Config::from_path`] (which validates the input exists),
/// then chain the `with_*` builders. Pass to [`crate::Codec::run`].
///
/// # Defaults
///
/// - `max_l`    → `None` → `min(W, H) / 2` (Nyquist)
/// - `rmse_tol` → `1e-3`
///
/// # Example
///
/// ```no_run
/// # use spherical_harmonic_codec::Config;
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let cfg = Config::from_path("sky.png")?
///     .with_rmse_tol(1e-3)
///     .with_max_l(64);
/// # Ok(()) }
/// ```
pub struct Config {
    /// Input image path. Checked for existence by [`Config::from_path`].
    pub input: PathBuf,
    /// Cap on spherical-harmonic band `L`. `None` means `min(W, H) / 2`.
    pub max_l: Option<usize>,
    /// Parseval-RMSE target for band selection. `0` disables early
    /// termination (all bands up to `max_l` are kept).
    pub rmse_tol: f64,
}

impl Config {
    /// Create a config from a path, validating that the file exists.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InputMissing`] if `p` doesn't resolve to an existing
    /// file.
    pub fn from_path(p: impl Into<PathBuf>) -> Result<Self, Error> {
        let input = p.into();
        if !input.exists() {
            return Err(Error::InputMissing(input));
        }
        Ok(Self {
            input,
            max_l: None,
            rmse_tol: 1e-3,
        })
    }

    /// Override the maximum band `L`.
    ///
    /// The actual `L` is still chosen by Parseval's residual-energy
    /// criterion, capped at this value.
    pub fn with_max_l(mut self, l: usize) -> Self {
        self.max_l = Some(l);
        self
    }

    /// Override the Parseval-RMSE tolerance.
    ///
    /// `0` disables early termination — all bands up to `max_l` are kept.
    pub fn with_rmse_tol(mut self, t: f64) -> Self {
        self.rmse_tol = t;
        self
    }

    /// The input file's stem (no directory, no extension), for default
    /// output prefixes. Falls back to `"output"` if the stem isn't valid
    /// UTF-8.
    pub fn stem(&self) -> String {
        self.input
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output")
            .to_owned()
    }
}
