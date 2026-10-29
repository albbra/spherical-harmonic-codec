//! Error type for all fallible operations in the crate.

use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("input image not found: {0}")]
    InputMissing(PathBuf),
    #[error("failed to decode '{path}': {source}")]
    Decode {
        path: PathBuf,
        #[source]
        source: image::ImageError,
    },
    #[error("failed to write '{path}': {source}")]
    Encode {
        path: PathBuf,
        #[source]
        source: image::ImageError,
    },
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("invalid configuration: {0}")]
    Config(String),
}
