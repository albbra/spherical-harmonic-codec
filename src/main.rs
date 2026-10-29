use clap::Parser;
use spherical_harmonic_codec::{Codec, Config};

#[derive(Parser)]
#[command(name = "shc", version, about = "Spherical Harmonic Codec")]
struct Cli {
    /// Input equirectangular PNG.
    input: std::path::PathBuf,
    /// Max spherical-harmonic band (default: min(W,H)/2).
    #[arg(long)]
    max_l: Option<usize>,
    /// Parseval-RMSE target for L selection (0 = keep all bands).
    #[arg(long, default_value_t = 1e-3)]
    rmse_tol: f64,
    /// Output prefix (default: <input>_sh).
    #[arg(long)]
    out_prefix: Option<String>,
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    let mut cfg = Config::from_path(&cli.input)?.with_rmse_tol(cli.rmse_tol);
    if let Some(l) = cli.max_l {
        cfg = cfg.with_max_l(l);
    }
    let prefix = cli
        .out_prefix
        .unwrap_or_else(|| format!("{}_sh", cfg.stem()));

    let codec = Codec::run(cfg)?;
    codec.print_summary();
    codec.write_images(&prefix)?;
    codec.write_reports(&prefix)?;
    Ok(())
}
