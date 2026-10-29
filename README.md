# spherical-harmonic-codec

Fast, parallel real spherical-harmonic (SH) analysis and synthesis for
equirectangular RGB images, with Parseval-based rate control and full
round-trip quality metrics.

Rust port of an earlier C++ tool (`sh_decompose`) — same math, same
output filenames, same `.csv`/`.log` formats.

## What it does

The crate runs the full analysis → select-L → synthesis → evaluate loop:

1. **Fourier over φ** — per row θ, take a real Fourier transform over the
   azimuth φ. Parallel over `m`.
2. **Legendre accumulation** — project onto fully-normalised associated
   Legendre functions to produce SH coefficients `a_lm`. Parallel with
   per-thread accumulators, merged via a reduce tree.
3. **Select L** — smallest band whose Parseval residual energy is below
   the caller's RMSE tolerance. No synthesis required for this step.
4. **Synthesis** — inverse transform back to pixels at the chosen band.
   Parallel over rows.
5. **Evaluate** — RMSE, MAE, PSNR, band RMSE (top 5% / middle 90% /
   bottom 5% of rows), and per-channel max error.

Unlike a pure "decomposer", the tool reconstructs the image and emits
quality metrics — that's what makes it a **codec** (coder + decoder).

## Install

From source:

    git clone https://github.com/albbra/spherical-harmonic-codec
    cd spherical-harmonic-codec
    cargo build --release

The binary is `shc`. The repository's `.cargo/config.toml` sets
`-C target-cpu=native` for the release profile.

## CLI usage

    shc <image.png> [options]

    Options:
      --max-l <N>         Max SH band (default: min(W,H)/2, the Nyquist limit)
      --rmse-tol <F>      Parseval-RMSE target for L selection (default: 1e-3;
                          0 keeps all bands up to --max-l)
      --out-prefix <STR>  Output prefix (default: <input-stem>_sh)
      -h, --help          Print help
      -V, --version       Print version

### Examples

    # Default: choose L automatically to hit ~1e-3 weighted RMSE
    shc sky.png

    # Cap at L=64 (useful for very high-res inputs where Nyquist is huge)
    shc sky.png --max-l 64

    # Keep all coefficients up to Nyquist (largest file, best reconstruction)
    shc sky.png --rmse-tol 0

    # Custom output prefix
    shc sky.png --out-prefix builds/sky_v2

## Outputs

For a chosen band `L`, with prefix `P`:

| File | Contents |
|---|---|
| `P_L###.png`           | Reconstruction from SH coefficients only |
| `P_L###_diff4x.png`    | 4× amplified \|error\| heatmap |
| `P_L###.log`           | Human-readable coefficient dump + metrics |
| `P_L###.csv`           | Machine-readable `l, m, a_r, a_g, a_b` |

The `.csv` is stable and machine-friendly; the `.log` is for reading.
Both include the SH convention in a header comment.

## Library usage

    use spherical_harmonic_codec::{Codec, Config};

    fn main() -> Result<(), Box<dyn std::error::Error>> {
        let cfg = Config::from_path("sky.png")?
            .with_rmse_tol(1e-3)
            .with_max_l(64);

        let codec = Codec::run(cfg)?;

        println!(
            "chosen L = {} ({} coefficients/channel), PSNR = {:.2} dB",
            codec.hit_l,
            (codec.hit_l + 1) * (codec.hit_l + 1),
            codec.stats.psnr_db,
        );

        codec.write_images("sky_sh")?;   // .png + _diff4x.png
        codec.write_reports("sky_sh")?;  // .log + .csv

        Ok(())
    }

The `Codec` struct exposes all intermediate state:

- `a_lm` — flat coefficient array, indexed `l * aw + (m + lmax)`
- `aw` — `2 * lmax + 1`, the width of the `m` axis
- `original`, `recon` — input and reconstructed pixels in `[0, 1]`
- `stats` — the full `RoundTripStats`

## Conventions

- `theta ∈ [0, π]` — polar angle, `0 = north pole`
- `phi ∈ [0, 2π)` — azimuth
- Real spherical harmonics, orthonormal on the unit sphere
- **No Condon–Shortley phase**

For `scipy.sph_harm` compatibility, multiply `m > 0` coefficients by
`(-1)^m`. The `.log` file repeats this note in its header.

## Numerical notes

**Nyquist aliasing.** The default band cap `lmax = min(W, H) / 2` is the
Nyquist limit. On a discrete equirectangular grid, reconstructing at
Nyquist aliases — the discrete SH basis is not exactly orthonormal there.
For arbitrary images this doesn't matter (tolerance-based L selection
lands well below Nyquist). For tests that want near-perfect
reconstruction, cap `max_l` well below `min(W, H) / 2`.

**Parseval clamps to zero.** The Parseval residual
`‖f‖² − Σ|a_lm|²` can go slightly negative from floating-point
cancellation. The estimator clamps to `0.0`. A reported
`parseval_rmse == 0.0` means "below numerical resolution", not "exact".

**8-bit quantization floor.** An 8-bit input has an irreducible RMSE
around 1.1e-3 purely from quantization. Functional tests use `1e-2` as
their threshold; `1e-3` and below would need float-valued input (not
currently supported by the loader).

## Project layout

    spherical-harmonic-codec/
    ├── Cargo.toml
    ├── .cargo/config.toml         # target-cpu=native for release
    ├── src/
    │   ├── lib.rs                 # crate root + load_rgb
    │   ├── main.rs                # thin CLI shell (clap)
    │   ├── config.rs              # Config builder
    │   ├── error.rs               # Error (thiserror)
    │   ├── vec3.rs                # Vec3 pixel/coefficient sample
    │   ├── legendre.rs            # normalised assoc. Legendre
    │   ├── fourier.rs             # step 1: Fourier over φ
    │   ├── accumulate.rs          # step 2: Legendre → a_lm
    │   ├── synthesis.rs           # step 3: inverse transform
    │   ├── metrics.rs             # RMSE / MAE / PSNR / band RMSE
    │   ├── report.rs              # .log + .csv writers
    │   └── codec.rs               # orchestration + RoundTripStats
    ├── tests/
    │   ├── roundtrip.rs           # 4 integration tests, no fixtures
    │   └── sphrs_agreement.rs     # #[ignore]d cross-check scaffold
    ├── examples/
    │   └── dump_coeffs.rs         # programmatic-use demo
    └── benches/
        └── transform.rs           # placeholder for criterion

The algorithm lives in the **library**; the binary is a 25-line shell
around `clap` → `Config` → `Codec::run(cfg)`. That split is what makes
`tests/`, `examples/`, and `benches/` possible without duplicating I/O
code.

## Performance

- Parallelised with `rayon`; per-thread accumulators merged via a reduce
  tree — no locks in the hot path
- Trig tables for `cos(m·phi)`, `sin(m·phi)` precomputed once
- Coarse-grained parallelism over rows (synthesis) or over `m` (forward
  Fourier) — minimal synchronization
- Compiled with `-C target-cpu=native` via `.cargo/config.toml`

## Tests

    cargo test --release

Runs four integration tests (synthetic images, no fixtures needed):

- `roundtrip_tight_tolerance` — tolerance-based L selection
- `roundtrip_full_lmax` — `--rmse-tol 0` keeps all bands
- `parseval_and_true_rmse_both_small_for_bandlimited` — both small for
  an exactly-representable signal
- `determinism` — repeated runs produce identical summary statistics

One additional test is `#[ignore]`d:
`tests/sphrs_agreement.rs::a_lm_matches_sphrs` — a scaffold for
cross-validation against the `sphrs` crate. Enable with
`cargo test -- --ignored` after filling in the body.

## Contributing

Issues and PRs welcome. Before submitting:

    cargo fmt --check
    cargo clippy --all-targets -- -D warnings
    cargo test --release