use image::{Rgb, RgbImage};
use std::sync::atomic::{AtomicUsize, Ordering};

use spherical_harmonic_codec::{Codec, Config};

static COUNTER: AtomicUsize = AtomicUsize::new(0);

fn make_input(w: u32, h: u32) -> (std::path::PathBuf, TempGuard) {
    let id = COUNTER.fetch_add(1, Ordering::Relaxed);
    let path = std::env::temp_dir().join(format!(
        "shc_test_{}_{}_{}_{}.png",
        std::process::id(),
        id,
        w,
        h
    ));
    let mut img = RgbImage::new(w, h);
    for y in 0..h {
        for x in 0..w {
            let theta = (y as f32 + 0.5) / h as f32 * std::f32::consts::PI;
            let phi = (x as f32 + 0.5) / w as f32 * 2.0 * std::f32::consts::PI;
            let r = 0.5 + 0.25 * theta.sin() * phi.cos();
            let g = 0.5 + 0.25 * theta.cos();
            let b = 0.5;
            img.put_pixel(
                x,
                y,
                Rgb([(r * 255.0) as u8, (g * 255.0) as u8, (b * 255.0) as u8]),
            );
        }
    }
    img.save(&path).expect("failed to write test image");
    let guard = TempGuard(path.clone());
    (path, guard)
}

struct TempGuard(std::path::PathBuf);
impl Drop for TempGuard {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

#[test]
fn roundtrip_tight_tolerance() {
    let (path, _guard) = make_input(64, 32);
    let cfg = Config::from_path(&path).unwrap().with_rmse_tol(1e-2);
    let codec = Codec::run(cfg).unwrap();

    assert!(
        codec.stats.rmse_weighted < 1e-2,
        "weighted RMSE {} exceeded 1e-2",
        codec.stats.rmse_weighted
    );
    assert!(
        codec.stats.psnr_db > 30.0,
        "PSNR too low: {:.2} dB",
        codec.stats.psnr_db
    );
    assert!(
        codec.stats.parseval_rmse < 1e-2,
        "Parseval RMSE {} did not meet tolerance",
        codec.stats.parseval_rmse
    );
    assert!(
        codec.hit_l <= 3,
        "Parseval should have chosen L<=3 for a low-order SH test image, got {}",
        codec.hit_l
    );
}

#[test]
fn roundtrip_full_lmax() {
    // 64×32 grid: Nyquist at L=16. Test image is pure L=1, so a
    // max_l of 1 is safely oversampled — no aliasing in the inverse.
    // (A 32×16 grid at lmax=8 is *at* Nyquist and will alias; do not
    // confuse that for a codec bug.)
    let (path, _guard) = make_input(64, 32);
    let cfg = Config::from_path(&path)
        .unwrap()
        .with_max_l(1)
        .with_rmse_tol(0.0);
    let codec = Codec::run(cfg).unwrap();

    assert_eq!(codec.hit_l, 1, "rmse_tol=0 should keep all L");
    assert!(
        codec.stats.rmse_weighted < 1e-2,
        "full-L reconstruction of a bandlimited signal should be clean, got {}",
        codec.stats.rmse_weighted
    );
}

#[test]
fn parseval_and_true_rmse_both_small_for_bandlimited() {
    // For a signal that IS an SH mode, both Parseval (energy-domain)
    // and the true spatial RMSE should collapse to near quantization
    // floor. Parseval may clamp to 0 from numerical overshoot —
    // that's expected; we assert "small", not "proportional".
    let (path, _guard) = make_input(64, 32);
    let cfg = Config::from_path(&path)
        .unwrap()
        .with_max_l(1)
        .with_rmse_tol(1e-6);
    let codec = Codec::run(cfg).unwrap();

    assert!(
        codec.hit_l <= 2,
        "should pick a small L for a bandlimited image, got {}",
        codec.hit_l
    );
    assert!(
        codec.stats.parseval_rmse < 1e-3,
        "Parseval should collapse for bandlimited signal, got {:.3e}",
        codec.stats.parseval_rmse
    );
    assert!(
        codec.stats.rmse_weighted < 1e-2,
        "true RMSE should be small for bandlimited signal, got {:.3e}",
        codec.stats.rmse_weighted
    );
}

#[test]
fn determinism() {
    let (path, _guard) = make_input(32, 16);
    let cfg = || Config::from_path(&path).unwrap().with_rmse_tol(1e-3);
    let a = Codec::run(cfg()).unwrap();
    let b = Codec::run(cfg()).unwrap();
    assert_eq!(a.hit_l, b.hit_l);
    assert!((a.stats.rmse_weighted - b.stats.rmse_weighted).abs() < 1e-12);
}
