//! Usage: cargo run --example dump_coeffs -- sky.png
use spherical_harmonic_codec::{Codec, Config};

fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1).expect("usage: dump_coeffs <image>");
    let codec = Codec::run(Config::from_path(path)?)?;
    println!(
        "chosen L = {}, {} coefficients/channel",
        codec.hit_l,
        (codec.hit_l + 1) * (codec.hit_l + 1)
    );
    for l in 0..=codec.hit_l.min(3) {
        for m in -(l as i32)..=(l as i32) {
            let c = codec.a_lm[l * codec.aw + (m + codec.lmax as i32) as usize];
            println!(
                "  l={} m={:+}  R={:+.6}  G={:+.6}  B={:+.6}",
                l, m, c.r, c.g, c.b
            );
        }
    }
    Ok(())
}
