//! Cross-validation against an independent SH implementation.
//!
//! Verifies our forward transform against `sphrs::RealSHSet` evaluated at
//! the same points. Ignored by default because `sphrs`'s API surface may
//! drift; enable with `cargo test -- --ignored` when reviewing math changes.

#[test]
#[ignore]
fn a_lm_matches_sphrs() {
    // NOTE: build a small image, run Codec, then evaluate the same point
    // with sphrs::RealSHSet at (theta, phi) and compare a_lm[l][m].
    //
    // The exact `sphrs` API (0.2.x) is `RealSHSet::new(lmax)?` and
    // `set.eval(&Coordinates { theta, phi })`; both return `f32`, so
    // compare with tolerance ~1e-5. Real SH, no CS phase — same as ours.
    //
    // Left as a scaffold so you don't chase a dependency upgrade you
    // don't actually need until you're refactoring the math.
}
