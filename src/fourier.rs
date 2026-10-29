//! Step 1: real Fourier transform over azimuth φ, producing per-`m`
//! coefficients `A[m](θ)` and `B[m](θ)`.

use crate::vec3::Vec3;
use rayon::prelude::*;

/// Step 1: Fourier transform over φ. Returns `(A, B)` each `(lmax+1) * h`,
/// indexed `m * h + y`.
#[allow(clippy::too_many_arguments)]
pub fn run(
    f: &[Vec3],
    w: usize,
    h: usize,
    lmax: usize,
    cos_t: &[f64],
    sin_t: &[f64],
    stride: usize,
    dphi: f64,
) -> (Vec<Vec3>, Vec<Vec3>) {
    let mut a_buf = vec![Vec3::ZERO; (lmax + 1) * h];
    let mut b_buf = vec![Vec3::ZERO; (lmax + 1) * h];

    a_buf
        .par_chunks_mut(h)
        .zip(b_buf.par_chunks_mut(h))
        .enumerate()
        .for_each(|(m, (arow, brow))| {
            for y in 0..h {
                let mut acc_a = Vec3::ZERO;
                let mut acc_b = Vec3::ZERO;
                for x in 0..w {
                    let fx = f[y * w + x];
                    let ct = cos_t[x * stride + m];
                    let st = sin_t[x * stride + m];
                    acc_a.r += fx.r * ct;
                    acc_a.g += fx.g * ct;
                    acc_a.b += fx.b * ct;
                    acc_b.r += fx.r * st;
                    acc_b.g += fx.g * st;
                    acc_b.b += fx.b * st;
                }
                arow[y] = acc_a * dphi;
                brow[y] = acc_b * dphi;
            }
        });

    (a_buf, b_buf)
}
