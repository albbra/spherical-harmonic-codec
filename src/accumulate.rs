//! Step 2: project the per-`m` Fourier coefficients onto the
//! fully-normalised associated Legendre basis to produce the SH
//! coefficients `a_lm`.

use crate::legendre::legendre_norm;
use crate::vec3::Vec3;
use rayon::prelude::*;

pub const SQRT2: f64 = std::f64::consts::SQRT_2;

/// Step 2: Legendre accumulation → flat `a_lm` of length `(lmax+1) * aw`,
/// indexed `l * aw + (m + lmax)`. Returns `(a_lm, aw)` where `aw = 2*lmax+1`.
#[allow(clippy::too_many_arguments)]
pub fn run(
    a_buf: &[Vec3],
    b_buf: &[Vec3],
    h: usize,
    lmax: usize,
    costheta: &[f64],
    sintheta: &[f64],
    dtheta: f64,
) -> (Vec<Vec3>, usize) {
    let aw = 2 * lmax + 1;
    let make_acc = || vec![Vec3::ZERO; (lmax + 1) * aw];
    let nthreads = rayon::current_num_threads();

    let chunk = h.div_ceil(nthreads).max(1);

    let a = (0..h)
        .collect::<Vec<_>>()
        .par_chunks(chunk)
        .map(|rows| {
            let mut local = make_acc();
            let mut q = vec![vec![0.0f64; lmax + 1]; lmax + 1];

            for &y in rows {
                legendre_norm(lmax, costheta[y], sintheta[y], &mut q);
                let st = sintheta[y] * dtheta;

                let a0 = a_buf[y] * st;
                for (l, row) in q.iter().enumerate().take(lmax + 1) {
                    let ql0 = row[0];
                    let idx = l * aw + lmax;
                    local[idx].r += ql0 * a0.r;
                    local[idx].g += ql0 * a0.g;
                    local[idx].b += ql0 * a0.b;
                }

                for m in 1..=lmax {
                    let am = a_buf[m * h + y] * (st * SQRT2);
                    let bm = b_buf[m * h + y] * (st * SQRT2);
                    for (l, row) in q.iter().enumerate().take(lmax + 1).skip(m) {
                        let qm = row[m];
                        let idx_p = l * aw + m + lmax;
                        let idx_n = l * aw + lmax - m;
                        local[idx_p].r += qm * am.r;
                        local[idx_p].g += qm * am.g;
                        local[idx_p].b += qm * am.b;
                        local[idx_n].r += qm * bm.r;
                        local[idx_n].g += qm * bm.g;
                        local[idx_n].b += qm * bm.b;
                    }
                }
            }
            local
        })
        .reduce(make_acc, |mut acc, other| {
            for (a, b) in acc.iter_mut().zip(other) {
                *a += b;
            }
            acc
        });

    (a, aw)
}
