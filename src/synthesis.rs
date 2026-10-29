//! Step 3: inverse transform from `a_lm` back to pixels at a given band `L`.

use rayon::prelude::*;

use crate::legendre::legendre_norm;
use crate::vec3::Vec3;

pub const SQRT2: f64 = std::f64::consts::SQRT_2;

/// Step 3: inverse transform at band `hit_l`.
#[allow(clippy::too_many_arguments)]
pub fn run(
    a_lm: &[Vec3],
    w: usize,
    h: usize,
    aw: usize,
    costheta: &[f64],
    sintheta: &[f64],
    cos_t: &[f64],
    sin_t: &[f64],
    stride: usize,
    hit_l: usize,
) -> Vec<Vec3> {
    let mut recon = vec![Vec3::ZERO; w * h];
    let lr = hit_l;
    let lmax = stride - 1;

    recon.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        let mut q = vec![vec![0.0f64; lr + 1]; lr + 1];
        legendre_norm(lr, costheta[y], sintheta[y], &mut q);

        let mut gp = vec![Vec3::ZERO; lr + 1];
        let mut gn = vec![Vec3::ZERO; lr + 1];

        for m in 0..=lr {
            let sm = if m == 0 { 1.0 } else { SQRT2 };
            let mut gpr = 0.0;
            let mut gpg = 0.0;
            let mut gpb = 0.0;
            let mut gnr = 0.0;
            let mut gng = 0.0;
            let mut gnb = 0.0;
            for l in m..=lr {
                let qm = q[l][m] * sm;
                let ap = a_lm[l * aw + m + lmax];
                gpr += ap.r * qm;
                gpg += ap.g * qm;
                gpb += ap.b * qm;
                if m > 0 {
                    let an = a_lm[l * aw + lmax - m];
                    gnr += an.r * qm;
                    gng += an.g * qm;
                    gnb += an.b * qm;
                }
            }
            gp[m] = Vec3 {
                r: gpr,
                g: gpg,
                b: gpb,
            };
            gn[m] = Vec3 {
                r: gnr,
                g: gng,
                b: gnb,
            };
        }

        for x in 0..w {
            let ct = &cos_t[x * stride..x * stride + lr + 1];
            let st = &sin_t[x * stride..x * stride + lr + 1];
            let mut acc = gp[0];
            for m in 1..=lr {
                acc.r += gp[m].r * ct[m] + gn[m].r * st[m];
                acc.g += gp[m].g * ct[m] + gn[m].g * st[m];
                acc.b += gp[m].b * ct[m] + gn[m].b * st[m];
            }
            row[x] = acc;
        }
    });

    recon
}
