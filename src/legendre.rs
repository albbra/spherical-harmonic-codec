//! Fully-normalised associated Legendre functions, evaluated by stable
//! three-term recurrence (no Condon–Shortley phase).

pub const INV_SQRT_4PI: f64 = 0.282_094_791_773_878_14;

/// Fill `q[l][m]` (0 <= m <= l <= lmax) with fully-normalised associated
/// Legendre functions for given cos(theta)=x, sin(theta)=sth.
#[allow(clippy::needless_range_loop)]
pub fn legendre_norm(lmax: usize, x: f64, sth: f64, q: &mut [Vec<f64>]) {
    q[0][0] = INV_SQRT_4PI;

    for m in 1..=lmax {
        let mf = m as f64;
        q[m][m] = ((2.0 * mf + 1.0) / (2.0 * mf)).sqrt() * sth * q[m - 1][m - 1];
    }

    for m in 0..lmax {
        let mf = m as f64;
        q[m + 1][m] = (2.0 * mf + 3.0).sqrt() * x * q[m][m];
    }

    for m in 0..=lmax {
        let mf = m as f64;
        for l in (m + 2)..=lmax {
            let lf = l as f64;
            let aa = (((2.0 * lf - 1.0) * (2.0 * lf + 1.0)) / ((lf - mf) * (lf + mf))).sqrt();
            let bb = (((2.0 * lf + 1.0) * (lf + mf - 1.0) * (lf - mf - 1.0))
                / ((2.0 * lf - 3.0) * (lf - mf) * (lf + mf)))
                .sqrt();
            q[l][m] = aa * x * q[l - 1][m] - bb * q[l - 2][m];
        }
    }
}
