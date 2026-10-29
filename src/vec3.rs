//! Three-channel sample value (RGB), used for both pixel data and SH
//! coefficients.

#[derive(Copy, Clone, Default, Debug, PartialEq)]
pub struct Vec3 {
    pub r: f64,
    pub g: f64,
    pub b: f64,
}

impl Vec3 {
    pub const ZERO: Vec3 = Vec3 {
        r: 0.0,
        g: 0.0,
        b: 0.0,
    };

    #[inline]
    pub fn dot(self, o: Self) -> f64 {
        self.r * o.r + self.g * o.g + self.b * o.b
    }

    #[inline]
    pub fn norm_sq(self) -> f64 {
        self.dot(self)
    }
}

impl std::ops::AddAssign for Vec3 {
    #[inline]
    fn add_assign(&mut self, o: Self) {
        self.r += o.r;
        self.g += o.g;
        self.b += o.b;
    }
}

impl std::ops::Mul<f64> for Vec3 {
    type Output = Vec3;
    #[inline]
    fn mul(self, s: f64) -> Vec3 {
        Vec3 {
            r: self.r * s,
            g: self.g * s,
            b: self.b * s,
        }
    }
}
