use std::fmt;
use num::complex::Complex;

mod consts;

/// Calculates the angle between two Vectors in radians. Use `x.to_degrees()` to get degrees.
pub fn angle_between_vectors(u: &Vector, v: &Vector) -> f64 {
    (dot_product(u, v) / (u.abs() * v.abs())).acos()
}

pub fn dot_product(u: &Vector, v: &Vector) -> f64 {
    u.x*v.x + u.y*v.y + u.z*v.z
}

#[derive(Debug)]
pub struct Vector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    theta: Option<f64>,
    phi: Option<f64>,
}

impl fmt::Display for Vector {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Vector ({}, {}, {})", self.x, self.y, self.z)
    }
}

impl Vector {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z, theta: None, phi: None }
    }

    pub fn x() -> Self { Self { x: 1.0, y: 0.0, z: 0.0, theta: Some(0.0), phi: Some(0.0) } }
    pub fn y() -> Self { Self { x: 0.0, y: 1.0, z: 0.0, theta: Some(0.0), phi: Some(0.0) } }
    pub fn z() -> Self { Self { x: 0.0, y: 0.0, z: 1.0, theta: Some(0.0), phi: Some(0.0) } }

    pub fn dsq(&self) -> f64 {
        self.x*self.x + self.y*self.y + self.z*self.z
    }

    pub fn abs(&self) -> f64 { self.dsq().sqrt() }

    pub fn theta(&mut self) -> f64 {
        match self.theta {
            Some(theta) => theta,
            None => {
                let theta = angle_between_vectors(self, &Vector::z());
                self.theta = Some(theta);
                theta
            }
        }
    }

    pub fn phi(&mut self) -> f64 {
        match self.phi {
            Some(phi) => phi,
            None => {
                let proj_vector = Self { x: self.x, y: self.y, z: 0.0, theta: None, phi: None };
                let phi = match self.y >= 0.0 {
                    true  =>   angle_between_vectors(&proj_vector, &Vector::x()),
                    false => - angle_between_vectors(&proj_vector, &Vector::x()),
                };
                self.phi = Some(phi);
                phi
            }
        }
    }
}

pub fn spherical_harmonic(l: &i8, m: &i8, r: &mut Vector) -> Result<Complex<f64>, &'static str> {
    let theta = r.theta();
    let phi = r.phi();
    let sh_r = match (l, m) {
        (3, -3) => - consts::SH3_3 * theta.sin().powi(3),
        (3, -2) =>   consts::SH3_2 * theta.sin().powi(2) * theta.cos(),
        (3, -1) => - consts::SH3_1 * theta.sin() * (5.0*theta.cos().powi(2) - 1.0),
        (3,  0) =>   consts::SH3_0 * (5.0*theta.cos().powi(3) - 3.0*theta.cos()),
        (3,  1) =>   consts::SH3_1 * theta.sin() * (5.0*theta.cos().powi(2) - 1.0),
        (3,  2) =>   consts::SH3_2 * theta.sin().powi(2) * theta.cos(),
        (3,  3) =>   consts::SH3_3 * theta.sin().powi(3),

        (4, -4) =>   consts::SH4_4 * theta.sin().powi(4),
        (4, -3) => - consts::SH4_3 * theta.sin().powi(3) * theta.cos(),
        (4, -2) =>   consts::SH4_2 * theta.sin().powi(2) * (7.0*theta.cos().powi(2) - 1.0),
        (4, -1) => - consts::SH4_1 * theta.sin() * (7.0*theta.cos().powi(3) - 3.0*theta.cos()),
        (4,  0) =>   consts::SH4_0 * (35.0*theta.cos().powi(4) - 30.0*theta.cos().powi(2) + 3.0),
        (4,  1) =>   consts::SH4_1 * theta.sin() * (7.0*theta.cos().powi(3) - 3.0*theta.cos()),
        (4,  2) =>   consts::SH4_2 * theta.sin().powi(2) * (7.0*theta.cos().powi(2) - 1.0),
        (4,  3) =>   consts::SH4_3 * theta.sin().powi(3) * theta.cos(),
        (4,  4) =>   consts::SH4_4 * theta.sin().powi(4),

        (6, -6) =>   consts::SH6_6 * theta.sin().powi(6),
        (6, -5) => - consts::SH6_5 * theta.sin().powi(5) * theta.cos(),
        (6, -4) =>   consts::SH6_4 * theta.sin().powi(4) * (11.0 * theta.cos().powi(2) - 1.0),
        (6, -3) => - consts::SH6_3 * theta.sin().powi(3) * (11.0 * theta.cos().powi(3) - 3.0*theta.cos()),
        (6, -2) =>   consts::SH6_2 * theta.sin().powi(2) * (33.0 * theta.cos().powi(4) - 18.0*theta.cos().powi(2) + 1.0),
        (6, -1) => - consts::SH6_1 * theta.sin() * (33.0*theta.cos().powi(5) - 30.0*theta.cos().powi(3) + 5.0*theta.cos()),
        (6,  0) =>   consts::SH6_0 * (231.0*theta.cos().powi(6) - 315.0*theta.cos().powi(4) + 105.0*theta.cos().powi(2) - 5.0),
        (6,  1) =>   consts::SH6_1 * theta.sin() * (33.0*theta.cos().powi(5) - 30.0*theta.cos().powi(3) + 5.0*theta.cos()),
        (6,  2) =>   consts::SH6_2 * theta.sin().powi(2) * (33.0 * theta.cos().powi(4) - 18.0*theta.cos().powi(2) + 1.0),
        (6,  3) =>   consts::SH6_3 * theta.sin().powi(3) * (11.0 * theta.cos().powi(3) - 3.0*theta.cos()),
        (6,  4) =>   consts::SH6_4 * theta.sin().powi(4) * (11.0 * theta.cos().powi(2) - 1.0),
        (6,  5) =>   consts::SH6_5 * theta.sin().powi(5) * theta.cos(),
        (6,  6) =>   consts::SH6_6 * theta.sin().powi(6),

        _ => return Err("non-implemented spherical harmonic required"),
    };
    Ok(Complex::from_polar(sh_r, (*m as f64)*phi))
}