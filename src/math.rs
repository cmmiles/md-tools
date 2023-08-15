use std::fmt;
use std::f32::consts::{PI, TAU};
use num::complex::Complex;

/// Calculates the angle between two Vectors in radians. Use `x.to_degrees()` to get degrees.
pub fn angle_between_vectors(u: &Vector, v: &Vector) -> f32 {
    (dot_product(u, v) / (u.abs() * v.abs())).acos()
}

pub fn dot_product(u: &Vector, v: &Vector) -> f32 {
    u.x*v.x + u.y*v.y + u.z*v.z
}

#[derive(Debug)]
pub struct Vector {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    theta: Option<f32>,
    phi: Option<f32>,
}

impl fmt::Display for Vector {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Vector ({}, {}, {})", self.x, self.y, self.z)
    }
}

impl Vector {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x, y, z, theta: None, phi: None }
    }

    pub fn x() -> Self { Self { x: 1.0, y: 0.0, z: 0.0, theta: Some(0.0), phi: Some(0.0) } }
    pub fn y() -> Self { Self { x: 0.0, y: 1.0, z: 0.0, theta: Some(0.0), phi: Some(0.0) } }
    pub fn z() -> Self { Self { x: 0.0, y: 0.0, z: 1.0, theta: Some(0.0), phi: Some(0.0) } }

    pub fn dsq(&self) -> f32 {
        self.x*self.x + self.y*self.y + self.z*self.z
    }

    pub fn abs(&self) -> f32 { self.dsq().sqrt() }

    pub fn theta(&mut self) -> f32 {
        match self.theta {
            Some(theta) => theta,
            None => {
                let theta = angle_between_vectors(self, &Vector::z());
                self.theta = Some(theta);
                theta
            }
        }
    }

    pub fn phi(&mut self) -> f32 {
        match self.phi {
            Some(phi) => phi,
            None => {
                let proj_vector = Self { x: self.x, y: self.y, z: 0.0, theta: None, phi: None };
                let phi = match self.y >= 0.0 {
                    true => angle_between_vectors(&proj_vector, &Vector::x()),
                    false => TAU - angle_between_vectors(&proj_vector, &Vector::x()),
                };
                self.phi = Some(phi);
                phi
            }
        }
    }
}

pub fn spherical_harmonic(l: &i8, m: &i8, r: &mut Vector) -> Result<Complex<f32>, &'static str> {
    let theta = r.theta();
    let phi = r.phi();
    let sh_r = match (l, m) {
        (3, -3) => 0.125 * (35.0/PI).sqrt() * theta.sin().powi(3),
        (3, -2) => 0.25 * (105.0/TAU).sqrt() * theta.sin().powi(2) * theta.cos(),
        (3, -1) => 0.125 * (21.0/PI).sqrt() * theta.sin() * (5.0*(theta.cos().powi(2)) - 1.0),
        (3, 0) =>  0.25 * (7.0/PI).sqrt() * (5.0*theta.cos().powi(3) - 3.0*theta.cos()),
        (3, 1) => -0.125 * (21.0/PI).sqrt() * theta.sin() * (5.0*(theta.cos().powi(2)) - 1.0),
        (3, 2) =>  0.25 * (105.0/TAU).sqrt() * theta.sin().powi(2) * theta.cos(),
        (3, 3) => -0.125 * (35.0/PI).sqrt() * theta.sin().powi(3),
        _ => return Err("non-implemented spherical harmonic required"),
    };
    Ok(Complex::from_polar(sh_r, (*m as f32)*phi))
}