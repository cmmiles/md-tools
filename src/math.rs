//! Module to contain mathematical functions and structures.

use std::{
    fmt,
    f64::consts::FRAC_PI_2,
    ops::{Add, Sub, Mul, Div},
};
use num::complex::Complex;

mod consts;

mod coords;
pub use coords::*;

mod vector;
pub use vector::*;

mod adjacency_matrix;
pub use adjacency_matrix::*;

/// Calculates the angle between two [`Vector`](crate::math::Vector)s in radians.
/// Use [`f64::to_degrees()`](f64::to_degrees) to get degrees.
pub fn angle_between_vectors(u: &Vector, v: &Vector) -> f64 {
    (dot_product(u, v) / (u.radius() * v.radius())).acos()
}

/// Calculates the dot product of two [`Vector`](crate::math::Vector)s.
pub fn dot_product(u: &Vector, v: &Vector) -> f64 {
    u.x*v.x + u.y*v.y + u.z*v.z
}

/// Computes the spherical harmonic Y_l^m for a given [`Vector`](crate::math::Vector).
/// Only 3rd, 4th and 6th order implemented at present.
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