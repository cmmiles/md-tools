use super::*;

/// Three-dimensional (mathematical) vector.
/// Stored in cartesian coordinates but spherical coordinates θ and φ are also stored if they are required.
/// Spherical coordinate r can be obtained via [`Vector::radius()`](crate::math::Vector::radius).
/// Addition / subtraction of `Vector`s and multiplication / division by a scalar ([`f64`](f64)) are all implemented.
///
/// # Examples
///
/// ## Adding two `Vector`s together
/// ```
/// use md_tools::math::Vector;
///
/// let u = Vector::new(1.0, 2.0, 3.0);
/// let v = Vector::new(4.0, 5.0, 6.0);
///
/// assert_eq!(u + v, Vector::new(5.0, 7.0, 9.0));
/// ```
///
/// ## Multiplying a `Vector` by a scalar
/// ```
/// use md_tools::math::Vector;
///
/// let v = Vector::new(1.0, 2.0, 3.0);
///
/// assert_eq!(v * 1.5, Vector::new(1.5, 3.0, 4.5));
/// ```
#[derive(Debug, Clone)]
pub struct Vector {
    /// x-component of the `Vector`
    pub x: f64,
    /// y-component of the `Vector`
    pub y: f64,
    /// z-component of the `Vector`
    pub z: f64,
    /// Polar spherical coordinate θ ∈ [0, π]
    theta: Option<f64>,
    /// Asimuthal spherical coordinate φ ∈ (-π, π]
    phi: Option<f64>,
}

impl fmt::Display for Vector {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Vector ({}, {}, {})", self.x, self.y, self.z)
    }
}

impl PartialEq for Vector {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y && self.z == other.z
    }
}

impl Add for Vector {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
            theta: None, phi: None,
        }
    }
}

impl Sub for Vector {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
            theta: None, phi: None,
        }
    }
}

/// Multiplication by a scalar `(Vector * `[`f64`](f64)`) -> Vector`.
impl Mul<f64> for Vector {
    type Output = Self;
    
    fn mul(self, scalar: f64) -> Self {
        Self {
            x: self.x * scalar,
            y: self.y * scalar,
            z: self.z * scalar,
            theta: self.theta,
            phi: self.phi,
        }
    }
}

/// Division by a scalar `(Vector / `[`f64`](f64)`) -> Vector`.
impl Div<f64> for Vector {
    type Output = Self;
    
    fn div(self, scalar: f64) -> Self {
        Self {
            x: self.x / scalar,
            y: self.y / scalar,
            z: self.z / scalar,
            theta: self.theta,
            phi: self.phi,
        }
    }
}

impl Vector {
    /// Constructor function for `Vector`, accepting cartesian coordinates.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z, theta: None, phi: None }
    }

    /// Constructor for x-axis basis `Vector` (1, 0, 0).
    pub fn ex() -> Self { Self { x: 1.0, y: 0.0, z: 0.0, theta: Some(FRAC_PI_2), phi: Some(0.0) } }

    /// Constructor for y-axis basis `Vector` (0, 1, 0).
    pub fn ey() -> Self { Self { x: 0.0, y: 1.0, z: 0.0, theta: Some(FRAC_PI_2), phi: Some(FRAC_PI_2) } }

    /// Constructor for z-axis basis `Vector` (0, 0, 1).
    pub fn ez() -> Self { Self { x: 0.0, y: 0.0, z: 1.0, theta: Some(0.0), phi: Some(0.0) } }

    /// Calculates the square of the size (radius) of the `Vector`.
    pub fn rsq(&self) -> f64 {
        self.x*self.x + self.y*self.y + self.z*self.z
    }

    /// Calculates the size (radius) of the `Vector`.
    pub fn radius(&self) -> f64 { self.rsq().sqrt() }

    /// Returns the polar spherical coordinate θ ∈ [0, π], computing it if necessary and storing it for later use.
    pub fn theta(&mut self) -> f64 {
        match self.theta {
            Some(theta) => theta,
            None => {
                let theta = angle_between_vectors(self, &Vector::ez());
                self.theta = Some(theta);
                theta
            }
        }
    }

    /// Returns the azimuthal spherical coordinate φ ∈ (-π, π], computing it if necessary and storing it for later use.
    pub fn phi(&mut self) -> f64 {
        match self.phi {
            Some(phi) => phi,
            None => {
                let proj_vector = Self { x: self.x, y: self.y, z: 0.0, theta: None, phi: None };
                let phi = match self.y >= 0.0 {
                    true  =>   angle_between_vectors(&proj_vector, &Vector::ex()),
                    false => - angle_between_vectors(&proj_vector, &Vector::ex()),
                };
                self.phi = Some(phi);
                phi
            }
        }
    }
}