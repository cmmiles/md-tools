use super::*;

/// An object which has coordinates, e.g. Atom, centre of mass
pub trait Coords {
    fn coords(&self) -> &[f32; 3];

    /// Calculates the vector to another Coords
    fn vector_to_coord(&self, other: &dyn Coords, opt_pbc: &Option<[f32; 3]>) -> Vector {
        let self_coords = self.coords();
        let other_coords = other.coords();
        let mut dx = self_coords[0] - other_coords[0];
        let mut dy = self_coords[1] - other_coords[1];
        let mut dz = self_coords[2] - other_coords[2];
        if let Some(box_dimensions) = opt_pbc {
            let box_dimensions = box_dimensions;
            while dx < -0.5 * box_dimensions[0] { dx += box_dimensions[0]; }
            while dx > 0.5 * box_dimensions[0] { dx -= box_dimensions[0]; }
            while dy < -0.5 * box_dimensions[1] { dy += box_dimensions[1]; }
            while dy > 0.5 * box_dimensions[1] { dy -= box_dimensions[1]; }
            while dz < -0.5 * box_dimensions[2] { dz += box_dimensions[2]; }
            while dz > 0.5 * box_dimensions[2] { dz -= box_dimensions[2]; }
        }
        Vector::new(dx, dy, dz)
    }

    /// Calculates the squared distance to another Coords
    fn dsq(&self, other: &dyn Coords, opt_pbc: &Option<[f32; 3]>) -> f32 {
        self.vector_to_coord(other, opt_pbc).dsq()
    }
}

/// Generic Point struct for positional data alone, e.g. centre of mass
pub struct Point([f32; 3]);

impl Coords for Point {
    fn coords(&self) -> &[f32; 3] { &self.0 }
}

impl Point {
    /// Finds the centre of a set of coordinates !not taking mass into account!
    /// There is potential for the centre to be incorrect if the points are spaced across more than 50% of the simulation box
    /// and periodic boundary conditions are in use.
    pub fn centre(points: Vec<&dyn Coords>, opt_pbc: &Option<[f32; 3]>) -> Point {
        let n = points.len() as f32;
        let mut coords_iter = points.iter().map(|point| point.coords());
        let [mut cx, mut cy, mut cz] = {
            let coords = coords_iter.next().unwrap();
            [coords[0], coords[1], coords[2]]
        };

        let [min, max] = match opt_pbc {
            Some(box_dimensions) => {
                let box_dimensions = box_dimensions;
                [
                    [
                        cx - 0.5 * box_dimensions[0],
                        cy - 0.5 * box_dimensions[1],
                        cz - 0.5 * box_dimensions[2],
                    ],
                    [
                        cx + 0.5 * box_dimensions[0],
                        cy + 0.5 * box_dimensions[1],
                        cz + 0.5 * box_dimensions[2]
                    ],
                ]
            }
            None => [
                [f32::MIN, f32::MIN, f32::MIN],
                [f32::MAX, f32::MAX, f32::MAX],
            ],
        };

        for coords in coords_iter {
            let [mut x, mut y, mut z] = [coords[0], coords[1], coords[2]];
            if let Some(box_dimensions) = opt_pbc {
                let box_dimensions = box_dimensions;
                while x < min[0] { x += box_dimensions[0]; }
                while x > max[0] { x -= box_dimensions[0]; }
                while y < min[1] { y += box_dimensions[1]; }
                while y > max[1] { y -= box_dimensions[1]; }
                while z < min[2] { z += box_dimensions[2]; }
                while z > max[2] { z -= box_dimensions[2]; }
            }
            cx += x;
            cy += y;
            cz += z;
        }
        Point ([cx/n, cy/n, cz/n])
    }
}