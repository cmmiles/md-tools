use crate::structure::*;
use crate::math::*;
use std::f32::consts::PI;

mod steinhardt;

/// Calculates the water dipoles for a given universe at its current Frame
pub fn water_dipole(molecules: &Vec<Molecule>, ref_axis: &Vector, opt_pbc: &Option<[f32;3]>)
    -> Result<(Vec<u32>, Vec<f32>), &'static str>
{
    let n = molecules.len();
    let mut frame_indices: Vec<u32> = Vec::with_capacity(n);
    let mut frame_output: Vec<f32> = Vec::with_capacity(n);

    for molecule in molecules.iter() {
        let com_hydrogen = Point::centre(vec![molecule.atoms[1], molecule.atoms[2]], opt_pbc);
        let dipole_vector = molecule.atoms[0].vector_to_coord(&com_hydrogen, opt_pbc);
        let theta = angle_between_vectors(&dipole_vector, ref_axis).to_degrees();
        frame_indices.push(molecule.id);
        frame_output.push(theta);
    }
    Ok((frame_indices, frame_output))
}

pub fn steinhardt(l: i8, molecules: &Vec<Molecule>, opt_pbc: &Option<[f32;3]>)
    -> (Vec<u32>, Vec<f32>)
{
    let n = molecules.len();
    let mut frame_indices: Vec<u32> = Vec::with_capacity(n);
    let mut frame_output: Vec<f32> = Vec::with_capacity(n);

    let qlm_vec = steinhardt::qlm(&l, &molecules, &opt_pbc);
    for molecule in molecules.iter() { frame_indices.push(molecule.id); }
    for qlm in qlm_vec.iter() {
        let sum: f32 = qlm.iter().map(|x| x.norm_sqr()).sum();
        frame_output.push(((4.0*PI)/((2*l+1) as f32) * sum).sqrt());
    }

    (frame_indices, frame_output)
}