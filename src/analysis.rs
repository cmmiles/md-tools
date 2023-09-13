//! Module for analysis functions.

use crate::structure::*;
use crate::math::*;
use std::f64::consts::PI;

mod steinhardt;

/// Calculates the water dipoles for a given universe at its current Frame.
pub fn water_dipole(molecules: &Vec<Molecule>, ref_axis: &Vector, opt_pbc: &Option<[f32;3]>)
    -> Result<(Vec<u32>, Vec<f64>), &'static str>
{
    let n = molecules.len();
    let mut frame_indices: Vec<u32> = Vec::with_capacity(n);
    let mut frame_output: Vec<f64> = Vec::with_capacity(n);

    for molecule in molecules.iter() {
        let com_hydrogen = Point::centre(vec![molecule.atoms[1], molecule.atoms[2]], opt_pbc);
        let dipole_vector = molecule.atoms[0].vector_to_coord(&com_hydrogen, opt_pbc);
        let theta = angle_between_vectors(&dipole_vector, ref_axis).to_degrees();
        frame_indices.push(molecule.id);
        frame_output.push(theta);
    }
    Ok((frame_indices, frame_output))
}

/// Calculates the Steinhardt order parameters for a given universe at its current Frame.
pub fn steinhardt(l: i8, cutoff: &f64, molecules: &Vec<Molecule>, opt_pbc: &Option<[f32;3]>, opt_ndx: Option<&Vec<usize>>)
-> (Vec<u32>, Vec<f64>) {
    let n = molecules.len();
    let frame_indices = match &opt_ndx {
        None => {
            let mut frame_indices: Vec<u32> = Vec::with_capacity(n);
            for molecule in molecules.iter() { frame_indices.push(molecule.id); }
            frame_indices
        }
        Some(ndx) => {
            let mut frame_indices: Vec<u32> = Vec::with_capacity(n);
            for i in ndx.iter() { frame_indices.push(molecules[*i].id); }
            frame_indices
        }
    };

    let mut frame_output: Vec<f64> = Vec::with_capacity(n);
    for qlm in steinhardt::qlm(&l, cutoff, molecules, opt_pbc, opt_ndx).into_iter() {
        let sum: f64 = qlm.iter().map(|x| x.norm_sqr()).sum();
        frame_output.push(((4.0*PI)/((2*l+1) as f64) * sum).sqrt());
    }

    (frame_indices, frame_output)
}

/// Calculates local version of Steinhardt Q6 parameter,
/// as defined by [Lechner and Dellago](https://doi.org/10.1063/1.2977970).
pub fn local_steinhardt(l: i8, cutoff: &f64, molecules: &Vec<Molecule>, opt_pbc: &Option<[f32;3]>, opt_ndx: Option<&Vec<usize>>)
-> (Vec<u32>, Vec<f64>) {
    let n = molecules.len();
    let frame_indices = match &opt_ndx {
        None => {
            let mut frame_indices: Vec<u32> = Vec::with_capacity(n);
            for molecule in molecules.iter() { frame_indices.push(molecule.id); }
            frame_indices
        }
        Some(ndx) => {
            let mut frame_indices: Vec<u32> = Vec::with_capacity(n);
            for i in ndx.iter() { frame_indices.push(molecules[*i].id); }
            frame_indices
        }
    };

    let mut frame_output: Vec<f64> = Vec::with_capacity(n);

    let qlm_vec = steinhardt::local_qlm(&l, cutoff, &molecules, &opt_pbc, opt_ndx);
    for qlm in qlm_vec.iter() {
        let sum: f64 = qlm.iter().map(|x| x.norm_sqr()).sum();
        frame_output.push(((4.0*PI)/((2*l+1) as f64) * sum).sqrt());
    }

    (frame_indices, frame_output)
}

/// Calculates local version of Steinhardt Q6 parameter,
/// as defined by [Li et al.](https://doi.org/10.1039/C1CP22167A).
pub fn local_steinhardt_2(l: i8, cutoff: &f64, molecules: &Vec<Molecule>, opt_pbc: &Option<[f32;3]>, opt_ndx: Option<&Vec<usize>>)
-> (Vec<u32>, Vec<f64>) {
    let n = molecules.len();
    let frame_indices = match &opt_ndx {
        None => {
            let mut frame_indices: Vec<u32> = Vec::with_capacity(n);
            for molecule in molecules.iter() { frame_indices.push(molecule.id); }
            frame_indices
        }
        Some(ndx) => {
            let mut frame_indices: Vec<u32> = Vec::with_capacity(n);
            for i in ndx.iter() { frame_indices.push(molecules[*i].id); }
            frame_indices
        }
    };

    let cutoff_sq = cutoff*cutoff;
    let mut frame_output: Vec<f64> = Vec::with_capacity(n);
    let qlm_vec = steinhardt::qlm(&l, cutoff, molecules, &opt_pbc, None);

    match &opt_ndx {
        None => {
            for (i, molecule) in molecules.iter().enumerate() {
                let mut n_neighbours: u8 = 0;
                let mut sum = 0.0;
                let qi_norm = qlm_vec[i].iter().map(|q| q.norm_sqr()).sum::<f64>().sqrt();
                for (j, neighbour) in molecules.iter().enumerate() {
                    if i != j && molecule.dsq(neighbour, &opt_pbc) <= cutoff_sq {
                        let mut pair_sum = 0.0;
                        n_neighbours += 1;
                        for m in 0..(2*l+1) as usize { pair_sum += (qlm_vec[i][m] * qlm_vec[j][m].conj()).re; }
                        let qj_norm = qlm_vec[j].iter().map(|q| q.norm_sqr()).sum::<f64>().sqrt();
                        sum += pair_sum / (qi_norm * qj_norm);
                    }
                }
                
                frame_output.push(sum / (n_neighbours as f64));
            }
        }
        Some(ndx) => {
            for i in ndx.iter() {
                let mut n_neighbours: u8 = 0;
                let mut sum = 0.0;
                let qi_norm = qlm_vec[*i].iter().map(|q| q.norm_sqr()).sum::<f64>().sqrt();
                for (j, neighbour) in molecules.iter().enumerate() {
                    if *i != j && molecules[*i].dsq(neighbour, &opt_pbc) <= cutoff_sq {
                        let mut pair_sum = 0.0;
                        n_neighbours += 1;
                        for m in 0..(2*l+1) as usize { pair_sum += (qlm_vec[*i][m] * qlm_vec[j][m].conj()).re; }
                        let qj_norm = qlm_vec[j].iter().map(|q| q.norm_sqr()).sum::<f64>().sqrt();
                        sum += pair_sum / (qi_norm * qj_norm);
                    }
                }
                
                frame_output.push(sum / (n_neighbours as f64));
            }
        }
    };

    (frame_indices, frame_output)
}

/// Computes the number of hydrogen bonds for a list of water molecules
/// (first atom must be oxygen, second and third atoms must be hydrogen)
/// given O--O distance cutoff `r_cut` \[nm\] and O--H--O angle cutoff `a_cut` [degrees].
pub fn get_n_hbonds<'a>(
    molecules: &Vec<Molecule<'a>>,
    r_cut: f64,
    a_cut: f64,
    opt_pbc:&Option<[f32;3]>,
) -> Vec<u8> {
    let r_cut_sq = r_cut*r_cut;
    let mut n_hbonds: Vec<u8> = (0..molecules.len()).map(|_| 0).collect();
    for (i, mol_i) in molecules.iter().enumerate() {
        for (j, mol_j) in molecules.iter().enumerate().skip(i+1) {
            if check_hbond_exists(&mol_i, &mol_j, r_cut_sq, a_cut, opt_pbc) {
                n_hbonds[i] += 1;
                n_hbonds[j] += 1;
            }
        }
    }

    n_hbonds
}

fn check_hbond_exists(
    mol1: &Molecule,
    mol2: &Molecule,
    r_cut_sq: f64,
    a_cut_deg: f64,
    opt_pbc:&Option<[f32;3]>,
) -> bool {
    if mol1.atoms[0].dsq(mol2.atoms[0], opt_pbc) <= r_cut_sq {
        let u = mol1.atoms[1].vector_to_coord(mol1.atoms[0], opt_pbc);
        let v = mol1.atoms[1].vector_to_coord(mol2.atoms[0], opt_pbc);
        if angle_between_vectors(&u, &v).to_degrees() >= a_cut_deg { return true; }
        let u = mol1.atoms[2].vector_to_coord(mol1.atoms[0], opt_pbc);
        let v = mol1.atoms[2].vector_to_coord(mol2.atoms[0], opt_pbc);
        if angle_between_vectors(&u, &v).to_degrees() >= a_cut_deg { return true; }
        let u = mol2.atoms[1].vector_to_coord(mol1.atoms[0], opt_pbc);
        let v = mol2.atoms[1].vector_to_coord(mol2.atoms[0], opt_pbc);
        if angle_between_vectors(&u, &v).to_degrees() >= a_cut_deg { return true; }
        let u = mol2.atoms[2].vector_to_coord(mol1.atoms[0], opt_pbc);
        let v = mol2.atoms[2].vector_to_coord(mol2.atoms[0], opt_pbc);
        if angle_between_vectors(&u, &v).to_degrees() >= a_cut_deg { return true; }
    }
    false
}