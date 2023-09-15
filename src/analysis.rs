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

/// Calculates the [Steinhardt bond order parameters](https://doi.org/10.1103/PhysRevB.28.784)[^cite]
/// for a given universe at its current Frame.
///
/// [^cite]: [P. J. Steinhardt, D. R. Nelson and M. Ronchetti, *Phys. Rev. B*, 1983, **28**, 784--805](https://doi.org/10.1103/PhysRevB.28.784).
pub fn steinhardt(l: i8, cutoff: &f64, molecules: &Vec<Molecule>, opt_pbc: &Option<[f32;3]>, opt_ndx: Option<&Vec<usize>>)
-> (Vec<u32>, Vec<f64>) {
    let (min_coord_number, max_coord_number) = (4, 4);
    let frame_indices = filter_indices(molecules, opt_ndx).iter().map(|mol| mol.id).collect();
    let frame_output = steinhardt::qlm(&l, cutoff, molecules, opt_pbc, opt_ndx, min_coord_number, max_coord_number)
        .into_iter().map(|qlm| qlm.iter().map(|x| x.norm_sqr()).sum::<f64>())
        .map(|sum| ((4.0*PI)/((2*l+1) as f64) * sum).sqrt()).collect();

    (frame_indices, frame_output)
}

/// Calculates local version of Steinhardt Q6 parameter,
/// as defined by [Lechner and Dellago](https://doi.org/10.1063/1.2977970).[^cite]
///
/// [^cite]: [W. Lechner and C. Dellago, *J. Chem. Phys.*, 2008, **129**, 114707.](https://doi.org/10.1063/1.2977970)
pub fn local_steinhardt(l: i8, cutoff: &f64, molecules: &Vec<Molecule>, opt_pbc: &Option<[f32;3]>, opt_ndx: Option<&Vec<usize>>)
-> (Vec<u32>, Vec<f64>) {
    let (min_coord_number, max_coord_number) = (4, 4);
    let frame_indices = filter_indices(molecules, opt_ndx).iter().map(|mol| mol.id).collect();
    let frame_output = steinhardt::local_qlm(&l, cutoff, molecules, opt_pbc, opt_ndx, min_coord_number, max_coord_number)
        .into_iter().map(|qlm| qlm.iter().map(|x| x.norm_sqr()).sum::<f64>())
        .map(|sum| ((4.0*PI)/((2*l+1) as f64) * sum).sqrt()).collect();

    (frame_indices, frame_output)
}

/// Calculates local version of Steinhardt Q6 parameter,
/// as defined by [Li et al.](https://doi.org/10.1039/C1CP22167A).[^cite]
///
/// [^cite]: [T. Li, D. Donadio, G. Russo and G. Galli, *Phys. Chem. Chem. Phys.*, 2011, **13**, 19807--19813.](https://doi.org/10.1063/1.2977970)
pub fn local_steinhardt_2(l: i8, cutoff: &f64, molecules: &Vec<Molecule>, opt_pbc: &Option<[f32;3]>, opt_ndx: Option<&Vec<usize>>)
-> (Vec<u32>, Vec<f64>) {
    let (min_coord_number, max_coord_number) = (4, 4);
    let filtered_molecules = filter_indices(molecules, opt_ndx);
    let frame_indices = filtered_molecules.iter().map(|mol| mol.id).collect();

    let frame_output = filtered_molecules.iter().enumerate().map(|(i, molecule)| {
        let cutoff_sq = cutoff*cutoff;
        let mut n_neighbours: usize = 0;
        let mut sum = 0.0;
        
        // Calculate qlm for every molecule
        let qlm_vec = steinhardt::qlm(&l, cutoff, molecules, &opt_pbc, None, 2, 8);

        // Calculate norm of all qlms for molecule i (molecule)
        let qi_norm = qlm_vec[i].iter().map(|q| q.norm_sqr()).sum::<f64>().sqrt();

        for (j, neighbour) in molecules.iter().enumerate() {
            if i != j && molecule.dsq(neighbour, &opt_pbc) <= cutoff_sq {
                let mut pair_sum = 0.0;
                n_neighbours += 1;
                for m in 0..(2*l+1) as usize { pair_sum += (qlm_vec[i][m] * qlm_vec[j][m].conj()).re; }
                // Calculate norm of all qlms for molecule j (neighbour)
                let qj_norm = qlm_vec[j].iter().map(|q| q.norm_sqr()).sum::<f64>().sqrt();
                sum += pair_sum / (qi_norm * qj_norm);
            }
        }
        // As with other parameters, output NaN if coordination number not within range
        if n_neighbours >= min_coord_number && n_neighbours <= max_coord_number { sum / (n_neighbours as f64) }
        else { f64::NAN }
    }).collect();


    (frame_indices, frame_output)
}

/// Computes the number of hydrogen bonds for a list of water molecules
/// (first atom must be oxygen, second and third atoms must be hydrogen)
/// given O--O distance cutoff `r_cut` \[nm\] and O--H--O angle cutoff `a_cut` \[degrees\].
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

/// Checks whether there exists a hydrogen bond between two water molecules
/// Molecules must be structured [O, H, H (..)]
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

/// Takes a list of molecules and a list of indices (which should match the molecules list)
/// and outputs a list of just those molecules.
fn filter_indices<'a>(molecules: &'a Vec<Molecule>, opt_ndx: Option<&'a Vec<usize>>) -> Vec<Molecule<'a>> {
    match &opt_ndx {
        None => molecules.clone(),
        Some(ndx) => ndx.into_iter().map(|i| molecules[*i].clone()).collect(),
    }
}