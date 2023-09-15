//! Handles subtasks for Steinhardt parameter calculations.

use num::complex::Complex;
use std::{thread, sync::mpsc};

use crate::structure::*;
use crate::math::*;

/// Calculates q_l^m parameters for the whole system at the current frame
pub fn qlm(
    &l: &i8,
    cutoff: &f64,
    molecules: &Vec<Molecule>,
    opt_pbc: &Option<[f32;3]>,
    opt_ndx: Option<&Vec<usize>>,
    min_coord_number: usize,
    max_coord_number: usize,
) -> Vec<Vec<Complex<f64>>>
{
    let nt = 1; // number of threads to use (experimental)

    // If no indices are specified, we calculate for all molecules
    // Otherwise, we calculate for only the specified molecules (but still need the other molecules for the calculation)
    let filtered_molecules = &super::filter_indices(molecules, opt_ndx);

    let n = filtered_molecules.len();

    let mut results: Vec<Vec<Complex<f64>>> = Vec::with_capacity(n);

    thread::scope(|s| {
        let mut receivers = Vec::with_capacity(nt);
        for i in 0..nt {
            let (tx, rx) = mpsc::channel();
            receivers.push(rx);
            s.spawn(move || {
                tx.send(thread_qlm(
                    &l,
                    cutoff,
                    &filtered_molecules[(i*n/nt)..((i+1)*n/nt)],
                    &molecules,
                    opt_pbc,
                    min_coord_number,
                    max_coord_number
                )).unwrap();
            });
        }
        for rx in receivers.iter() {
            let mut thread_results = rx.recv().unwrap();
            results.append(&mut thread_results)
        }
    });

    results
}

/// Calculates q_l^m parameters for a slice of the system at the current frame
/// Requires a slice of the molcules to calculate for, and also the full molecule list for the calculations
fn thread_qlm(
    &l: &i8,
    cutoff: &f64,
    mol_slice: &[Molecule],
    molecules: &Vec<Molecule>,
    opt_pbc: &Option<[f32;3]>,
    min_coord_number: usize,
    max_coord_number: usize,
) -> Vec<Vec<Complex<f64>>>
{
    let mut results: Vec<Vec<Complex<f64>>> = Vec::with_capacity(mol_slice.len());
    for molecule in mol_slice.iter() {
        let neighbours = molecule.coord_shell(&molecules, cutoff, opt_pbc);
        let n_neighbours = neighbours.len();

        // Only include molecules with the right sized coordination shell
        // Min should be at least 2 -- if n_neighbours==1, this will output 1.0 for all l,m
        if n_neighbours >= min_coord_number && n_neighbours <= max_coord_number {
            // Instantiate empty 2d vec for the spherical harmonics
            let mut spherical_harmonics: Vec<Vec<Complex<f64>>> = (0..(2*l + 1))
                .map(|_| Vec::with_capacity(n_neighbours)).collect();

            // Create iterator of vectors to each of the first coordination shell molecules
            // compute all the l-th order spherical harmonics for each such vector
            let r_iter = neighbours.iter().map(|x| molecule.atoms[0].vector_to_coord(x.atoms[0], opt_pbc));
            for mut r in r_iter {
                for m in -l..=l {
                    spherical_harmonics[(m + l) as usize].push(spherical_harmonic(&l, &m, &mut r)
                        .expect("Spherical harmonic error"));
                }
            }

            // Take average of spherical harmonics to give qlm
            let qlm_atom = spherical_harmonics.into_iter()
                .map(|sh_vec| sh_vec.iter().sum::<Complex<f64>>().unscale(n_neighbours as f64)).collect();

            results.push(qlm_atom);

        } else {
            // If n_neighbours not within chosen range, output qlm = NaN for all m
            results.push((-l..l).map(|_| Complex::new(f64::NAN, f64::NAN)).collect());
        }
    }
    results
}

/// Calculates averaged q_l^m parameters over first coord shell for the whole system at the current frame
pub fn local_qlm(
    &l: &i8,
    cutoff: &f64,
    molecules: &Vec<Molecule>,
    opt_pbc: &Option<[f32;3]>,
    opt_ndx: Option<&Vec<usize>>,
    min_coord_number: usize,
    max_coord_number: usize,
) -> Vec<Vec<Complex<f64>>>
{
    let n = molecules.len();
    let mut results: Vec<Vec<Complex<f64>>> = Vec::with_capacity(n);
    let cutoff_sq = cutoff*cutoff;

    // First calculate regular qlm values for all molecules
    let qlm_vec = qlm(&l, cutoff, molecules, opt_pbc, None, 2, 8);

    // If no indices are specified, we calculate for all molecules
    // Otherwise, we calculate for only the specified molecules (but still need the other molecules for the calculation)
    let filtered_molecules = &super::filter_indices(molecules, opt_ndx);

    // For each of the molecules we are calculating for, average the regular qlm across the molecule *and* its first
    // coordination shell
    for (i, molecule) in filtered_molecules.iter().enumerate() {
        let mut lqlm_atom = qlm_vec[i].clone();
        let mut n_neighbours: usize = 1; // Starts at 1 to include molecule i
        for (j, neighbour) in molecules.iter().enumerate() {
            if i != j && molecule.dsq(neighbour, &opt_pbc) <= cutoff_sq {
                n_neighbours += 1;
                for k in 0..(2*l as usize + 1) {
                    lqlm_atom[k] += qlm_vec[j][k];
                }
            }
        }
        
        if n_neighbours >= min_coord_number + 1 && n_neighbours <= max_coord_number + 1 {
            // Only include molecules with the right sized coordination shell
            results.push(lqlm_atom.into_iter().map(|q| q/(n_neighbours as f64)).collect());
        } else {
            // If n_neighbours not within chosen range, output local_qlm = NaN for all m
            results.push((-l..l).map(|_| Complex::new(f64::NAN, f64::NAN)).collect());
        }
    };

    results
}