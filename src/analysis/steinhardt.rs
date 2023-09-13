use num::complex::Complex;
use std::{thread, sync::mpsc};

use crate::structure::*;
use crate::math::*;

/// Calculates q_l^m parameters for the whole system at the current frame
pub fn qlm(&l: &i8, cutoff: &f64, molecules: &Vec<Molecule>, &opt_pbc: &Option<[f32;3]>, opt_ndx: Option<&Vec<usize>>)
-> Vec<Vec<Complex<f64>>>
{
    let nt = 1;

    let filtered_molecules = &match &opt_ndx {
        None => molecules.clone(),
        Some(ndx) => ndx.into_iter().map(|i| molecules[*i].clone()).collect(),
    };

    let n = filtered_molecules.len();

    let mut results: Vec<Vec<Complex<f64>>> = Vec::with_capacity(n);

    thread::scope(|s| {
        let mut receivers = Vec::with_capacity(nt);
        for i in 0..nt {
            let (tx, rx) = mpsc::channel();
            receivers.push(rx);
            s.spawn(move || {
                tx.send(thread_qlm(&l, cutoff, &filtered_molecules[(i*n/nt)..((i+1)*n/nt)], &molecules, &opt_pbc)).unwrap();
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
fn thread_qlm(&l: &i8, cutoff: &f64, mol_slice: &[Molecule], molecules: &Vec<Molecule>, opt_pbc: &Option<[f32;3]>)
    -> Vec<Vec<Complex<f64>>>
{
    let mut results: Vec<Vec<Complex<f64>>> = Vec::with_capacity(mol_slice.len());
    for molecule in mol_slice.iter() {
        let neighbours = molecule.coord_shell(&molecules, cutoff, opt_pbc);
        let n_neighbours = neighbours.len();
        let r_iter = neighbours.iter().map(|x| molecule.atoms[0].vector_to_coord(x.atoms[0], opt_pbc));
        let mut spherical_harmonics: Vec<Vec<Complex<f64>>> = (0..(2*l + 1)).map(|_| Vec::with_capacity(n_neighbours)).collect();
        for mut r in r_iter {
            for m in -l..=l {
                spherical_harmonics[(m + l) as usize].push(spherical_harmonic(&l, &m, &mut r).expect("Spherical harmonic error"));
            }
        }
        let mut qlm_atom = Vec::with_capacity((2*l + 1) as usize);
        for sh_vec in spherical_harmonics.iter() {
            qlm_atom.push(sh_vec.iter().sum::<Complex<f64>>().unscale(n_neighbours as f64));
        }
        results.push(qlm_atom);
    }
    results
}

/// Calculates averaged q_l^m parameters over first coord shell for the whole system at the current frame
pub fn local_qlm(&l: &i8, cutoff: &f64, molecules: &Vec<Molecule>, &opt_pbc: &Option<[f32;3]>, opt_ndx: Option<&Vec<usize>>)
-> Vec<Vec<Complex<f64>>>
{
    let n = molecules.len();
    let mut results: Vec<Vec<Complex<f64>>> = Vec::with_capacity(n);
    let qlm_vec = qlm(&l, cutoff, molecules, &opt_pbc, None);
    let cutoff_sq = cutoff*cutoff;

    match &opt_ndx {
        None => {
            for (i, molecule) in molecules.iter().enumerate() {
                let mut lqlm_atom = qlm_vec[i].clone();
                let mut n_neighbours: f64 = 1.0; // Starts at 1 to include molecule i
                for (j, neighbour) in molecules.iter().enumerate() {
                    if i != j && molecule.dsq(neighbour, &opt_pbc) <= cutoff_sq {
                        n_neighbours += 1.0;
                        for k in 0..(2*l as usize + 1) {
                            lqlm_atom[k] += qlm_vec[j][k];
                        }
                    }
                }
                results.push(lqlm_atom.into_iter().map(|q| q/n_neighbours).collect());
            }
        }
        Some(ndx) => {
            for i in ndx.iter() {
                let mut lqlm_atom = qlm_vec[*i].clone();
                let mut n_neighbours: f64 = 1.0; // Starts at 1 to include molecule i
                for (j, neighbour) in molecules.iter().enumerate() {
                    if *i != j && molecules[*i].dsq(neighbour, &opt_pbc) <= cutoff_sq {
                        n_neighbours += 1.0;
                        for k in 0..(2*l as usize + 1) {
                            lqlm_atom[k] += qlm_vec[j][k];
                        }
                    }
                }
                results.push(lqlm_atom.into_iter().map(|q| q/n_neighbours).collect());
            }
        }
    };

    results
}