use num::complex::Complex;
use std::{thread, sync::mpsc};

use crate::structure::*;
use crate::math::*;

pub fn qlm(&l: &i8, molecules: &Vec<Molecule>, &opt_pbc: &Option<[f32;3]>)
-> Vec<Vec<Complex<f32>>>
{
    let nt = 1;
    let n = molecules.len();
    let mut results: Vec<Vec<Complex<f32>>> = Vec::with_capacity(n);

    thread::scope(|s| {
        let mut receivers = Vec::with_capacity(nt);
        for i in 0..nt {
            let (tx, rx) = mpsc::channel();
            receivers.push(rx);
            s.spawn(move || { tx.send(thread_qlm(&l, &molecules[(i*n/nt)..((i+1)*n/nt)], &molecules, &opt_pbc)).unwrap(); });
        }
        for rx in receivers.iter() {
            let mut thread_results = rx.recv().unwrap();
            results.append(&mut thread_results)
        }
    });

    results
}

fn thread_qlm(&l: &i8, mol_slice: &[Molecule], molecules: &Vec<Molecule>, opt_pbc: &Option<[f32;3]>)
    -> Vec<Vec<Complex<f32>>>
{
    let mut results: Vec<Vec<Complex<f32>>> = Vec::with_capacity(mol_slice.len());
    for molecule in mol_slice.iter() {
        let neighbours = molecule.coord_shell(&molecules, 0.35, opt_pbc);
        let n_neighbours = neighbours.len();
        let r_iter = neighbours.iter().map(|x| molecule.atoms[0].vector_to_coord(x.atoms[0], opt_pbc));
        let mut spherical_harmonics: Vec<Vec<Complex<f32>>> = (0..(2*l + 1)).map(|_| Vec::with_capacity(n_neighbours)).collect();
        for mut r in r_iter {
            for m in -l..=l {
                spherical_harmonics[(m + l) as usize].push(spherical_harmonic(&l, &m, &mut r).expect("Spherical harmonic error"));
            }
        }
        let mut qlm_atom = Vec::with_capacity((2*l + 1) as usize);
        for sh_vec in spherical_harmonics.iter() {
            qlm_atom.push(sh_vec.iter().sum::<Complex<f32>>().unscale(n_neighbours as f32));
        }
        results.push(qlm_atom);
    }
    results
}