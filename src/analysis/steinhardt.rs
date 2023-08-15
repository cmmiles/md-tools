use num::complex::Complex;

use crate::structure::*;
use crate::math::*;

pub fn qlm(&l: &i8, molecules: &Vec<Molecule>, opt_pbc: &Option<[f32;3]>)
-> Result<Vec<Vec<Complex<f32>>>, &'static str>
{
    let mut results: Vec<Vec<Complex<f32>>> = Vec::with_capacity(molecules.len());
    for molecule in molecules.iter() {
        let neighbours = molecule.coord_shell(molecules, 0.35, opt_pbc);
        let n_neighbours = neighbours.len();
        let r_iter = neighbours.iter().map(|x| molecule.atoms[0].vector_to_coord(x.atoms[0], opt_pbc));
        let mut spherical_harmonics: Vec<Vec<Complex<f32>>> = (0..(2*l + 1)).map(|_| Vec::with_capacity(n_neighbours)).collect();
        for mut r in r_iter {
            for m in -l..=l {
                spherical_harmonics[(m + l) as usize].push(spherical_harmonic(&l, &m, &mut r)?);
            }
        }
        let mut qlm_atom = Vec::with_capacity((2*l + 1) as usize);
        for sh_vec in spherical_harmonics.iter() {
            qlm_atom.push(sh_vec.iter().sum::<Complex<f32>>().unscale(n_neighbours as f32));
        }
        results.push(qlm_atom);
    }
    Ok(results)
}