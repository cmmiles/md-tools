use super::*;
use structure::{Atom, Universe};
use std::path::PathBuf;
use math;
use utils;

// Test GRO file reading
#[test]
fn read_gro_test() {
    let u = Universe::new(Some(PathBuf::from("tests/read_gro_test.gro")), None, true);
    let atoms = vec![
        Atom::new(1, "OW", "SOL", 1, [0.181, 0.043, 0.027]),
        Atom::new(2, "HW1", "SOL", 1, [0.086, 0.054, 0.035]),
        Atom::new(3, "HW2", "SOL", 1, [0.209, 0.016, 0.115]),
    ];
    let pbc = [0.3, 0.2, 0.3];
    assert!(u.compare_gro(atoms, pbc));
}

// 3 invalid GRO files and one which is nonexistent
#[test]
#[should_panic]
fn invalid_gro_test_1() {
    Universe::new(Some(PathBuf::from("tests/invalid_gro_1.gro")), None, true);
}

#[test]
#[should_panic]
fn invalid_gro_test_2() {
    Universe::new(Some(PathBuf::from("tests/invalid_gro_2.gro")), None, true);
}

#[test]
#[should_panic]
fn invalid_gro_test_3() {
    Universe::new(Some(PathBuf::from("tests/invalid_gro_3.gro")), None, true);
}

#[test]
#[should_panic]
fn missing_gro_test() {
    Universe::new(Some(PathBuf::from("tests/missing.gro")), None, true);
}

// Test PDB to GRO conversion
#[test]
fn pdb_to_gro_line() {
    let pdb_line = "ATOM      1  OW  SOL X   1      29.390  24.890 110.650  0.00  0.00              ";
    let gro_line = "    1SOL     OW    1   2.939   2.489  11.065";
    assert_eq!(Atom::from_pdb(pdb_line).unwrap().to_gro(), gro_line);
}

// Test GRO to PDB conversion
#[test]
fn gro_to_pdb_line() {
    let pdb_line = "ATOM      1  OW  SOL X   1      29.390  24.890 110.650  0.00  0.00              ";
    let gro_line = "    1SOL     OW    1   2.939   2.489  11.065";
    assert_eq!(Atom::from_gro(gro_line).unwrap().to_pdb(), pdb_line);
}

// Testing writing a GRO file with indices avove the max of 99,999
#[test]
fn write_gro_oversize() {
    let oversize_atom = Atom::new(108_056, "OW", "SOL", 127_001, [1.0, 1.0, 1.0]);
    let gro_line = "27001SOL     OW 8056   1.000   1.000   1.000";
    assert_eq!(oversize_atom.to_gro(), gro_line);
}

// Testing writing a PDB file with indices avove the max of 99,999 and 9,999
#[test]
fn write_pdb_oversize() {
    let oversize_atom = Atom::new(108_056, "OW", "SOL", 127_001, [1.0, 1.0, 1.0]);
    let pdb_line = "ATOM   8056  OW  SOL X7001      10.000  10.000  10.000  0.00  0.00              ";
    assert_eq!(oversize_atom.to_pdb(), pdb_line);
}

// Testing math functions give correct outputs
#[test]
fn angle_between_perpendicular_vectors_test() {
    let u = Vector::new(0.0, 0.0, 1.0);
    let v = Vector::new(1.0, 0.0, 0.0);
    assert_eq!(math::angle_between_vectors(&u, &v).to_degrees(), 90.0);
}

#[test]
fn angle_between_parallel_vectors_test() {
    let u = Vector::new(89.0, 0.0, 0.0);
    let v = Vector::new(0.0001, 0.0, 0.0);
    assert_eq!(math::angle_between_vectors(&u, &v).to_degrees(), 0.0);
}

// Testing PDB atom name justification
#[test]
fn justify_atom_names_test() {
    assert_eq!(utils::justify_atom_name("OW"), String::from(" OW "));
    assert_eq!(utils::justify_atom_name("HW1"), String::from(" HW1"));
    assert_eq!(utils::justify_atom_name("NA"), String::from("NA  "));
    assert_eq!(utils::justify_atom_name("H10R"), String::from("H10R"));
}