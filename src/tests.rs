use super::*;
use structure::{Atom, Universe};
use std::path::PathBuf;
use math;

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