use std::{
    fmt,
    rc::Rc,
};
use xdrfile::*;
use crate::math::*;

pub mod generic;
pub use generic::*;

pub mod atom;
pub use atom::*;

pub mod universe;
pub use universe::*;

/// Container for Trajectory Iterators. Currently only XTC/TRR supported but can be expanded in the future.
/// Additional trajectory formats will have to be converted into the `Frame` type from the `xdrfile` crate.
enum MultiTrajIter {
    XTC(TrajectoryIterator<XTCTrajectory>),
    TRR(TrajectoryIterator<TRRTrajectory>),
}

impl Iterator for MultiTrajIter {
    type Item = Result<Rc<Frame>>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::XTC(traj_iter) => traj_iter.next(),
            Self::TRR(traj_iter) => traj_iter.next(),
        }
    }
}