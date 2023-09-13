//! # [`md-tools`](crate)
//! A collection of analysis tools for molecular simulations in .xtc and .trr file formats.
//! 
//! ## How to run:
//! [`md-tools`](crate) can be run from the command line as follows:
//! ```text
//! md-tools steinhardt --q6 -s md-gro -t traj.xtc -start 0 --end 1000
//! ```
//! where `md-tools` is the executable. The first argument after `md-tools` is the analysis to be run, options:
//! * `nframes` for outputting the total number of frames in the trajectory,[^note1]
//! * `convert` for converting into .xtc .trr or a series of .gro files, this can be used to shorten the trajectory,
//! * `dipole` for calculating water dipole orientations with respect to the *z* axis,
//! * `steinhardt` for calculating Steinhardt bond order parameters,
//! * `cluster` for calculating the largest ice-like clusters in each frame (work in progress).
//!
//! ### Additional arguments
//! * `-s` or `--sfile` \[\<.gro\>\] structure file,
//! * `-t` or `--tfile` \[\<.xtc/.trr\>\] trajectory file,
//! * `-o` or `--output` output file,[^note2]
//! * `--start` \<int\> (0) analysis start time \[ps\],[^note3]
//! * `--end` \<int\> (+∞) analysis end time \[ps\],[^note3]
//! * `--stride` \<int\> (1) analysis time step \[ps\],[^note3] should be a multiple of the trajectory time step,
//! * `--no-pbc` disable periodic boundary conditions for analysis.
//! ### Arguments for `md-tools steinhardt`
//! * `--q3` enable 3rd order Steinhardt analysis (for use with `md-tools steinhardt`),
//! * `--q4` enable 4th order Steinhardt analysis (for use with `md-tools steinhardt`),
//! * `--q6` enable 6th order Steinhardt analysis (for use with `md-tools steinhardt`).
//! ### Arguments for `md-tools cluster`
//! * `--n_cls` \<int\> (1) number of clusters to output for `md-tools cluster`.
//! ### Unused arguments
//! * `--th` enable water dipole analysis (for future use, currently does nothing).
//!
//! [^note1]: `md-tools nframes` only accepts the arguments `-s`/`--sfile` and `-t`/`--tfile`.
//!
//! [^note2]: Depending on the analysis, this may be extended, e.g. `"md-tools.out" -> "md-tools.out.q6"`.
//! 
//! [^note3]: Same unit as in trajectory file -- likely picoseconds.

use structure::Universe;
use config::Task;
use math::Vector;

pub mod structure;
pub mod config;
pub mod analysis;
pub mod math;
pub mod output;

#[cfg(test)] mod tests;

/// Main function for running md-tools.
pub fn run(config: config::Config) -> Result<(), &'static str> {
    //println!("{:#?}", config);

    let mut u = Universe::new(config.sfile, config.tfile, config.pbc);
    //println!("{}", u);

    //let firstcs = u.first_coord_shell(0, "OW");
    //println!("{} {:?}", firstcs.len(), firstcs);

    //println!("{:#?}", u.get_molecules());

    /*let start = Instant::now();
    let duration = start.elapsed();

    println!("Time elapsed: {:?}", duration);*/
    
    match config.task {
        Task::NumFrames => Ok(println!("{}", u.get_nframes()?)),
        Task::Convert => u.convert(&config.mint, &config.maxt, &config.stride, &config.outfile),
        Task::Dipole => u.water_dipole("SOL", Vector::ez(), &config.mint, &config.maxt, &config.stride, &config.outfile),
        Task::Steinhardt => u.steinhardt("SOL", &config.op_list, &0.32, &config.mint, &config.maxt, &config.stride, &config.outfile),
        Task::Cluster => u.q6_clustering("SOL", &0.32, &config.mint, &config.maxt, &config.stride, &config.n_cls, &config.outfile),
    }
}