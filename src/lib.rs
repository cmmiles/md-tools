use structure::Universe;
use config::Task;
use math::Vector;

pub mod structure;
pub mod config;
pub mod analysis;
pub mod math;
pub mod output;

#[cfg(test)] mod tests;

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
        Task::Dipole => u.water_dipole("SOL", Vector::z(), &config.mint, &config.maxt, &config.stride, &config.outfile),
        Task::Steinhardt => u.steinhardt("SOL", &config.op_list, &config.mint, &config.maxt, &config.stride, &config.outfile),
    }
}