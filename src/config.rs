//! Module for parsing input arguments.

use std::path::PathBuf;

/// Holds config setup for the CLI app.
#[derive(Debug)]
pub struct Config {
    /// Analysis procedure
    pub task: Task,
    /// Structure file
    pub sfile: Option<PathBuf>,
    /// Trajectory file
    pub tfile: Option<PathBuf>,
    /// Output file
    pub outfile: Option<PathBuf>,
    /// Order parameters to calculate for `md-tools steinhardt`
    pub op_list: Vec<OrderParameter>,
    /// Number of clusters to output for `md-tools cluster` (1)
    pub n_cls: usize,
    /// Start time for analysis (0)
    pub mint: u32,
    /// End time for analysis (+∞)
    pub maxt: u32,
    /// Time step for analysis (1)
    pub stride: u32,
    /// Whether to use periodic boundary conditions (yes)
    pub pbc: bool,
}

impl Config {
    /// Build the [`Config`](Config) from command line arguments.
    pub fn build(mut args: impl Iterator<Item = String>) -> Result<Self, &'static str> {
        args.next();

        let task = match args.next() {
            Some(task) => match &task[..] {
                "nframes" => Task::NumFrames,
                "convert" => Task::Convert,
                "dipole" => Task::Dipole,
                "steinhardt" => Task::Steinhardt,
                "cluster" => Task::Cluster,
                _ => return Err("invalid task"),
            }
            None => return Err("no arguments passed"),
        };

        let mut config = Self { task,
            sfile: None,
            tfile: None,
            outfile: None,
            op_list: Vec::new(),
            n_cls: 1,
            mint: 0,
            maxt: u32::MAX,
            stride: 1,
            pbc: true,
        };

        while let Some(result) = Self::read_next_arg(&mut args, &mut config) {
            if let Err(err) = result { return Err(err) }
        };

        Ok(config)
    }

    fn read_next_arg(args: &mut impl Iterator<Item = String>, config: &mut Self) -> Option<Result<(), &'static str>> {
        match &args.next()?[..] {
            "-s" | "--sfile" => match args.next() {
                Some(sfile) => config.sfile = Some(PathBuf::from(sfile)),
                None => return Some(Err("missing structure file declaration")),
            }
            "-t" | "--tfile" => match args.next() {
                Some(tfile) => config.tfile = Some(PathBuf::from(tfile)),
                None => return Some(Err("missing trajectory file declaration")),
            }
            "-o" | "--output" => match args.next() {
                Some(outfile) => config.outfile = Some(PathBuf::from(outfile)),
                None => return Some(Err("missing output file declaration")),
            }
            "--th" => if !config.op_list.contains(&OrderParameter::Theta) { config.op_list.push(OrderParameter::Theta) },
            "--q3" => if !config.op_list.contains(&OrderParameter::Q3) { config.op_list.push(OrderParameter::Q3) },
            "--q4" => if !config.op_list.contains(&OrderParameter::Q4) { config.op_list.push(OrderParameter::Q4) },
            "--q6" => if !config.op_list.contains(&OrderParameter::Q6) { config.op_list.push(OrderParameter::Q6) },
            "--n_cls" => match args.next() {
                Some(n_cls) => config.n_cls = match n_cls.parse() {
                    Ok(n_cls) => n_cls,
                    Err(_) => return Some(Err("invalid number of cluster"))
                },
                None => return Some(Err("missing number of clusters declaration")),
            }
            "--start" => match args.next() {
                Some(mint) => config.mint = match mint.parse() {
                    Ok(mint) => mint,
                    Err(_) => return Some(Err("invalid start time")),
                },
                None => return Some(Err("missing start time declaration")),
            }
            "--end" => match args.next() {
                Some(maxt) => config.maxt = match maxt.parse() {
                    Ok(maxt) => maxt,
                    Err(_) => return Some(Err("invalid end time"))
                },
                None => return Some(Err("missing start time declaration")),
            }
            "--stride" => match args.next() {
                Some(stride) => config.stride = match stride.parse() {
                    Ok(stride) => stride,
                    Err(_) => return Some(Err("invalid stride value")),
                },
                None => return Some(Err("missing start time declaration")),
            }
            "--no-pbc" => config.pbc = false,
            _ => return Some(Err("invalid input argument")),
        };
        Some(Ok(()))
    }
}

/// Type of analysis to be run.
#[derive(Debug)]
pub enum Task {
    NumFrames,
    Convert,
    Dipole,
    Steinhardt,
    Cluster,
}

/// Available order parameters.
#[derive(Debug)]
#[derive(PartialEq)]
pub enum OrderParameter {
    Theta, Q3, Q4, Q6,
}

impl OrderParameter {
    /// Title line for output files. E.g.
    /// ```text
    /// Steinhardt q3 bond order parameters outputted from md-tools
    /// ```
    pub const fn title_line(&self) -> &'static str {
        match self {
            OrderParameter::Theta => "Water dipole orientations outputted from md-tools",
            OrderParameter::Q3 => "Steinhardt q3 bond order parameters outputted from md-tools",
            OrderParameter::Q4 => "Steinhardt q3 bond order parameters outputted from md-tools",
            OrderParameter::Q6 => "Steinhardt q3 bond order parameters outputted from md-tools",
        }
    }
}