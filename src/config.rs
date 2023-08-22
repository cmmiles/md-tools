use std::path::PathBuf;

#[derive(Debug)]
pub struct Config {
    pub task: Task,
    pub sfile: Option<PathBuf>,
    pub tfile: Option<PathBuf>,
    pub outfile: Option<PathBuf>,
    pub op_list: Vec<OrderParameter>,
    pub mint: u32,
    pub maxt: u32,
    pub stride: u32,
    pub pbc: bool,
}

impl Config {
    pub fn build(mut args: impl Iterator<Item = String>) -> Result<Self, &'static str> {
        args.next();

        let task = match args.next() {
            Some(task) => match &task[..] {
                "nframes" => Task::NumFrames,
                "convert" => Task::Convert,
                "dipole" => Task::Dipole,
                "steinhardt" => Task::Steinhardt,
                _ => return Err("invalid task"),
            }
            None => return Err("no arguments passed"),
        };

        let mut config = Self { task,
            sfile: None,
            tfile: None,
            outfile: None,
            op_list: Vec::new(),
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
            "--q3" => if !config.op_list.contains(&OrderParameter::Q3) { config.op_list.push(OrderParameter::Q3) },
            "--q4" => if !config.op_list.contains(&OrderParameter::Q4) { config.op_list.push(OrderParameter::Q4) },
            "--q6" => if !config.op_list.contains(&OrderParameter::Q6) { config.op_list.push(OrderParameter::Q6) },
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

#[derive(Debug)]
pub enum Task {
    NumFrames,
    Convert,
    Dipole,
    Steinhardt,
}

#[derive(Debug)]
#[derive(PartialEq)]
pub enum OrderParameter {
    Q3, Q4, Q6,
}