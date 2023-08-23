use std::{
    fmt, fs,
    path::{PathBuf, Path},
    io::Write,
    thread,
    sync::{mpsc, Arc, Mutex},
};
use xdrfile::*;
use closure::closure;
use crate::analysis;
use crate::output;
use crate::config::OrderParameter;
use super::*;

/// Contains the data for the whole system.
/// This is intended to be a mutable structure, with the coordinates mutated for each step in the trajectory.
#[derive(Debug)]
pub struct Universe {
    /// Contains the Atoms for the whole system
    pub atoms: Vec<Atom>,
    /// Dimensions of the simulation box as a 3D array
    pub box_dimensions: [f32; 3],
    /// Location of the trajectory file (optional)
    pub traj: Option<PathBuf>,
    /// Whether to use periodic boundary conditions (default: true)
    pub pbc: bool,
    /// Current time
    pub time: f32,
}

impl fmt::Display for Universe {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let box_dimensions = self.box_dimensions;
        write!(f,
            "Universe{{{} Atoms, box dimensions: [{} {} {}]}}",
            self.atoms.len(), box_dimensions[0], box_dimensions[1], box_dimensions[2]
        )
    }
}

impl Universe {
    /// Constructor method for Universe, reads a structure file (currently only gro files supported).
    /// Optional trajectory file specified (xtc and trr formats supported) which is not read at time of construction.
    pub fn new(sfile: Option<PathBuf>, xtc_file: Option<PathBuf>, pbc: bool) -> Self {
        let (atoms, box_dimensions) = match sfile {
            Some(sfile) => match (&sfile).extension().expect("structure file must have an extension").to_str().unwrap() {
                "gro" => Universe::from_gro(&*sfile),
                ext => panic!("I don't know how to read .{} files", ext),
            }
            None => panic!("no structure file specified"),
        };

        let traj = xtc_file;
        Self { atoms, box_dimensions, traj, pbc, time: 0.0 }
    }

    /// Read a gro file.
    fn from_gro(gro_file: &Path) -> (Vec<Atom>, [f32; 3]) {
        let gro_file = fs::read_to_string(gro_file).unwrap();
        let mut gro_iter = gro_file.lines();

        match gro_iter.next() {
            Some(_) => {},
            None => panic!("empty gro input")
        };

        let natoms: u32 = gro_iter.next().expect("invalid / missing atom count line in gro")
            .trim().parse().expect("invalid / missing atom count line in gro");

        let mut atoms: Vec<Atom> = Vec::new();
        for i in 0..natoms {
            atoms.push(match Atom::from_gro(gro_iter.next().unwrap()) {
                Ok(atom) => atom,
                Err(err) => panic!("error in gro line {}: {}", i+3, err),
            });
        }

        let box_dimensions_str = gro_iter.next().expect("missing box dimensions in gro");
        let mut box_dimensions_iter = box_dimensions_str.split(" ").filter(|&s| s != "").map(|s| s.parse::<f32>()
            .expect("invalid box dimensions in gro"));

        let box_dimensions = [
            box_dimensions_iter.next().expect("invalid box dimensions in gro"),
            box_dimensions_iter.next().expect("invalid box dimensions in gro"),
            box_dimensions_iter.next().expect("invalid box dimensions in gro"),
        ];

        ( atoms, box_dimensions )
    }

    /// Reads the trajectory if it exists (xtc and trr files supported).
    fn read_traj(&self) -> Result<MultiTrajIter, &'static str> {
        match &self.traj {
            Some(tfile) => match (&tfile).extension().expect("trajectory file must have an extension").to_str().unwrap() {
                "xtc" => Self::read_xtc(&*tfile, self.atoms.len()),
                "trr" => Self::read_trr(&*tfile, self.atoms.len()),
                ext => panic!("I don't know how to read .{} files", ext),
            }
            None => Err("no trajectory file specified"),
        }
    }

    /// Reads an xtc file
    fn read_xtc(xtc_file: &Path, natoms: usize) -> Result<MultiTrajIter, &'static str> {
        let mut trj = match XTCTrajectory::open_read(xtc_file) {
            Ok(trj) => trj,
            Err(_) => return Err("error reading xtc file"),
        };
        if trj.get_num_atoms().unwrap() != natoms {
            Err("structure and trajectory files do not match")
        } else {
            Ok(MultiTrajIter::XTC(trj.into_iter()))
        }
    }

    /// Reads a trr file
    fn read_trr(trr_file: &Path, natoms: usize) -> Result<MultiTrajIter, &'static str> {
        let mut trj = match TRRTrajectory::open_read(trr_file) {
            Ok(trj) => trj,
            Err(_) => return Err("error reading trr file"),
        };
        if trj.get_num_atoms().unwrap() != natoms {
            Err("structure and trajectory files do not match")
        } else {
            Ok(MultiTrajIter::TRR(trj.into_iter()))
        }
    }

    /// Writes the current frame to a gro file.
    pub fn write_gro(&self, gro_file: &str) {
        let err_msg = &format!("error writing to file {}", gro_file);
        let mut f = fs::File::create(gro_file).expect(err_msg);
        let box_dimensions = self.box_dimensions;
        write!(f, "Structure outputted from md-tools\n").expect(err_msg);
        write!(f, "{:>5}\n", &self.atoms.len()).expect(err_msg);
        for atom in self.atoms.iter() {
            write!(f, "{}\n", atom.to_gro()).expect(err_msg);
        }
        write!(f, "{:10.5}{:10.5}{:10.5}", box_dimensions[0], box_dimensions[1], box_dimensions[2]).expect(err_msg);
    }

    /// Returns a vector of molecules, assumes molecules are grouped in order in the structure file.
    pub fn get_molecules(&self) -> Vec<Molecule> {
        let mut molecules: Vec<Molecule> = Vec::new();
        for atom in &self.atoms {
            match &mut molecules.last_mut() {
                None => molecules.push(Molecule { id: atom.resid, name: &atom.resname[..], atoms: vec!(&atom) }),
                Some(mol) if mol.id == atom.resid => mol.atoms.push(&atom),
                Some(_) => molecules.push(Molecule { id: atom.resid, name: &atom.resname[..], atoms: vec!(&atom) }),
            }
        }
        molecules
    }

    /// Returns the number of frames in the trajectory.
    pub fn get_nframes(&self) -> Result<usize, &'static str> {
        let traj_iter = self.read_traj()?;
        Ok(traj_iter.count())
    }

    /// Convert the Universe to the desired format.
    pub fn convert(&mut self, mint: &u32, maxt: &u32, tstep: &u32, outfile: &Option<PathBuf>) -> Result<(), &'static str> {
        match outfile {
            Some(outfile) => match (&outfile).extension().expect("output file must have an extension").to_str().unwrap() {
                "gro" => self.write_gro_range(mint, maxt, tstep, &outfile.with_extension("").to_str().unwrap()),
                "xtc" => self.write_xtc(mint, maxt, tstep, &outfile.with_extension("").to_str().unwrap()),
                "trr" => self.write_trr(mint, maxt, tstep, &outfile.with_extension("").to_str().unwrap()),
                ext => panic!("I don't know how to write .{} files", ext),
            }
            None => Err("no output file specified"),
        }
    }

    /// Write a range of gro files, from `mint` to `maxt` with stride `tstep`.
    fn write_gro_range(&mut self, mint:&u32, maxt: &u32, tstep: &u32, outfile_str: &str) -> Result<(), &'static str> {
        let traj_iter = self.read_traj()?.filter_map(|x| x.ok())
            .filter(|x| x.time as u32 >= *mint && x.time as u32 <= *maxt && x.time as u32 % *tstep == 0);

        for frame in traj_iter {
            self.load_frame(&frame);
            self.write_gro(&format!("outputs/{}_{}.gro", outfile_str, frame.time));
        }
        Ok(())
    }

    /// Write an xtc (sub)trajectory, from `mint` to `maxt` with stride `tstep`.
    fn write_xtc(&self, mint: &u32, maxt: &u32, tstep: &u32, outfile_str: &str) -> Result<(), &'static str> {
        let traj_iter = self.read_traj()?.filter_map(|x| x.ok())
            .filter(|x| x.time as u32 >= *mint && x.time as u32 <= *maxt && x.time as u32 % *tstep == 0);
        
        let mut trj = match XTCTrajectory::open_write(format!("outputs/{}.xtc", outfile_str)) {
            Ok(trj) => trj,
            Err(_) => return Err("error opening output file"),
        };
        for frame in traj_iter {
            if let Err(_) = trj.write(&frame) {return Err("error writing to xtc file")}
        }

        Ok(())
    }

    /// Write a trr (sub)trajectory, from `mint` to `maxt` with stride `tstep`.
    fn write_trr(&self, mint: &u32, maxt: &u32, tstep: &u32, outfile_str: &str) -> Result<(), &'static str> {
        let traj_iter = self.read_traj()?.filter_map(|x| x.ok())
            .filter(|x| x.time as u32 >= *mint && x.time as u32 <= *maxt && x.time as u32 % *tstep == 0);
        
        let mut trj = match TRRTrajectory::open_write(format!("outputs/{}.trr", outfile_str)) {
            Ok(trj) => trj,
            Err(_) => return Err("error opening output file"),
        };
        for frame in traj_iter {
            if let Err(_) = trj.write(&frame) {return Err("error writing to trr file")}
        }

        Ok(())
    }

    /// Mutates the Universe to the coordinates in the given Frame.
    fn load_frame(&mut self, frame: &Frame) {
        let mut frame_iter = frame.coords.iter();
        self.time = frame.time;
        for atom in &mut self.atoms {
            atom.coords = *frame_iter.next().unwrap();
        }
        self.box_dimensions = [
            frame.box_vector[0][0],
            frame.box_vector[1][1],
            frame.box_vector[2][2],
        ];
    }

    fn request_pbc(&self) -> Option<[f32;3]> {
        match self.pbc {
            true => Some(self.box_dimensions.clone()),
            false => None
        }
    }

    pub fn first_coord_shell<'a>(&'a self, id: usize, atom_name: &'a str) -> Vec<&'a Atom> {
        let opt_pbc = self.request_pbc();
        self.atoms[id].coord_shell(self, 0.3, atom_name, &opt_pbc)
    }

    /// Calculates the water dipoles for the frames specified in the Config
    pub fn water_dipole(
        &mut self,
        resname: &str,
        ref_axis: Vector,
        mint: &u32,
        maxt: &u32,
        tstep: &u32,
        outfile: &Option<PathBuf>,
    )
        -> Result<(), &'static str>
    {
        let outfile = match outfile {
            Some(outfile) => &outfile.to_str().unwrap(),
            None => "outputs/dipole.out",
        };
        let err_msg = &format!("error creating file {}", outfile);
        let opt_pbc = self.request_pbc();

        let mut f = fs::File::create(outfile).expect(err_msg);
        write!(f, "Water dipole orientations outputted from md-tools").expect(err_msg);
        let f = Arc::new(Mutex::new(f));

        let traj_iter = self.read_traj()?.filter_map(|x| x.ok())
            .filter(|x| x.time as u32 >= *mint && x.time as u32 % *tstep == 0);
        
        for frame in traj_iter {
            if frame.time as u32 > *maxt { break; }
            self.load_frame(&frame);
            // Full list of molecules being analysed -- there may be additional framewise filtering
            let molecules: Vec<Molecule> = self.get_molecules().into_iter()
                .filter(|mol| mol.name == resname && mol.atoms.len() > 2)
                .collect();
            let (frame_indices, frame_output) = analysis::water_dipole(&molecules, &ref_axis, &opt_pbc)?;
            let f_clone = f.clone();
            let time = self.time;
            thread::spawn(move || output::water_dipole(frame_indices, frame_output, &time, f_clone));
        }

        Ok(())
    }

    pub fn steinhardt(
        &mut self, 
        resname: &str,
        op_list: &Vec<OrderParameter>,
        mint: &u32,
        maxt: &u32,
        tstep: &u32,
        outfile: &Option<PathBuf>,
    )
        -> Result<(), &'static str>
    {
        if op_list.len() == 0 { panic!("No Steinhardt parameters selected.") }
        thread::scope(|analysis_scope| -> Result<(), &'static str> {
            let (tx, rx) = mpsc::channel();

            analysis_scope.spawn(move || {
                let outfile = match &outfile {
                    Some(outfile) => outfile.to_str().unwrap(),
                    None => "outputs/steinhardt.out",
                };
                let mut f_list = Vec::new();
                for op in op_list.iter() {
                    let mut f = match op {
                        OrderParameter::Q3 => fs::File::create(String::from(outfile) + ".q3"),
                        OrderParameter::Q4 => fs::File::create(String::from(outfile) + ".q4"),
                        OrderParameter::Q6 => fs::File::create(String::from(outfile) + ".q6"),
                    }.expect("error creating output file");
                    write!(f, "Steinhardt {:?} bond order parameters outputted from md-tools", op)
                        .expect("error writing to output file");
                    f_list.push(f);
                }

                for (op_index, time, frame_indices, frame_output) in rx {
                    output::steinhardt(frame_indices, frame_output, time, &mut f_list[op_index])
                        .expect("error writing to output file");
                }
            });
        
            let traj_iter = self.read_traj()?.filter_map(|x| x.ok())
                .filter(|x| x.time as u32 >= *mint && x.time as u32 % *tstep == 0);

            for frame in traj_iter {
                if frame.time as u32 > *maxt { break; }
                self.load_frame(&frame);
                let time = self.time;
                let opt_pbc = self.request_pbc();
                // Full list of molecules being analysed -- there may be additional framewise filtering
                let molecules: Vec<Molecule> = self.get_molecules().into_iter()
                    .filter(|mol| mol.name == resname && mol.atoms.len() > 2)
                    .collect();
                
                thread::scope(|frame_scope| {
                    for (i, op) in op_list.iter().enumerate() {
                        let tx_clone = tx.clone();
                        frame_scope.spawn(closure!(ref molecules, ref opt_pbc, || {
                            let (frame_indices, frame_output) = match op {
                                OrderParameter::Q3 => analysis::steinhardt(3, &molecules, &opt_pbc),
                                OrderParameter::Q4 => analysis::steinhardt(4, &molecules, &opt_pbc),
                                OrderParameter::Q6 => analysis::steinhardt(6, &molecules, &opt_pbc),
                            };
                            tx_clone.send((i, time, frame_indices, frame_output)).expect("Error sending write data");
                        }));
                    }
                });
            }
            Ok(())
        })?;
        
        Ok(())
    }

    #[cfg(test)]
    pub fn compare_gro(&self, atoms: Vec<Atom>, box_dimensions: [f32; 3]) -> bool {
        self.atoms == atoms && self.box_dimensions == box_dimensions
    }
}