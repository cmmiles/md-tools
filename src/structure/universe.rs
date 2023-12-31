use std::{
    fmt, fs,
    path::{PathBuf, Path},
    io::Write,
    thread,
    sync::mpsc,
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
    /// Constructor method for Universe, reads a structure file (currently only GRO and PDB files supported).
    /// Optional trajectory file specified (xtc and trr formats supported) which is not read at time of construction.
    pub fn new(sfile: Option<PathBuf>, xtc_file: Option<PathBuf>, pbc: bool) -> Self {
        let (atoms, box_dimensions) = match sfile {
            Some(sfile) => match (&sfile).extension().expect("structure file must have an extension").to_str().unwrap() {
                "gro" => Self::from_gro(&*sfile),
                "pdb" => Self::from_pdb(&*sfile),
                ext => panic!("I don't know how to read .{} files", ext),
            }
            None => panic!("no structure file specified"),
        };

        let traj = xtc_file;
        Self { atoms, box_dimensions, traj, pbc, time: 0.0 }
    }

    /// Read a GRO file.
    fn from_gro(gro_file: &Path) -> (Vec<Atom>, [f32; 3]) {
        let gro_file = fs::read_to_string(gro_file).unwrap();
        let mut gro_iter = gro_file.lines();

        gro_iter.next().expect("empty GRO input");

        let natoms: u32 = gro_iter.next().expect("invalid / missing atom count line in GRO")
            .trim().parse().expect("invalid / missing atom count line in gro");

        let mut atoms = Vec::new();
        for i in 0..natoms {
            atoms.push(match Atom::from_gro(gro_iter.next().unwrap()) {
                Ok(atom) => atom,
                Err(err) => panic!("error in GRO line {}: {}", i+3, err),
            });
        }

        let box_dimensions_str = gro_iter.next().expect("missing box dimensions in GRO");
        let mut box_dimensions_iter = box_dimensions_str.split(" ").filter(|&s| s != "").map(|s| s.parse::<f32>()
            .expect("invalid box dimensions in GRO"));

        let box_dimensions = [
            box_dimensions_iter.next().expect("invalid box dimensions in GRO"),
            box_dimensions_iter.next().expect("invalid box dimensions in GRO"),
            box_dimensions_iter.next().expect("invalid box dimensions in GRO"),
        ];

        ( atoms, box_dimensions )
    }

    /// Read a PDB file, ignores all lines that are not CRYST1, ATOM or END.
    fn from_pdb(pdb_file: &Path) -> (Vec<Atom>, [f32; 3]) {
        let pdb_file = fs::read_to_string(pdb_file).unwrap();
        let pdb_iter = pdb_file.lines();
        let mut box_dimensions = None;
        let mut atoms = Vec::new();
        for line in pdb_iter {
            if line.len() >= 3 && &line[..3] == "END" {
                break;
            } else if line.len() >= 6 && &line[..6] == "CRYST1" {
                if let Some(_) = box_dimensions { panic!("multiple CRYST1 records found"); }
                if line.len() < 70 { panic!("PDB: invalid CRYST1 record"); }
                if &line[33..70] != "  90.00  90.00  90.00 P 1           1" {
                    panic!("CRYST1: angles must be 90 degrees, space group P 1, Z=1");
                }

                box_dimensions = Some([
                    line[6..15].replace(" ", "").parse::<f32>().expect("invalid box dimensions in PDB") / 10.0,
                    line[15..24].replace(" ", "").parse::<f32>().expect("invalid box dimensions in PDB") / 10.0,
                    line[24..33].replace(" ", "").parse::<f32>().expect("invalid box dimensions in PDB") / 10.0,
                ]);
            } else if line.len() >= 4 && &line[..4] == "ATOM" {
                atoms.push(match Atom::from_pdb(line) {
                    Ok(atom) => atom,
                    Err(err) => panic!("error reading PDB: {}", err),
                });
            }
        }
        if atoms.len() == 0 { panic! ("empty PDB input"); }
        if let None = box_dimensions { panic! ("PDB: missing CRYST1 record"); }
        (atoms, box_dimensions.unwrap())
    }

    /// Read the trajectory if it exists (xtc and trr files supported).
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

    /// Read an xtc file.
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

    /// Read a trr file.
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

    /// Write the current frame to a GRO file.
    pub fn write_gro(&self, gro_file: &str) {
        let err_msg = &format!("error writing to file {}", gro_file);
        let mut f = fs::File::create(gro_file).expect(err_msg);
        writeln!(f, "Structure outputted from md-tools").expect(err_msg);
        writeln!(f, "{:>5}", &self.atoms.len()).expect(err_msg);
        for atom in self.atoms.iter() { writeln!(f, "{}", atom.to_gro()).expect(err_msg); }
        write!(f,
            "{:10.5}{:10.5}{:10.5}",
            self.box_dimensions[0], self.box_dimensions[1], self.box_dimensions[2]
        ).expect(err_msg);
    }

    /// Write the current frame to a PDB file.
    pub fn write_pdb(&self, pdb_file: &str) {
        let err_msg = &format!("error writing to file {}", pdb_file);
        let mut f = fs::File::create(pdb_file).expect(err_msg);
        writeln!(f, "Structure outputted from md-tools").expect(err_msg);
        writeln!(f,
            "CRYST1{:>9.3}{:>9.3}{:>9.3}  90.00  90.00  90.00 P 1           1",
            self.box_dimensions[0] * 10.0, self.box_dimensions[1] * 10.0, self.box_dimensions[2] * 10.0
        ).expect(err_msg);
        for atom in self.atoms.iter() { writeln!(f, "{}", atom.to_pdb()).expect(err_msg); }
        write!(f, "END").expect(err_msg);
    }

    /// Return a vector of molecules, assumes molecules are grouped in order in the structure file.
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

    /// Return the number of frames in the trajectory.
    pub fn get_nframes(&self) -> Result<usize, &'static str> {
        let traj_iter = self.read_traj()?;
        Ok(traj_iter.count())
    }

    /// Convert the Universe to the desired format.
    pub fn convert(&mut self, mint: &u32, maxt: &u32, tstep: &u32, outfile: &Option<PathBuf>) -> Result<(), &'static str> {
        match outfile {
            Some(outfile) => match (&outfile).extension().expect("output file must have an extension").to_str().unwrap() {
                "gro" => self.convert_sfile(mint, maxt, tstep, &outfile.with_extension("").to_str().unwrap(), "gro"),
                "pdb" => self.convert_sfile(mint, maxt, tstep, &outfile.with_extension("").to_str().unwrap(), "pdb"),
                "xtc" => self.write_xtc(mint, maxt, tstep, &outfile.with_extension("").to_str().unwrap()),
                "trr" => self.write_trr(mint, maxt, tstep, &outfile.with_extension("").to_str().unwrap()),
                ext => panic!("I don't know how to write .{} files", ext),
            }
            None => Err("no output file specified"),
        }
    }

    /// Write a range of GRO files, from `mint` to `maxt` with stride `tstep`.
    /// If no trajectory file is supplied, write a single GRO file.
    fn convert_sfile(&mut self, mint:&u32, maxt: &u32, tstep: &u32, outfile_str: &str, format: &str)
    -> Result<(), &'static str> {
        match self.traj {
            Some(_) => {
                let traj_iter = self.read_traj()?.filter_map(|x| x.ok())
                    .filter(|x| x.time as u32 >= *mint && x.time as u32 <= *maxt && x.time as u32 % *tstep == 0);

                for frame in traj_iter {
                    self.load_frame(&frame);
                    match format {
                        "gro" => self.write_gro(&format!("{}_{}.gro", outfile_str, frame.time)),
                        "pdb" => self.write_pdb(&format!("{}_{}.pdb", outfile_str, frame.time)),
                        ext => panic!("I don't know how to write .{} files", ext),
                    }
                }
            }
            None => {
                match  format {
                    "gro" => self.write_gro(&format!("{}.gro", outfile_str)),
                    "pdb" => self.write_pdb(&format!("{}.pdb", outfile_str)),
                    ext => panic!("I don't know how to write .{} files", ext),
                }
            }
        }
        
        Ok(())
    }

    /// Write an xtc (sub)trajectory, from `mint` to `maxt` with stride `tstep`.
    fn write_xtc(&self, mint: &u32, maxt: &u32, tstep: &u32, outfile_str: &str) -> Result<(), &'static str> {
        let traj_iter = self.read_traj()?.filter_map(|x| x.ok())
            .filter(|x| x.time as u32 >= *mint && x.time as u32 <= *maxt && x.time as u32 % *tstep == 0);
        
        let mut trj = match XTCTrajectory::open_write(format!("{}.xtc", outfile_str)) {
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
        
        let mut trj = match TRRTrajectory::open_write(format!("{}.trr", outfile_str)) {
            Ok(trj) => trj,
            Err(_) => return Err("error opening output file"),
        };
        for frame in traj_iter {
            if let Err(_) = trj.write(&frame) {return Err("error writing to trr file")}
        }

        Ok(())
    }

    /// Mutate the Universe to the coordinates in the given Frame.
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

    /// Returns None if pbcs are off, returns box dimensions if pbcs are on.
    fn request_pbc(&self) -> Option<[f32;3]> {
        match self.pbc {
            true => Some(self.box_dimensions.clone()),
            false => None
        }
    }

    /// Calculate order parameters for the frames specified in the Config.
    pub fn order(
        &mut self, 
        resname: &str,
        op_list: &Vec<OrderParameter>,
        cutoff: &f64,
        ref_axis: Vector,
        mint: &u32,
        maxt: &u32,
        tstep: &u32,
        outfile: &Option<PathBuf>,
    )
        -> Result<(), &'static str>
    {
        if op_list.len() == 0 { panic!("No Steinhardt parameters selected.") }
        let molecules: Vec<Molecule> = self.get_molecules().into_iter()
            .filter(|mol| mol.name == resname && mol.atoms.len() > 2)
            .collect();
        if molecules.len() == 0 { panic!("No molecules in selection!"); }

        // Read the output location if it exists, set it to md-tools.out if not.
        // Check whether to output a regular output file or a csv
        let (outfile, ext) = match &outfile {
            Some(outfile) => (
                outfile.parent().unwrap().join(outfile.file_stem().unwrap()),
                match outfile.extension() {
                    Some(ext) => String::from(".") + ext.to_str().unwrap(),
                    None => String::new(),
                }
            ),
            None => (PathBuf::from("md-tools"), String::from(".out")),
        };

        let csv_flag = ext == ".csv";

        let (csv_header, csv_indices) = if csv_flag {
            let mut csv_header = format!("{}:{}", molecules[0].name, molecules[0].id);
            for molecule in molecules[1..].iter() { csv_header += &format!(",{}:{}", molecule.name, molecule.id); }
            (csv_header, molecules.iter().map(|mol| mol.id).collect())
        } else {
            (String::new(), Vec::new())
        };

        thread::scope(|analysis_scope| -> Result<(), &'static str> {
            let (tx, rx) = mpsc::channel();

            analysis_scope.spawn(move || {
                let mut f_list = Vec::new();
                for op in op_list.iter() {
                    let mut f = match op {
                        OrderParameter::Q3 => fs::File::create(String::from(outfile.to_str().unwrap()) + ".q3" + &ext),
                        OrderParameter::Q4 => fs::File::create(String::from(outfile.to_str().unwrap()) + ".q4" + &ext),
                        OrderParameter::Q6 => fs::File::create(String::from(outfile.to_str().unwrap()) + ".q6" + &ext),
                        OrderParameter::Theta => fs::File::create(String::from(outfile.to_str().unwrap()) + ".theta" + &ext),
                    }.expect("error creating output file");

                    if csv_flag { write!(f, "{}", csv_header).expect("error writing to output file"); }
                    else { write!(f, "{}", op.title_line()).expect("error writing to output file"); }
                    f_list.push(f);
                }

                for (op_index, time, frame_indices, frame_output) in rx {
                    if csv_flag {
                        output::csv(&op_list[op_index], frame_indices, frame_output, &csv_indices, &mut f_list[op_index])
                            .expect("error writing to csv file");
                    } else {
                        output::order(&op_list[op_index], frame_indices, frame_output, time, &mut f_list[op_index])
                            .expect("error writing to output file");
                    }
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
                        frame_scope.spawn(closure!(ref molecules, ref opt_pbc, ref ref_axis, || {
                            let (frame_indices, frame_output) = match op {
                                OrderParameter::Q3 => analysis::steinhardt(3, cutoff, &molecules, &opt_pbc, None),
                                OrderParameter::Q4 => analysis::steinhardt(4, cutoff, &molecules, &opt_pbc, None),
                                OrderParameter::Q6 => analysis::steinhardt(6, cutoff, &molecules, &opt_pbc, None),
                                OrderParameter::Theta => analysis::water_dipole(&molecules, &ref_axis, &opt_pbc),
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

    /// Largest icy cluster analysis protocol (work in progress, currently very hard-coded).
    pub fn q6_clustering(
        &mut self, 
        resname: &str,
        cutoff: &f64,
        mint: &u32,
        maxt: &u32,
        tstep: &u32,
        n_clusters: &usize,
        outfile: &Option<PathBuf>,
    )
        -> Result<(), &'static str>
    {
        const Q6_CUT: f64 = 0.5; // q6 cutoff for including molecules in clustering

        let outfile = match outfile {
            Some(outfile) => &outfile.to_str().unwrap(),
            None => "outputs/md-tools.out.cls",
        };

        let traj_iter = self.read_traj()?.filter_map(|x| x.ok())
            .filter(|x| x.time as u32 >= *mint && x.time as u32 % *tstep == 0);
        
        let mut f = fs::File::create(outfile).expect("error creating  output file");
        write!(f, "Q6 clusters outputted from md-tools").expect("error writing to output file");

        for frame in traj_iter {
            if frame.time as u32 > *maxt { break; }
            self.load_frame(&frame);
            let time = self.time;
            let opt_pbc = self.request_pbc();
            // Full list of molecules being analysed -- there may be additional framewise filtering
            let molecules: Vec<Molecule> = self.get_molecules().into_iter()
                .filter(|mol| mol.name == resname && mol.atoms.len() > 2)
                .collect();

            // Compute number of hydrogen bonds, generate list of indices of molecules with 4
            let n_hbonds = analysis::get_n_hbonds(&molecules, 0.32, 150.0, &opt_pbc);
            let ndx: Vec<usize> = n_hbonds.iter().enumerate()
                .filter_map(|(i, n)| match n {
                    4 => Some(i),
                    _ => None,
                }).collect();

            // Create list of molecules not included in opt_ndx
            let mut extra_molecules: Vec<Molecule> = n_hbonds.iter().enumerate()
                .filter_map(|(i, n)| match n {
                    4 => None,
                    _ => Some(molecules[i].clone()),
                }).collect();

            // Compute Steinhardt local q6 parameter
            let (_, frame_output) = analysis::local_steinhardt_2(6, cutoff, &molecules, &opt_pbc, Some(&ndx));

            // Generate list of molecules to those with q6 >= Q6_CUT and four neighbours and add those with q6 < Q6_CUT
            // to extra_molecules
            let mut filtered_molecules: Vec<Molecule> = Vec::with_capacity(frame_output.len());
            extra_molecules.reserve(frame_output.len());
            for (i, ndx) in ndx.into_iter().enumerate() {
                if frame_output[i] >= Q6_CUT { filtered_molecules.push(molecules[ndx].clone()); }
                else { extra_molecules.push(molecules[ndx].clone()); }
            }

            let clusters: Vec<Vec<u32>> = if filtered_molecules.len() > 0 {
                // Build the adjacency matrix for filtered molecules, run DFS to get the clusters and truncate to n+5 largest
                // Includes a few extra clusters because surface molecules may affect the order of sizes
                let adj_matrix = AdjacencyMatrix::build_from_coords(&filtered_molecules, cutoff, &opt_pbc);
                let mut clusters = adj_matrix.cluster();
                clusters.truncate(*n_clusters+5);

                // Replace indices (from within filtered_molecules) with corresponding Molecule structs
                let mut clusters: Vec<Vec<Molecule>> = clusters.into_iter().map( |cls|
                    cls.into_iter().map(|i| filtered_molecules[i].clone()).collect()
                ).collect();

                // Add surface molecules to clusters
                for cluster in clusters.iter_mut() {
                    let mut surface_mols: Vec<Molecule> = Vec::new();
                    for extra_mol in extra_molecules.iter() {
                        for cluster_mol in cluster.iter() {
                            if cluster_mol.dsq(extra_mol, &opt_pbc) <= cutoff*cutoff {
                                surface_mols.push(extra_mol.clone());
                                break;
                            }
                        }
                    cluster.append(&mut surface_mols);
                    }
                }

                // Reorder (in case surface atoms changed the order of sizes) and truncate to only n largest
                clusters.sort_unstable_by(|a, b| b.len().cmp(&a.len()));
                clusters.truncate(*n_clusters);

                // Replace molecules with resids for output
                clusters.into_iter().map(|cls| cls.into_iter().map(|mol| mol.id).collect()).collect()
            } else {
                (0..*n_clusters).into_iter().map(|_| vec![]).collect()
            };

            output::q6_clustering(clusters, time, &mut f)?;
        }

        Ok(())
    }

    /// Testing function to compare two Universe snapshots.
    #[cfg(test)]
    pub fn compare_gro(&self, atoms: Vec<Atom>, box_dimensions: [f32; 3]) -> bool {
        self.atoms == atoms && self.box_dimensions == box_dimensions
    }
}