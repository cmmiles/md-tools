use std::{
    fmt, fs,
    path::{PathBuf, Path},
    io::Write,
    rc::Rc,
    thread,
    sync::{mpsc, Arc, Mutex},
};
use xdrfile::*;
use crate::math::*;
use crate::analysis;
use crate::output;

/// An object which has coordinates, e.g. Atom, centre of mass
pub trait Coords {
    fn coords(&self) -> &[f32; 3];

    /// Calculates the vector to another Coords
    fn vector_to_coord(&self, other: &dyn Coords, opt_pbc: &Option<[f32; 3]>) -> Vector {
        let self_coords = self.coords();
        let other_coords = other.coords();
        let mut dx = self_coords[0] - other_coords[0];
        let mut dy = self_coords[1] - other_coords[1];
        let mut dz = self_coords[2] - other_coords[2];
        if let Some(box_dimensions) = opt_pbc {
            let box_dimensions = box_dimensions;
            while dx < -0.5 * box_dimensions[0] { dx += box_dimensions[0]; }
            while dx > 0.5 * box_dimensions[0] { dx -= box_dimensions[0]; }
            while dy < -0.5 * box_dimensions[1] { dy += box_dimensions[1]; }
            while dy > 0.5 * box_dimensions[1] { dy -= box_dimensions[1]; }
            while dz < -0.5 * box_dimensions[2] { dz += box_dimensions[2]; }
            while dz > 0.5 * box_dimensions[2] { dz -= box_dimensions[2]; }
        }
        Vector::new(dx, dy, dz)
    }

    /// Calculates the squared distance to another Coords
    fn dsq(&self, other: &dyn Coords, opt_pbc: &Option<[f32; 3]>) -> f32 {
        self.vector_to_coord(other, opt_pbc).dsq()
    }
}

/// Contains data for a single atom.
#[derive(PartialEq)]
pub struct Atom {
    /// Atom index
    pub id: u32,
    /// Atom name, e.g. OW
    pub name: String,
    /// Residue name, e.g. SOL
    pub resname: String,
    /// Residue index
    pub resid: u32,
    /// Cartesian position (nm)
    pub coords: [f32; 3],
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,
            "Atom #{}{{{}:{}}}",
            self.id, self.resname, self.name
        )
    }
}

impl fmt::Debug for Atom {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,
            "Atom #{}{{{}:{}}}",
            self.id, self.resname, self.name
        )
    }
}

impl Coords for Atom {
    fn coords(&self) -> &[f32; 3] { &self.coords }
}

impl Atom {
    /// Constructor method for Atom, intended for testing purposes only.
    #[cfg(test)]
    pub fn new(id: u32, name: &str, resname: &str, resid: u32, coords: [f32; 3]) -> Self {
        Self {
            id, resid, coords,
            name: String::from(name),
            resname: String::from(resname),
        }
    }

    /// Create an Atom from a line of a gro file.
    pub fn from_gro(line: &str) -> Result<Self, &'static str> {
        if line.len() < 44 { return Err("invalid line") }
        let resid: u32 = match &line[0..5].trim().parse() {
            Ok(val) => *val,
            Err(_) => return Err("invalid residue id"),
        };
        let id: u32 = match &line[15..20].trim().parse() {
            Ok(val) => *val,
            Err(_) => return Err("invalid atom id"),
        };
        let coords: [f32; 3] = [
            match &line[20..28].trim().parse() {
                Ok(val) => *val,
                Err(_) => return Err("invalid x coordinate"),
            },
            match &line[28..36].trim().parse() {
                Ok(val) => *val,
                Err(_) => return Err("invalid y coordinate"),
            },
            match &line[36..44].trim().parse() {
                Ok(val) => *val,
                Err(_) => return Err("invalid z coordinate"),
            }
        ];

        Ok(Self {
            id, resid, coords,
            name: String::from(*&line[10..15].trim()),
            resname: String::from(*&line[5..10].trim()),
        })
    }

    /// Write a line for the Atom in gro format.
    pub fn to_gro(&self) -> String {
        let coords = self.coords;
        format!("{:>5}{:<5.5}{:>5.5}{:>5}{:8.3}{:8.3}{:8.3}",
            self.resid % 100_000, self.resname, self.name, self.id % 100_000,
            coords[0], coords[1], coords[2])
    }

    /// Calculates the coordination shell of the Atom
    pub fn coord_shell<'a>(&self, universe: &'a Universe, cutoff: f32, atom_name: &'a str, opt_pbc: &Option<[f32; 3]>) -> Vec<&'a Atom> {
        let u_iter = universe.atoms.iter()
            .filter(|&x| &x.name[..] == atom_name)
            .filter(|&x| 0.0 < self.dsq(x, opt_pbc) && self.dsq(x, opt_pbc) <= cutoff*cutoff);
        u_iter.collect()
    }
}

/// Generic Point struct for positional data alone, e.g. centre of mass
pub struct Point([f32; 3]);

impl Coords for Point {
    fn coords(&self) -> &[f32; 3] { &self.0 }
}

impl Point {
    /// Finds the centre of a set of coordinates !not taking mass into account!
    /// There is potential for the centre to be incorrect if the points are spaced across more than 50% of the simulation box
    /// and periodic boundary conditions are in use.
    pub fn centre(points: Vec<&dyn Coords>, opt_pbc: &Option<[f32; 3]>) -> Point {
        let n = points.len() as f32;
        let mut coords_iter = points.iter().map(|point| point.coords());
        let [mut cx, mut cy, mut cz] = {
            let coords = coords_iter.next().unwrap();
            [coords[0], coords[1], coords[2]]
        };

        let [min, max] = match opt_pbc {
            Some(box_dimensions) => {
                let box_dimensions = box_dimensions;
                [
                    [
                        cx - 0.5 * box_dimensions[0],
                        cy - 0.5 * box_dimensions[1],
                        cz - 0.5 * box_dimensions[2],
                    ],
                    [
                        cx + 0.5 * box_dimensions[0],
                        cy + 0.5 * box_dimensions[1],
                        cz + 0.5 * box_dimensions[2]
                    ],
                ]
            }
            None => [
                [f32::MIN, f32::MIN, f32::MIN],
                [f32::MAX, f32::MAX, f32::MAX],
            ],
        };

        for coords in coords_iter {
            let [mut x, mut y, mut z] = [coords[0], coords[1], coords[2]];
            if let Some(box_dimensions) = opt_pbc {
                let box_dimensions = box_dimensions;
                while x < min[0] { x += box_dimensions[0]; }
                while x > max[0] { x -= box_dimensions[0]; }
                while y < min[1] { y += box_dimensions[1]; }
                while y > max[1] { y -= box_dimensions[1]; }
                while z < min[2] { z += box_dimensions[2]; }
                while z > max[2] { z -= box_dimensions[2]; }
            }
            cx += x;
            cy += y;
            cz += z;
        }
        Point ([cx/n, cy/n, cz/n])
    }
}

#[derive(Debug)]
/// Contains references to atoms for a single molcule, along with molecule details
pub struct Molecule<'a> {
    /// Molecule index
    pub id: u32,
    /// Molecule name
    pub name: &'a str,
    /// Atoms within molecule
    pub atoms: Vec<&'a Atom>,
}

impl<'a> Molecule<'a> {
    /// Calculates the coordination shell of the Molecule, using the first atom as the centre
    pub fn coord_shell<'b>(&'b self, molecules: &'b Vec<Molecule>, cutoff: f32, opt_pbc: &Option<[f32; 3]>)
        -> Vec<&'b Molecule>
    {
        molecules.iter()
            .filter(|&x| 0.0 < self.atoms[0].dsq(x.atoms[0], opt_pbc) && self.atoms[0].dsq(x.atoms[0], opt_pbc) <= cutoff*cutoff)
            .collect()
    }
}

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
        mint: &u32,
        maxt: &u32,
        tstep: &u32,
        outfile: &Option<PathBuf>,
    )
        -> Result<(), &'static str>
    {
        thread::scope(|s| -> Result<(), &'static str> {
            let (tx, rx) = mpsc::channel();

            s.spawn(move || {
                let outfile = match &outfile {
                    Some(outfile) => outfile.to_str().unwrap(),
                    None => "outputs/steinhardt.out",
                };
                let err_msg = &format!("error writing to file {}", outfile);
                let mut f = fs::File::create(outfile).expect(err_msg);
                write!(f, "Steinhardt bond order parameters outputted from md-tools").expect(err_msg);

                for (time, frame_indices, frame_output) in rx {
                    output::steinhardt(frame_indices, frame_output, time, &mut f).expect(err_msg);
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
                let (frame_indices, frame_output) = analysis::steinhardt(3, &molecules, &opt_pbc)?;
                tx.send((time, frame_indices, frame_output)).expect("Error sending write data");
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