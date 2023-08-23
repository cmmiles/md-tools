use super::*;

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
    pub fn coord_shell<'a>(&self, universe: &'a Universe, cutoff: f64, atom_name: &'a str, opt_pbc: &Option<[f32; 3]>) -> Vec<&'a Atom> {
        let u_iter = universe.atoms.iter()
            .filter(|&x| &x.name[..] == atom_name)
            .filter(|&x| 0.0 < self.dsq(x, opt_pbc) && self.dsq(x, opt_pbc) <= cutoff*cutoff);
        u_iter.collect()
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
    pub fn coord_shell<'b>(&'b self, molecules: &'b Vec<Molecule>, cutoff: f64, opt_pbc: &Option<[f32; 3]>)
        -> Vec<&'b Molecule>
    {
        molecules.iter()
            .filter(|&x| 0.0 < self.atoms[0].dsq(x.atoms[0], opt_pbc) && self.atoms[0].dsq(x.atoms[0], opt_pbc) <= cutoff*cutoff)
            .collect()
    }
}