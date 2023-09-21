//! Module to contain small utility functions.

/// Justify atom name by PDB copnventions. Where possible, the element symbol (e.g. `H` (Hydrogen), `NA` (sodium)) should be
/// right-aligned in the first two columns. Examples: ` OW `, ` HW1`, `NA  `.
pub fn justify_atom_name(atom_name: &str) -> String {

    // Array of element symbols for justify_atom_name
    const ELEMENT_SYMBOLS: [&str; 104] = [
        "HE", "LI", "BE", "NE", "NA", "MG", "AL", "SI", "CL", "AR", "CA", "SC", "TI", "CR", "MN", "FE", "CO", "NI",
        "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG",
        "CD", "IN", "SN", "SB", "TE", "XE", "CS", "BA", "HF", "TA", "RE", "OS", "IR", "PT", "AU", "HG", "TI", "PB",
        "BI", "PO", "AT", "RN", "FR", "RA", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG", "CN", "NH", "FL", "MC",
        "LV", "TS", "OG", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU",
        "AC", "TH", "PA", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR"
    ];

    match atom_name.len() {
        1 => format!(" {}  ", atom_name),
        2 => {
            match ELEMENT_SYMBOLS.contains(&atom_name) {
                true => format!("{}  ", atom_name),
                false => format!(" {} ", atom_name),
            }
        }
        3 => {
            match ELEMENT_SYMBOLS.contains(&&atom_name[..2]) {
                true => format!("{} ", atom_name),
                false => format!(" {}", atom_name),
            }
        }
        4 => String::from(atom_name),
        _ => format!("{:>4.4}", atom_name),
    }
}