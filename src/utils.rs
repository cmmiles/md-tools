//! Module to contain small utility functions.

/// Justify atom name by PDB conventions. Where possible, the element symbol (e.g. `O`, `H`, `NA`) should be
/// right-aligned in the first two columns. Examples: ` OW `, ` HW1`, `NA  `.
pub fn justify_atom_name(atom_name: &str) -> String {

    // Array of two character element symbols for justify_atom_name
    // If the first two characters of the atom name don't match any of these then the first character should be the
    // element symbol, and should be justified in the second position (unless atom name is >= 4 chars).
    // Exception: if first character is numeral, it should be justified in the first position.
    const ELEMENT_SYMBOLS: [&str; 104] = [
        "HE", "LI", "BE", "NE", "NA", "MG", "AL", "SI", "CL", "AR", "CA", "SC", "TI", "CR", "MN", "FE", "CO", "NI",
        "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG",
        "CD", "IN", "SN", "SB", "TE", "XE", "CS", "BA", "HF", "TA", "RE", "OS", "IR", "PT", "AU", "HG", "TI", "PB",
        "BI", "PO", "AT", "RN", "FR", "RA", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG", "CN", "NH", "FL", "MC",
        "LV", "TS", "OG", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU",
        "AC", "TH", "PA", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR"
    ];

    // Array of numerals for justify_atom_name
    const NUMERALS: [&str; 10] = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"];

    match atom_name.len() {
        1 => format!(" {}  ", atom_name),
        2 => {
            match ELEMENT_SYMBOLS.contains(&atom_name) {
                true => format!("{}  ", atom_name),
                false => match NUMERALS.contains(&&atom_name[..1]) {
                    true => format!("{}  ", atom_name),
                    false => format!(" {} ", atom_name),
                }
            }
        }
        3 => {
            match ELEMENT_SYMBOLS.contains(&&atom_name[..2]) {
                true => format!("{} ", atom_name),
                false => match NUMERALS.contains(&&atom_name[..1]) {
                    true => format!("{} ", atom_name),
                    false => format!(" {}", atom_name),
                },
            }
        }
        _ => format!("{:>4.4}", atom_name),
    }
}