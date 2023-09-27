//! Module to handle outputting analysis results to files.

use std::{
    fs,
    io::Write,
};

use ordinal::Ordinal;
use crate::config::OrderParameter;

/// Writes output to CSV file.
pub fn csv(op: &OrderParameter, frame_indices: Vec<u32>, frame_output: Vec<f64>, csv_indices: &Vec<u32>, f: &mut fs::File)
    -> Result<(), &'static str>
{
    // Each line contains max 9n characters, all ASCII:
    // 8 per molecule for the values (or 0 for NaN) and 1 per molecule for the comma (or newline)
    let mut csv_line = String::with_capacity(9*csv_indices.len());

    // frame_indices are just what's computed for the frame, csv_indices are all molecules which may be included
    // i is current frame_indices index
    let mut i = if csv_indices[0] == frame_indices[0] {
        if frame_output[0].is_nan() { csv_line +="\n"; }
        else { csv_line += &format!("\n{0:.1$}", frame_output[0], 7 - op.magnitude()); }
        1
    } else {
        csv_line += "\n";
        0
    };
    for mol_id in csv_indices[1..].iter() {
        if *mol_id == frame_indices[i] {
            if frame_output[i].is_nan() { csv_line += ","; }
            else { csv_line += &format!(",{0:.1$}", frame_output[i], 7 - op.magnitude()); }
            i += 1
        } else {
            csv_line += ",";
        }
    }

    match write!(f, "{csv_line}") {
        Ok(_) => Ok(()),
        Err(_) => return Err("error writing to CSV output file"),
    }
}

/// Writes order parameter output to file.
pub fn order(op: &OrderParameter, frame_indices: Vec<u32>, frame_output: Vec<f64>, time: f32, f: &mut fs::File)
    -> Result<(), &'static str>
{
    let mut frame_string = String::with_capacity(26 + 18*frame_indices.len());
    frame_string += &format!("\nt={:<10}", time);
    for id in frame_indices.iter() { frame_string += &format!(" {id:>8}"); };
    frame_string += &format!("\n            ");
    for theta in frame_output.iter() { frame_string += &format!(" {theta:8.*}", 7 - op.magnitude()); };

    match write!(f, "{frame_string}") {
        Ok(_) => Ok(()),
        Err(_) => return Err("error writing to order parameter output file"),
    }
}

/// Writes q6 clustering output to file.
pub fn q6_clustering(clusters: Vec<Vec<u32>>, time: f32, f: &mut fs::File)
    -> Result<(), &'static str>
{   
    let max_string_len: usize = match clusters.len() {
        1 => 60 + 9*clusters[0].len(),
        _ => 40 + 25*clusters.len() + 9*clusters.iter().map(|cls| cls.len()).sum::<usize>(),
    };
    
    let mut frame_string = String::with_capacity(max_string_len);
    if clusters.len() > 1 {
        frame_string += "\nSize of 1st";
        for i in 1..clusters.len() { frame_string += &format!(", {}", Ordinal(i+1)); }
        frame_string += &format!(" largest clusters at t={}: {}", time, clusters[0].len());
        for cluster in clusters[1..].iter() { frame_string += &format!(", {}", cluster.len()); }
    } else {
        frame_string += &format!("\nSize of largest cluster at t={}: {}", time, clusters[0].len());
    }
    for cluster in clusters.into_iter() {
        frame_string += "\nMolecules: ";
        for id in cluster {
            frame_string += &format!(" {id:>8}");
        }
    }

    match write!(f, "{frame_string}") {
        Ok(_) => Ok(()),
        Err(_) => return Err("error writing to clustering output file"),
    }
}