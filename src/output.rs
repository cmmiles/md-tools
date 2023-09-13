//! Module to handle outputting analysis results to files.

use std::{
    fs,
    io::Write,
    sync::{Arc, Mutex},
};

use ordinal::Ordinal;

/// Writes water dipole output to file. N.B. to be combined with `steinhardt()`.
pub fn water_dipole(frame_indices: Vec<u32>, frame_output: Vec<f64>, time: &f32, f: Arc<Mutex<fs::File>>)
    -> Result<(), &'static str>
{
    let mut frame_string = String::with_capacity(26 + 22*frame_indices.len());
    frame_string.push_str(&format!("\nt={:<10}", time));
    for id in frame_indices.iter() { frame_string.push_str(&format!(" {id:>10}")); };
    frame_string.push_str(&format!("\n            "));
    for theta in frame_output.iter() { frame_string.push_str(&format!(" {theta:10.4}")); };

    match write!(f.lock().unwrap(), "{frame_string}") {
        Ok(_) => Ok(()),
        Err(_) => return Err("error writing to dipole output file"),
    }
}

/// Writes Steinhardt parameter output to file. N.B. to be combined with `water_dipole()`.
pub fn steinhardt(frame_indices: Vec<u32>, frame_output: Vec<f64>, time: f32, f: &mut fs::File)
    -> Result<(), &'static str>
{
    let mut frame_string = String::with_capacity(26 + 18*frame_indices.len());
    frame_string.push_str(&format!("\nt={:<10}", time));
    for id in frame_indices.iter() { frame_string.push_str(&format!(" {id:>8}")); };
    frame_string.push_str(&format!("\n            "));
    for theta in frame_output.iter() { frame_string.push_str(&format!(" {theta:8.6}")); };

    match write!(f, "{frame_string}") {
        Ok(_) => Ok(()),
        Err(_) => return Err("error writing to Steinhardt output file"),
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
        frame_string.push_str("\nSize of 1st");
        for i in 1..clusters.len() { frame_string.push_str(&format!(", {}", Ordinal(i+1))); }
        frame_string.push_str(&format!(" largest clusters at t={}: {}", time, clusters[0].len()));
        for cluster in clusters[1..].iter() { frame_string.push_str(&format!(", {}", cluster.len())); }
    } else {
        frame_string.push_str(&format!("\nSize of largest cluster at t={}: {}", time, clusters[0].len()));
    }
    for cluster in clusters.into_iter() {
        frame_string.push_str("\nMolecules: ");
        for id in cluster {
            frame_string.push_str(&format!(" {id:>8}"));
        }
    }

    match write!(f, "{frame_string}") {
        Ok(_) => Ok(()),
        Err(_) => return Err("error writing to clustering output file"),
    }
}