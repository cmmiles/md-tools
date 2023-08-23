use std::{
    fs,
    io::Write,
    sync::{Arc, Mutex},
};

pub fn water_dipole (frame_indices: Vec<u32>, frame_output: Vec<f64>, time: &f32, f: Arc<Mutex<fs::File>>)
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

pub fn steinhardt (frame_indices: Vec<u32>, frame_output: Vec<f64>, time: f32, f: &mut fs::File)
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