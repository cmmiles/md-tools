use std::{process, env};
use md_tools::config::Config;

fn main() {

    let config = Config::build(env::args()).unwrap_or_else(|err| {
        eprintln!("Issue parsing arguments: {err}.");
        process::exit(1);
    });

    if let Err(e) = md_tools::run(config) {
        eprintln!("Application error: {e}.");
        process::exit(1);
    }

}
