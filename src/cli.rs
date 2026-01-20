use clap::Parser;

use crate::commands::Commands;

#[derive(Parser)]
#[command(
    name = "biology-ru",
    color = clap::ColorChoice::Always,
    author = "Biology CLI",
    version = "1.0.0",
    about = "Biology CLI for various biological data processing tasks",
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}


pub fn run_cli() {
    let _ = Cli::parse();
}
