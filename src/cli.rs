use clap::{Parser, Subcommand};

use crate::commands;

#[derive(Subcommand, Debug)]
pub enum Commands {
    #[command(subcommand)]
    Uniprot(commands::uniprot::Commands),
    #[command(subcommand)]
    Uaspire(commands::uaspire::Commands),
}

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
