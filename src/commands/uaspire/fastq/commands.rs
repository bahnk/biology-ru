use clap::Subcommand;

use super::process;

#[derive(Debug, Subcommand)]
pub enum Commands {
    #[command(name = "process")]
    Process(process::args::Args),
}
