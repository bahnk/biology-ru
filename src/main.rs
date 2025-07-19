use clap::Parser;

use biology_ru::cli::{Cli, Commands};
use biology_ru::commands;

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Uniprot(cmd) => {
            commands::uniprot::command(cmd);
        }
        Commands::Uaspire(cmd) => {
            commands::uaspire::command(cmd);
        }
    }
}
