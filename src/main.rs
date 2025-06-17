use clap::Parser;

use biology_ru::commands;
use biology_ru::cli::{Cli, Commands};

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Uniprot(cmd) => { commands::uniprot::command(cmd); },
        Commands::Uaspire(cmd) => { commands::uaspire::command(cmd); },
    }
}
