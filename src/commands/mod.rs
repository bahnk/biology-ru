use clap::Subcommand;

pub mod uaspire;
pub mod uniprot;

#[derive(Subcommand, Debug)]
pub enum Commands {
    #[command(subcommand)]
    Uaspire(uaspire::commands::Commands),
    #[command(subcommand)]
    Uniprot(uniprot::commands::Commands),
}
