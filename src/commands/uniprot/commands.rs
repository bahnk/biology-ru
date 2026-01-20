use clap::Subcommand;

use super::similar_proteins;

#[derive(Debug, Subcommand)]
pub enum Commands {
    #[command(name = "similar-proteins")]
    SimilarProtein(similar_proteins::args::Args),
}
