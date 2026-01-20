use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(name = "similar_proteins")]
pub struct Args {
    #[arg(short, long, default_value = "assets/config")]
    pub config: PathBuf,
    #[arg(
        short,
        long,
        default_value = "http://www.uniprot.org/docs/similar.txt"
    )]
    pub url: String,
    #[arg(short, long, default_value = "similar_proteins")]
    pub table_name: String,
}
