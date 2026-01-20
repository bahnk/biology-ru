use clap::Subcommand;

use super::fastq;

#[derive(Debug, Subcommand)]
pub enum Commands {
    #[command(subcommand)]
    Fastq(fastq::commands::Commands),
}
