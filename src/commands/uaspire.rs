
use clap::{Parser, Subcommand};
use tracing_subscriber;
use rayon::ThreadPoolBuilder;

use crate::uaspire::fastq::process_fastq;

#[derive(Subcommand, Debug, Clone)]
pub enum Commands {
    #[command(name = "parse-fastq")]
    ParseFastq(ParseFastqCommand),
}

#[derive(Parser, Debug, Clone)]
pub struct ParseFastqCommand {
    #[arg()]
    read1: std::path::PathBuf,
    #[arg()]
    read2: std::path::PathBuf,
    #[arg(long, short, default_value = "./output")]
    output_dir: std::path::PathBuf,
}


pub fn command(cmds: Commands) {
    match cmds {
        Commands::ParseFastq(cmd) => {

            tracing_subscriber::fmt()
                .compact()
                .with_max_level(tracing::Level::INFO)
                .init();

            ThreadPoolBuilder::new()
                .num_threads(10)
                .build_global()
                .expect("Failed to build thread pool");

            process_fastq(
                &cmd.read1.to_string_lossy(),
                &cmd.read2.to_string_lossy(),
                1_000_000,
                &cmd.output_dir.to_string_lossy(),
            );
        }
    }
}
