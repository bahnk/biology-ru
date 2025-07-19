use clap::{Parser, Subcommand};
use rayon::ThreadPoolBuilder;
use tracing_subscriber;

use crate::uaspire::fastq::process_fastq;

#[derive(Subcommand, Debug, Clone)]
pub enum Commands {
    #[command(name = "process-sample")]
    ParseFastq(ParseFastqCommand),
}

#[derive(Parser, Debug, Clone)]
pub struct ParseFastqCommand {
    // Input FASTQ files
    #[arg()]
    read1: std::path::PathBuf,
    #[arg()]
    read2: std::path::PathBuf,

    // Sample name
    #[arg(long, short)]
    sample_name: String,

    // Output directory
    #[arg(long, short, default_value = "./output")]
    output_dir: std::path::PathBuf,

    // Chunk and parquet sizes
    #[arg(long, short, default_value = "10000")]
    chunk_size: usize,
    #[arg(long, short, default_value = "10000")]
    parquet_size: usize,
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
                &cmd.sample_name,
                &cmd.output_dir.to_string_lossy(),
                cmd.chunk_size,
                cmd.parquet_size,
            );
        }
    }
}
