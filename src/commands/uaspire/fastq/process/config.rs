use std::path::PathBuf;

use super::args::Args;

#[derive(Debug)]
pub struct Config {
    pub read1: PathBuf,
    pub read2: PathBuf,
    pub sample_name: String,
    pub output_dir: PathBuf,
    pub chunk_size: usize,
    pub parquet_size: usize,
    pub threads: Option<usize>,
    pub log_level: Option<String>,
}

impl From<Args> for Config {
    fn from(a: Args) -> Self {
        Self {
            read1: a.read1,
            read2: a.read2,
            sample_name: a.sample_name,
            output_dir: a.output_dir,
            chunk_size: a.chunk_size,
            parquet_size: a.parquet_size,
            threads: a.threads,
            log_level: a.log_level,
        }
    }
}
