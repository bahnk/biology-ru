use clap::Parser;

#[derive(Debug, Parser)]
pub struct Args {
    // Input FASTQ files
    #[arg()]
    pub read1: std::path::PathBuf,
    #[arg()]
    pub read2: std::path::PathBuf,
    
    // Sample name
    #[arg(long, short)]
    pub sample_name: String,
    
    // Output directory
    #[arg(long, short, default_value = "./output")]
    pub output_dir: std::path::PathBuf,
    
    // Chunk and parquet sizes
    #[arg(long, short, default_value = "10000")]
    pub chunk_size: usize,
    #[arg(long, short, default_value = "10000")]
    pub parquet_size: usize,

    // Optional threads and log level
    #[arg(long, short)]
    pub threads: Option<usize>,
    #[arg(long, short)]
    pub log_level: Option<String>,
}

#[cfg(test)]
mod tests {
    use super::Args;
    use clap::Parser;
    use std::path::PathBuf;

    #[test]
    fn parses_minimal_args() {
        let args = Args::try_parse_from([
            "prog",
            "r1.fastq",
            "r2.fastq",
            "--sample-name",
            "sample1",
        ])
        .unwrap();
    
        assert_eq!(args.read1, PathBuf::from("r1.fastq"));
        assert_eq!(args.read2, PathBuf::from("r2.fastq"));
        assert_eq!(args.sample_name, "sample1");
        assert_eq!(args.output_dir, PathBuf::from("./output"));
        assert_eq!(args.chunk_size, 10_000);
        assert_eq!(args.parquet_size, 10_000);
    }
}
