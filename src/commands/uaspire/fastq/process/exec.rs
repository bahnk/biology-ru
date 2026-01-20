use rayon::ThreadPoolBuilder;
use tracing_subscriber;

use crate::config::global::GlobalConfig;
use crate::uaspire::fastq::process_fastq;

use super::config::Config;

pub fn execute(cfg: Config, global: GlobalConfig) {

    let level_str = cfg
        .log_level
        .unwrap_or_else(|| global.log_level.clone());

    let level = level_str
        .parse::<tracing::Level>()
        .unwrap_or(tracing::Level::INFO);

    tracing_subscriber::fmt()
        .compact()
        .with_max_level(level)
        .init();

    let threads = cfg.threads.unwrap_or(global.threads);

    ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .expect("Failed to build thread pool");

    process_fastq(
        &cfg.read1.to_string_lossy(),
        &cfg.read2.to_string_lossy(),
        &cfg.sample_name,
        &cfg.output_dir.to_string_lossy(),
        cfg.chunk_size,
        cfg.parquet_size,
    );
}

#[cfg(test)]
mod tests {
    use tempfile::TempDir;
    use std::path::PathBuf;
    use super::execute;
    use super::super::config::Config;
    use crate::config::global::GlobalConfig;

    #[test]
    fn process_small_fastq() {
        let cfg = Config {
            read1: PathBuf::from("test/data/fastq/uaspire/example_R1.fastq.gz"),
            read2: PathBuf::from("test/data/fastq/uaspire/example_R2.fastq.gz"),
            sample_name: "TestSample".to_string(),
            output_dir: TempDir::new().unwrap().keep(),
            chunk_size: 10_000,
            parquet_size: 10_000,
            threads: None,
            log_level: None,
        };

        let global_cfg = GlobalConfig {
            threads: 1,
            log_level: "INFO".to_string(),
            hive_tables_dir: TempDir::new().unwrap().keep(),
        };

        execute(cfg, global_cfg);
    }
}
