use std::path::PathBuf;
use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct GlobalConfig {
    pub log_level: String,
    pub hive_tables_dir: PathBuf,
    pub threads: usize,
}
