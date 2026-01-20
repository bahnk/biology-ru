pub mod process;

use serde::Deserialize;
use std::path::PathBuf;

#[derive(Debug, Deserialize)]
pub struct Config {
    pub log_level: Option<String>,
    pub hive_tables_dir: Option<PathBuf>,
    pub process: Option<process::Config>,
}
