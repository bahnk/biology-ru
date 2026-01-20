pub mod fastq;

use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct Config {
    pub fastq: Option<fastq::Config>,
}
