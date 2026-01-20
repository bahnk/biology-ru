
use super::args::Args;

#[derive(Debug)]
pub struct Config {
    pub url: String,
    pub table_name: String,
}

impl From<Args> for Config {
    fn from(a: Args) -> Self {
        Self {
            url: a.url,
            table_name: a.table_name,
        }
    }
}
