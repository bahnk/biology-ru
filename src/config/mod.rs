pub mod global;
pub mod commands;
pub mod uaspire;
pub mod uniprot;

use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct Config {
    pub global: global::GlobalConfig,
    pub commands: commands::Commands,
    pub uaspire: uaspire::Config,
    pub uniprot: uniprot::Config,
}

#[cfg(test)]
mod tests {
    use toml;
    use super::Config;

    #[test]
    fn test_import_from_toml() {
        let config: Config = toml::from_str(r#"
        [global]
        log_level = "INFO"
        hive_tables_dir = "/tmp/tables"
        threads = 4

        [uniprot.similar_proteins]
        url = "http://www.uniprot.org/docs/similar.txt"
        species = ["HUMAN", "MOUSE"]
        "#).unwrap();
        println!("{config:?}");
    }
}
