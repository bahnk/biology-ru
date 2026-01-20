use config::builder::DefaultState;
use config::{Config, ConfigBuilder, Environment, File};
use dotenvy::dotenv;
use std::path::PathBuf;

#[derive(Debug, Default)]
struct PartialConfig {
    hive_tables_dir: String,
}


fn load_settings() {
    //dotenv().ok();

    println!("{}", std::env::var("HIVE_TABLES_DIR").unwrap());

    //let config_file = args.config.to_str().ok_or("Invalid config path")?;
    //let settings = ConfigBuilder::<DefaultState>::default()
    //    .add_source(File::with_name(config_file))
    //    .add_source(Environment::default())
    //    .build()?;

    //Ok(settings)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use std::fs;
    use tempfile::NamedTempFile;
    use std::io::Write;

    #[test]
    fn create_test_file() {
        let mut config_file = NamedTempFile::new().unwrap();
        writeln!(config_file, r#"
        [uniprot.similar]
        url = "http://www.uniprot.org/docs/similar.txt"
        species = ["HUMAN", "MOUSE"]
        "#);

        std::env::set_var("HIVE_TABLES_DIR", "/tmp/hive_tables");

        load_settings();
    }
}
