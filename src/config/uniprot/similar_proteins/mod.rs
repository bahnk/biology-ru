use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct Config {
    pub url: Option<String>,
    pub species: Option<Vec<String>>,
}

#[cfg(test)]
mod tests {

    use toml;
    use super::Config;

    #[test]
    fn test_import_from_toml() {
        let config: Config = toml::from_str(r#"
            url = "http://www.uniprot.org/docs/similar.txt"
            species = ["HUMAN", "MOUSE"]
        "#).unwrap();
        println!("{config:?}")
    }
}
