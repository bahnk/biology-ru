pub mod similar_proteins;

use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct Config {
    pub similar_proteins: Option<similar_proteins::Config>,
}

#[cfg(test)]
mod tests {
    use toml;
    use super::Config;

    #[test]
    fn test_import_from_toml() {
        let config: Config = toml::from_str(r#"
            [similar_proteins]
            url = "http://www.uniprot.org/docs/similar.txt"
            species = ["HUMAN", "MOUSE"]
        "#).unwrap();
        println!("{config:?}");
    }
}
