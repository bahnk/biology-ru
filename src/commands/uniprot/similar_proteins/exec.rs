use std::fs;
use polars::prelude::*;

use crate::config::global::GlobalConfig;

use crate::uniprot::similar_proteins::{
    fetch_similar_proteins,
    entries_to_dataframe,
};

use super::config::Config;

pub fn execute(cfg: Config, global_cfg: GlobalConfig) {
    let entries = fetch_similar_proteins(&cfg.url).unwrap();
    let mut df = entries_to_dataframe(&entries).unwrap();
    let output_dir = global_cfg.hive_tables_dir.join(cfg.table_name);
    let _ = fs::create_dir_all(&output_dir);
    let path = output_dir.join("part-0.parquet");
    let file = fs::File::create(path).unwrap();
    let _ = ParquetWriter::new(file).finish(&mut df);
}

#[cfg(test)]
mod tests {
    use super::execute;
    use super::Config;
    use super::GlobalConfig;

    #[test]
    fn test_execute() {
        let global_cfg = GlobalConfig{
            log_level: "info".into(),
            hive_tables_dir: "/tmp/uniprot/parquet".into(),
            threads: 4,
        };
        let cfg = Config{
            url: "http://www.uniprot.org/docs/similar.txt".into(),
            table_name: "similar_proteins".into(),
        };
        execute(cfg, global_cfg);
    }
}
