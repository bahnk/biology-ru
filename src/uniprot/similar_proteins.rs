/// Uniprot Protein Families List
use log::info;
use polars::prelude::*;
use regex::Regex;
use reqwest::blocking;
use thiserror::Error;

// ----------------------------------------------------------------------------
// Errors Definition
// ----------------------------------------------------------------------------

#[derive(Debug, Error)]
pub enum SimilarProtError {
    #[error("Cannot get URL")]
    Get,
    #[error("Cannot decode URL content")]
    Decode,
    #[error("Regex pattern error: {0}")]
    Regex(#[from] regex::Error),
    #[error("Insufficient data: fewer than header")]
    InsufficientData,
    #[error("Error parsing entry information")]
    Parse,
    #[error("Polars error: {0}")]
    Polars(#[from] PolarsError),
}

// ----------------------------------------------------------------------------
// Fetching Lines Function
// ----------------------------------------------------------------------------

fn fetch_lines(url: &str) -> Result<Vec<String>, SimilarProtError> {
    let response: blocking::Response;
    let text: String;

    info!("Fetching content of {}", url);

    match blocking::get(url) {
        Ok(resp) => response = resp,
        Err(_) => return Err(SimilarProtError::Get),
    }

    match response.text() {
        Ok(txt) => text = txt,
        Err(_) => return Err(SimilarProtError::Decode),
    }

    Ok(text.lines().map(|line| line.to_string()).collect())
}

// ---------------------------------------------------------------------------
// Fetch Similar Protein Entries
// ----------------------------------------------------------------------------

#[derive(Debug)]
pub struct Entry {
    #[allow(unused)]
    family: String,
    #[allow(unused)]
    entry_name: String,
    #[allow(unused)]
    accession_number: String,
}

pub fn fetch_similar_proteins(
    url: &str,
) -> Result<Vec<Entry>, SimilarProtError> {
    let lines = fetch_lines(url)?;

    let family_pattern = Regex::new(r"(^\S.*)")?;
    let entry_pattern =
        Regex::new(r"^(?P<entry_name>[^(]+)\((?P<accession_number>[^)]+)\)$")?;

    let mut entries: Vec<Entry> = Vec::new();
    let mut family = String::new();

    for line in lines.into_iter().skip(16) {
        if line.is_empty() {
            continue;
        }

        if let Some(caps) = family_pattern.captures(&line) {
            family = caps
                .get(0)
                .map(|m| m.as_str().to_string())
                .unwrap_or_else(|| "".to_string());

            continue;
        }

        for entry in line.replace(" ", "").split(",") {
            if entry.is_empty() {
                continue;
            }

            if let Some(caps) = entry_pattern.captures(entry) {
                let entry_name = caps
                    .name("entry_name")
                    .map(|m| m.as_str().to_string())
                    .unwrap_or_else(|| "".to_string());

                let accession_number = caps
                    .name("accession_number")
                    .map(|m| m.as_str().to_string())
                    .unwrap_or_else(|| "".to_string());

                entries.push(Entry {
                    family: family.clone(),
                    entry_name: entry_name.clone(),
                    accession_number: accession_number.clone(),
                });
            }
        }
    }

    info!("Fetched {} entries", entries.len());

    Ok(entries)
}

pub fn entries_to_dataframe(
    entries: &Vec<Entry>,
) -> Result<DataFrame, SimilarProtError> {
    let mut families: Vec<String> = Vec::new();
    let mut species: Vec<String> = Vec::new();
    let mut entry_names: Vec<String> = Vec::new();
    let mut accession_numbers: Vec<String> = Vec::new();

    let species_pattern = Regex::new(r"^[^_]+_(?P<species>[^_]+)$")?;

    for entry in entries.into_iter() {
        if let Some(caps) = species_pattern.captures(&entry.entry_name) {
            let species_code = caps
                .name("species")
                .map(|m| m.as_str().to_string())
                .unwrap_or_else(|| "".to_string());

            families.push(entry.family.clone());
            species.push(species_code);
            entry_names.push(entry.entry_name.clone());
            accession_numbers.push(entry.accession_number.clone());
        }
    }

    let df: DataFrame;

    match df!(
        "family" => families,
        "species" => species,
        "entry_name" => entry_names,
        "accession_number" => accession_numbers,
    ) {
        Ok(dataframe) => df = dataframe,
        Err(e) => return Err(SimilarProtError::Polars(e)),
    }

    Ok(df)
}

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::entries_to_dataframe;
    use super::fetch_lines;
    use super::fetch_similar_proteins;

    #[test]
    fn fetch_similarity_list_lines() {
        let url = "http://www.uniprot.org/docs/similar.txt";
        let lines = fetch_lines(url).unwrap();
        for line in lines.into_iter().take(20) {
            println!("{line}");
        }
    }

    #[test]
    fn fetch_similar_proteins_url() {
        let url = "http://www.uniprot.org/docs/similar.txt";
        let entries = fetch_similar_proteins(url).unwrap();
        println!("Number of entries: {}", entries.len());
        for entry in entries.into_iter().take(5) {
            println!("{entry:?}");
        }
    }

    #[test]
    fn create_entries_dataframe() {
        let url = "http://www.uniprot.org/docs/similar.txt";
        let entries = fetch_similar_proteins(url).unwrap();
        let df = entries_to_dataframe(&entries);
        println!("{df:?}");
    }
}
