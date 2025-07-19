use diesel::prelude::*;
use log::info;
use regex::Regex;
use reqwest::blocking::{get, Response};
use thiserror::Error;

use crate::schema::*;
use crate::uniprot::models::*;

#[derive(Error, Debug)]
pub enum EntryError {
    #[error("Regex pattern error: {0}")]
    RegexError(#[from] regex::Error),

    #[error("Insufficient data: fewer than header")]
    InsufficientData,

    #[error("Error parsing entry information")]
    ParseError,
}

fn fetch_and_parse(url: &str) -> Vec<String> {
    let lines: Vec<String> = Vec::new();
    let response: Response;
    let text: String;

    match get(url) {
        Ok(resp) => response = resp,
        Err(_) => return lines,
    }

    match response.text() {
        Ok(txt) => text = txt,
        Err(_) => return lines,
    }

    text.lines().map(|line| line.to_string()).collect()
}

fn find_first_line(_lines: &Vec<String>) -> usize {
    16 as usize
}

fn find_last_line(lines: &Vec<String>) -> usize {
    lines.len() - 5
}

fn get_line_range(lines: &Vec<String>) -> (usize, usize) {
    (find_first_line(lines), find_last_line(lines))
}

pub fn get_similar_entries(
    url: &str,
) -> Result<Vec<(UniprotFamily, UniprotEntry)>, EntryError> {
    let lines = fetch_and_parse(url);

    // Need to improve error handling
    if lines.len() < 17 {
        return Err(EntryError::InsufficientData);
    }

    let (first_line, last_line) = get_line_range(&lines);

    let family_pattern = Regex::new(r"(^\S.*)")?;
    let entry_pattern =
        Regex::new(r"^(?P<entry_name>[^(]+)\((?P<accession_number>[^)]+)\)$")?;

    let mut entries: Vec<(UniprotFamily, UniprotEntry)> = Vec::new();
    let mut family = String::new();

    for line in &lines[first_line..last_line] {
        if line.is_empty() {
            continue;
        }

        if let Some(caps) = family_pattern.captures(line) {
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

                let family = UniprotFamily {
                    name: family.clone(),
                };

                let entry = UniprotEntry {
                    entry_name: entry_name.clone(),
                    accession_number: accession_number.clone(),
                    mass: None,
                    seq_length: None,
                };
                entries.push((family, entry));
            }
        }
    }

    Ok(entries)
}

pub fn filter_by_species(
    entries: &[(UniprotFamily, UniprotEntry)],
    species: &[String],
) -> Result<Vec<(UniprotFamily, UniprotEntry)>, regex::Error> {
    let joined = species.join("|");
    let regex_pattern = format!("[^_]+_({})$", joined);
    let regex = Regex::new(&regex_pattern)?;

    let mut selected_entries: Vec<(UniprotFamily, UniprotEntry)> = Vec::new();

    for (family, entry) in entries {
        if regex.is_match(&entry.entry_name) {
            selected_entries.push((family.clone(), entry.clone()));
        }
    }

    Ok(selected_entries)
}

pub fn insert_entries(
    entries: &[(UniprotFamily, UniprotEntry)],
    connection: &mut SqliteConnection,
) -> Result<(), Box<dyn std::error::Error>> {
    info!("Starting to insert {} entries", entries.len());

    for (index, (family, entry)) in entries.iter().enumerate() {
        if index % 1000 == 0 {
            info!("Inserted {0}/{1}", index, entries.len());
        }

        diesel::insert_into(uniprot_sequence_similarity_families::table)
            .values(family)
            .on_conflict(uniprot_sequence_similarity_families::name)
            .do_update()
            .set(
                uniprot_sequence_similarity_families::name
                    .eq(family.name.clone()),
            )
            .execute(connection)?;

        diesel::insert_into(uniprot_entries::table)
            .values(entry)
            .on_conflict(uniprot_entries::accession_number)
            .do_update()
            .set((
                uniprot_entries::accession_number
                    .eq(entry.accession_number.clone()),
                uniprot_entries::entry_name.eq(entry.entry_name.clone()),
                uniprot_entries::mass.eq(entry.mass.clone()),
                uniprot_entries::seq_length.eq(entry.seq_length.clone()),
            ))
            .execute(connection)?;

        let belongs_to = BelongsToFamily {
            entry: entry.accession_number.clone(),
            family: family.name.clone(),
        };

        diesel::insert_into(
            belongs_to_uniprot_sequence_similarity_family::table,
        )
        .values(&belongs_to)
        .on_conflict((
            belongs_to_uniprot_sequence_similarity_family::entry,
            belongs_to_uniprot_sequence_similarity_family::family,
        ))
        .do_update()
        .set((
            belongs_to_uniprot_sequence_similarity_family::entry
                .eq(belongs_to.entry.clone()),
            belongs_to_uniprot_sequence_similarity_family::family
                .eq(belongs_to.family.clone()),
        ))
        .execute(connection)?;
    }

    info!("Finishing inserting all entries");
    Ok(())
}
