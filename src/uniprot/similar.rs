use thiserror::Error;
use regex::Regex;
//use reqwest::Error;
use reqwest::blocking::{get, Response};

#[derive(Error, Debug)]
pub enum EntryError {
    //#[error("Network error occured while fetching data: {0}")]
    //NetworkError(#[from] reqwest::Error),

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

pub fn get_similar_entries(url: &str) -> Result<Vec<(String, String, String)>, EntryError> {

    let lines = fetch_and_parse(url);

    // Need to improve error handling
    if lines.len() < 17 {
        return Err(EntryError::InsufficientData)
    }

    let (first_line, last_line) = get_line_range(&lines);

    let family_pattern = Regex::new(r"(^\S.*)")?;
    let entry_pattern = Regex::new(r"^(?P<entry_name>[^(]+)\((?P<accession_number>[^)]+)\)$")?;

    let mut entries: Vec<(String, String, String)> = Vec::new();
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

                entries.push((
                    family.clone(), entry_name.clone(), accession_number.clone()
                ));
            }
        }
    }

    Ok(entries)
}
