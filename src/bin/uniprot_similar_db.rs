
use biology_ru::uniprot::models::*;

use regex::Regex;
use reqwest::blocking::{get, Response};
use std::error::Error;

fn fetch_and_parse(url: &str) -> Result<Vec<String>, Box<dyn Error>> {

    let text = get(url)?.text()?;

    let lines: Vec<String> = text
        .lines()
        .map(|line| line.to_string())
        .collect();

    Ok(lines)
}

fn find_first_line(lines: &Vec<String>) -> usize {
    16 as usize
}

fn find_last_line(lines: &Vec<String>) -> usize {
    lines.len() - 5
}

fn get_line_range(lines: &Vec<String>) -> (usize, usize) {
    (find_first_line(lines), find_last_line(lines))
}

fn main() -> Result<(), Box<dyn Error>> {
    let url = "http://www.uniprot.org/docs/similar.txt";

    let lines = fetch_and_parse(url)?;

    let (first_line, last_line) = get_line_range(&lines);

    println!("first line: {}, and last line {}", first_line, last_line);

    let family_pattern = Regex::new(r"(^\S.*)").unwrap();
    let entry_pattern = Regex::new(
        r"^(?P<entry_name>[^(]+)\((?P<accession_number>[^)]+)\)$"
    ).unwrap();

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

                println!("{}, {}, {}", family, entry_name, accession_number);

                let entry = SimilarEntry {
                    family: family.clone(),
                    entry_name: entry_name.clone(),
                    accession_number: accession_number.clone(),
                };
            }
        }




    }

    Ok(())
}
