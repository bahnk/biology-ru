use biology_ru::uniprot::models::*;
use biology_ru::uniprot::similar::get_similar_entries;
use biology_ru::schema::uniprot_families;
use biology_ru::schema::uniprot_entries;

use std::env;
use env_logger;
use log::{info, warn};
use regex::Regex;
use dotenvy::dotenv;
use diesel::prelude::*;

pub fn establish_connection() -> SqliteConnection {
    dotenv().ok();

    let database_url = env::var("DATABASE_URL").expect("DATABASE_URL must be set");

    SqliteConnection::establish(&database_url)
        .unwrap_or_else(|_| panic!("Error connecting to {}", database_url))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();

    let url = "http://www.uniprot.org/docs/similar.txt";

    let entries = get_similar_entries(url)?;
    let connection = &mut establish_connection();

    let species_pattern = Regex::new(r"^[^_]+_(HUMAN|MOUSE)$")
        .expect("Error compiling species regex");

    for entry in entries {
        let (family, entry_name, accession_number) = entry;

        if ! species_pattern.is_match(&entry_name) {
            continue;
        }

        info!("Processing {entry_name}");

        let similar_family = SimilarFamily{ name: family.clone() };

        diesel::insert_into(uniprot_families::table)
            .values(&similar_family)
            .on_conflict(uniprot_families::name)
            .do_update()
            .set(uniprot_families::name.eq(similar_family.name.clone()))
            .execute(connection)
            .expect("Error inserting new family");

        let similar_entry = SimilarEntry{
            family: family.clone(),
            entry_name: entry_name.clone(),
            accession_number: accession_number.clone(),
        };

        diesel::insert_into(uniprot_entries::table)
            .values(&similar_entry)
            .on_conflict(uniprot_entries::accession_number)
            .do_update()
            .set(
                uniprot_entries::accession_number
                .eq(similar_entry.accession_number.clone())
            )
            //.set(entry_name.eq(similar_entry.entry_name.clone()))
            //.set(family.eq(similar_entry.family.clone()))
            .execute(connection)
            .expect("Error inserting new family");



    }

    Ok(())
}
