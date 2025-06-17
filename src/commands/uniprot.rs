
use clap::{Parser, Subcommand};
use config::builder::DefaultState;
use config::{Config, ConfigBuilder, Environment, File};
use diesel::prelude::*;
use dotenvy::dotenv;
use env_logger;
use std::path::PathBuf;

use crate::uniprot::similar::{
    filter_by_species,
    get_similar_entries,
    insert_entries
};

///////////////////////////////////////////////////////////////////////////////

#[derive(Subcommand, Debug)]
pub enum Commands {
    Uniprot(Args),
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    #[arg(short, long, default_value = "assets/config")]
    config: PathBuf,
}

///////////////////////////////////////////////////////////////////////////////

#[allow(dead_code)]
fn establish_connection(settings: &Config) -> Result<SqliteConnection, Box<dyn std::error::Error>> {
    let database_url: String = settings.get("DATABASE_URL")?;
    let connection = SqliteConnection::establish(&database_url)
        .map_err(|e| format!("Error connecting to {}: {}", database_url, e))?;
    Ok(connection)
}

#[allow(dead_code)]
fn load_settings(args: &Args) -> Result<Config, Box<dyn std::error::Error>> {
    dotenv().ok();
    let config_file = args.config.to_str().ok_or("Invalid config path")?;
    let settings = ConfigBuilder::<DefaultState>::default()
        .add_source(File::with_name(config_file))
        .add_source(Environment::default())
        .build()?;

    Ok(settings)
}

pub fn command(_cmds: Commands) {
    println!("Running Uniprot command");
}

fn run(args: &Args) -> Result<(), Box<dyn std::error::Error>> {
    // Configuration
    let settings = load_settings(&args)?;
    let url: String = settings.get("uniprot.similar.url")?;
    let species: Vec<String> = settings.get("uniprot.similar.species")?;
    let mut connection = establish_connection(&settings)?;

    // Process entries
    let all_entries = get_similar_entries(&url)?;
    let entries = filter_by_species(&all_entries, &species)?;
    insert_entries(&entries, &mut connection)?;

    Ok(())
}

#[allow(dead_code)]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let args = Args::parse();
    run(&args)?;
    Ok(())
}
