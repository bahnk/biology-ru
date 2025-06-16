use clap::{Parser, Subcommand};
use tracing_subscriber;
use rayon::ThreadPoolBuilder;

#[derive(Parser)]
#[command(name = "fastq")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    #[command(subcommand)]
    Uaspire(uaspire::Command),
}

mod uaspire {
///////////////////////////////////////////////////////////////////////////////

use std::fs::File;
use std::io::BufReader;
use std::sync::{Arc, Mutex};
use std::collections::{HashMap, HashSet};
use csv::Writer;
use bio::io::fastq;
use flate2::read::MultiGzDecoder;
use clap::Parser;
use clap::Subcommand;
use rayon::prelude::*;
use tracing::{info, error};

#[derive(Subcommand, Clone)]
pub enum Command {
    Create(CreateFastqCommand),
    Parse(ParseFastqCommand),
}

pub fn command(cmd: Command) {
    match cmd {
        Command::Create(cmd) => {
            println!("Creating FASTQ file with name: {}", cmd.path.display());
        }
        Command::Parse(cmd) => {
            process_fastq(
                &cmd.read1.to_string_lossy(),
                &cmd.read2.to_string_lossy(),
                1_000_000,
            );
        }
    }
}

fn count_records(path1: &str, path2: &str) -> usize {
    let file1 = File::open(path1).unwrap();
    let file2 = File::open(path2).unwrap();

    let reader1 = fastq::Reader::new(MultiGzDecoder::new(file1));
    let reader2 = fastq::Reader::new(MultiGzDecoder::new(file2));

    let count1 = reader1.records().count();
    let count2 = reader2.records().count();

    if count1 != count2 {
        error!("The number of records in the two files do not match: {} vs {}", count1, count2);
    }

    count1
}

fn process_pairs(
    rec1: &fastq::Record,
    rec2: &fastq::Record,
    results: &mut HashMap<(String, String), u64>,
    rbs_counts: &mut HashMap<(String, String), HashMap<String, u64>>,
) {
    let barcode_list: HashSet<&str> = [
        "ATCACG", "CGATGT", "CTTGTA", "GCCAAT", "ACAGTG", "ACTTGA",
    ].into();
    ///////////////////////////////////////////////////////////////////////////
    let barcode_len = 6;
    let discriminator_offset = 6;
    ///////////////////////////////////////////////////////////////////////////

    let id1 = rec1.id();
    let id2 = rec2.id();

    if id1 != id2 {
        error!("Record IDs do not match: {} vs {}", id1, id2);
    }

    ///////////////////////////////////////////////////////////////////////////
    let seq2 = std::str::from_utf8(rec2.seq()).unwrap();
    let start_pos = seq2.find("GAGCTCGCAT");
    if start_pos.is_none() {
        return;
    }
    let start_pos = start_pos.unwrap();
    if start_pos < barcode_len {
        return;
    }

    let barcode2 = (&seq2[start_pos - barcode_len..start_pos]).to_string();

    ///////////////////////////////////////////////////////////////////////////
    let seq1 = std::str::from_utf8(rec1.seq()).unwrap();

    let start_pos_nf = seq1.find("GGGTTTGTACCGTACAC");  // Non-flipped
    let start_pos_ff = seq1.find("GCCCGGATGATCCTGAC");  // Flipped

    if start_pos_nf.is_none() && start_pos_ff.is_none() {
        return;
    }

    let mut barcode1 = String::new();

    if ! start_pos_nf.is_none() {
        let start_pos_nf = start_pos_nf.unwrap();
        if start_pos_nf > (discriminator_offset + barcode_len) {
            let begin = start_pos_nf - barcode_len - discriminator_offset;
            let end = start_pos_nf - discriminator_offset;
            barcode1 = (&seq1[begin..end]).to_string();
        }
    }

    if ! start_pos_ff.is_none() {
        let start_pos_ff = start_pos_ff.unwrap();
        if start_pos_ff > (discriminator_offset + barcode_len) {
            let begin = start_pos_ff - barcode_len - discriminator_offset;
            let end = start_pos_ff - discriminator_offset;
            barcode1 = (&seq1[begin..end]).to_string();
        }
    }

    let pair = (barcode1.clone(), barcode2.clone());

    ///////////////////////////////////////////////////////////////////////////
    let rbs_length = 17;

    if (
        barcode_list.contains(barcode1.as_str()) &&
        barcode_list.contains(barcode2.as_str())
    ) {
        let rbs = (&seq2[start_pos+10..start_pos+10+rbs_length]).to_string();

        rbs_counts
            .entry(pair.clone())
            .or_insert_with(HashMap::new)
            .entry(rbs)
            .and_modify(|count| *count += 1)
            .or_insert(1);

    }
    ///////////////////////////////////////////////////////////////////////////

    if results.contains_key(&pair) {
        let count = results.get_mut(&pair).unwrap();
        *count += 1;
    } else {
        results.insert(pair, 1);
    }
}

fn process_records(
    rec1: &fastq::Record,
    rec2: &fastq::Record,
    results: &mut HashMap<String, u64>,
) {
    let id1 = rec1.id();
    let id2 = rec2.id();

    if id1 != id2 {
        error!("Record IDs do not match: {} vs {}", id1, id2);
    }

    let seq2 = std::str::from_utf8(rec2.seq()).unwrap();


    let start_pos = seq2.find("GAGCTCGCAT");

    if start_pos.is_none() {
        return;
    }


    let barcode_len = 6;

    let start_pos = start_pos.unwrap();

    if start_pos < barcode_len {
        return;
    }

    let barcode = &seq2[start_pos - barcode_len..start_pos];

    if results.contains_key(barcode) {
        let count = results.get_mut(barcode).unwrap();
        *count += 1;
    } else {
        results.insert(barcode.to_string(), 1);
    }
}

fn process_fastq(path1: &str, path2: &str, chunk_size: usize) {
    info!("Processing FASTQ files: {} and {}", path1, path2);

    let file1 = BufReader::new(File::open(path1).unwrap());
    let file2 = BufReader::new(File::open(path2).unwrap());

    let mut reader1 = fastq::Reader::new(MultiGzDecoder::new(file1));
    let mut reader2 = fastq::Reader::new(MultiGzDecoder::new(file2));

    let mut iter1 = reader1.records();
    let mut iter2 = reader2.records();

    let results = Arc::new(Mutex::new(HashMap::new()));
    let rbs_counts:
        Arc<Mutex<HashMap<(String, String), HashMap<String, u64>>>>
        = Arc::new(Mutex::new(HashMap::new()));

    let mut n = 0;

    loop {

        info!("Processing {}", n);

        let chunk1: Vec<_> = iter1
            .by_ref()
            .take(chunk_size)
            .collect();

        let chunk2: Vec<_> = iter2
            .by_ref()
            .take(chunk_size)
            .collect();


        if chunk1.is_empty() || chunk2.is_empty() {
            info!("No more records to process.");
            break;
        }

        chunk1
            .par_iter()
            .zip(chunk2.par_iter())
            .for_each(|(rec1, rec2)| {
                let mut results_lock = results.lock().unwrap();
                let mut rbs_counts_lock = rbs_counts.lock().unwrap();
                //process_records(
                //    &rec1.as_ref().unwrap(),
                //    &rec2.as_ref().unwrap(),
                //    &mut results_lock,
                //);
                process_pairs(
                    &rec1.as_ref().unwrap(),
                    &rec2.as_ref().unwrap(),
                    &mut results_lock,
                    &mut rbs_counts_lock,
                );

            });

        n += chunk1.len();
    }

    let mut sorted_results = results.lock().unwrap().clone().into_iter()
        .collect::<Vec<_>>();
    sorted_results.sort_by(|a, b| b.1.cmp(&a.1));

    //let mut wtr = Writer::from_path("results.csv")
    //    .expect("Failed to create output CSV file");
    //wtr.write_record(&["Sequence", "Count"])
    //    .expect("Failed to write CSV header");

    //for count in &sorted_results {
    //    wtr.write_record(&[&count.0.to_string(), &count.1.to_string()])
    //        .expect("Failed to write record to CSV");
    //}

    let mut wtr = Writer::from_path("results.csv")
        .expect("Failed to create output CSV file");
    wtr.write_record(&["Barcode 1", "Barcode 2", "Count"])
        .expect("Failed to write CSV header");

    for count in &sorted_results {
        wtr.write_record(&[
            &count.0.0.to_string(),
            &count.0.1.to_string(),
            &count.1.to_string()
            ])
            .expect("Failed to write record to CSV");
    }

    wtr.flush().expect("Failed to flush CSV writer");

    for ((barcode1, barcode2), rbs_map) in rbs_counts.lock().unwrap().iter() {
        let mut wtr = Writer::from_path(format!("tmp/rbs/{}_{}.csv", barcode1, barcode2))
            .expect("Failed to create output CSV file for RBS counts");
        wtr.write_record(&["RBS", "Count"])
            .expect("Failed to write RBS CSV header");

        for (rbs, count) in rbs_map {
            wtr.write_record(&[rbs, &count.to_string()])
                .expect("Failed to write RBS record to CSV");
        }

        wtr.flush().expect("Failed to flush RBS CSV writer");
    }


}

#[derive(Parser, Clone, Debug)]
pub struct CreateFastqCommand {
    #[arg(short, long)]
    path: std::path::PathBuf,
}

#[derive(Parser, Clone, Debug)]
pub struct ParseFastqCommand {
    #[arg()]
    read1: std::path::PathBuf,
    #[arg()]
    read2: std::path::PathBuf,
}

///////////////////////////////////////////////////////////////////////////////
}


fn main() {
    tracing_subscriber::fmt()
        .compact()
        .with_max_level(tracing::Level::INFO)
        .init();

    let cli = Cli::parse();

    ThreadPoolBuilder::new()
        .num_threads(10)
        .build_global()
        .expect("Failed to build thread pool");

    match &cli.command {
        Commands::Uaspire(cmd) => { uaspire::command(cmd.clone()); }
    }
}
