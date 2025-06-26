/// This module processes FASTQ files to count barcode pairs and RBS sequences.

use std::fs::File;
use std::path::Path;
use std::io::BufReader;
use std::sync::{Arc, Mutex};
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use csv::Writer;
use bio::io::fastq;
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use tracing::{info, error};

use crate::uaspire::constants;

/// Counts the number of reads that pass various quality checks.
#[derive(Debug)]
struct ReadCounts {
    basecall_failed: u64,
    constant_region_sequence_failed: u64,
    constant_region_position_failed: u64,
    barcode_1_failed: u64,
    barcode_2_failed: u64,
    discriminator_sequence_failed: u64,
    discriminator_position_failed: u64,
    valid: u64,
    total: u64,
}

/// Represents a sample with two barcodes.
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
struct Sample {
    barcode_1: String,
    barcode_2: String,
}

/// Counts the number of flips for a given pair of barcodes.
#[derive(Debug, Clone)]
struct FlipCounts {
    non_flipped_count: u64,
    flipped_count: u64,
}

/// A type alias for a HashMap that maps barcode pairs to their flip counts.
type RBSCounts<T, U> = HashMap<T, HashMap<U, FlipCounts>>;

/// Counts the number of 'N' characters in a sequence.
fn count_n(sequence: &str) -> u64 {
    sequence.chars().filter(|&c| c == 'N').count() as u64
}
    

fn process_pairs(
    rec1: &fastq::Record,
    rec2: &fastq::Record,
    reads_counts: &mut ReadCounts,
    rbs_counts: &mut RBSCounts<Sample, String>,
) {
    reads_counts.total += 1;

    let id1 = rec1.id();
    let id2 = rec2.id();

    if id1 != id2 {
        error!("Record IDs do not match: {} vs {}", id1, id2);
    }

    let seq1 = std::str::from_utf8(rec1.seq()).unwrap();
    let seq2 = std::str::from_utf8(rec2.seq()).unwrap();

    if count_n(seq1) + count_n(seq2) > constants::MAX_N_COUNT as u64 {
        reads_counts.basecall_failed += 1;
        return;
    }

    let (win_start, win_end) = constants::CONSTANT_REGION_WINDOW;
    let window = &seq2[(win_start as usize - 1)..(win_end as usize)];
    let local_start = window.find(constants::CONTSTANT_REGION);
    if local_start.is_none() {
        reads_counts.constant_region_sequence_failed += 1;
        return;
    }
    let const_reg_offset = local_start.unwrap() + win_start as usize - 1;
    if const_reg_offset < constants::BARCODE_LEN {
        reads_counts.constant_region_position_failed += 1;
        return;
    }
    if const_reg_offset + constants::BARCODE_LEN + constants::RBS_LEN > seq2.len() {
        reads_counts.constant_region_position_failed += 1;
        return;
    }

    let rbs_start = const_reg_offset + constants::CONTSTANT_REGION.len();
    let rbs = &seq2[(rbs_start)..(rbs_start + constants::RBS_LEN)];

    let barcode2 = &seq2[(const_reg_offset - constants::BARCODE_LEN)..const_reg_offset];

    if !constants::BARCODES_2.contains(&barcode2) {
        reads_counts.barcode_2_failed += 1;
        return;
    }

    let start_nf = seq1.find(constants::NON_FLIPPED_SEQ);
    let start_ff = seq1.find(constants::FLIPPED_SEQ);

    let discriminator_start: usize;
    let flipped: bool;
    if let Some(start_nf) = start_nf {
        discriminator_start = start_nf;
        flipped = false;
    } else if let Some(start_ff) = start_ff {
        discriminator_start = start_ff;
        flipped = true;
    } else {
        reads_counts.discriminator_sequence_failed += 1;
        return;
    }

    if discriminator_start < (constants::DISCRIMINATOR_OFFSET + constants::BARCODE_LEN) {
        reads_counts.discriminator_position_failed += 1;
        return;
    }

    let start_bcd_1 = discriminator_start - constants::BARCODE_LEN - constants::DISCRIMINATOR_OFFSET;
    let end_bcd_1 = discriminator_start - constants::DISCRIMINATOR_OFFSET;
    let barcode1 = &seq1[start_bcd_1..end_bcd_1];

    if !constants::BARCODES_1.contains(&barcode1) {
        reads_counts.barcode_1_failed += 1;
        return;
    }

    reads_counts.valid += 1;

    let sample = Sample {
        barcode_1: barcode1.to_string(),
        barcode_2: barcode2.to_string(),
    };

    let counts = rbs_counts
        .entry(sample.clone())
        .or_insert_with(|| HashMap::new());

    let flip_counts = counts
        .entry(rbs.to_string())
        .or_insert_with(|| FlipCounts {
            non_flipped_count: 0,
            flipped_count: 0,
        });

    if flipped {
        flip_counts.flipped_count += 1;
    } else {
        flip_counts.non_flipped_count += 1;
    }
}

pub fn process_fastq(path1: &str, path2: &str, chunk_size: usize, output_dir: &str) {

    info!("Creating output directories if they do not exist");

    if !std::path::Path::new(output_dir).exists() {
        std::fs::create_dir_all(output_dir).expect("Failed to create output directory");
    }
    let samples_dir = Path::new(output_dir).join("samples");
    if !samples_dir.exists() {
        std::fs::create_dir_all(&samples_dir).expect("Failed to create samples directory");
    }

    info!("Processing FASTQ files: {} and {}", path1, path2);

    let file1 = BufReader::new(File::open(path1).unwrap());
    let file2 = BufReader::new(File::open(path2).unwrap());

    let reader1 = fastq::Reader::new(MultiGzDecoder::new(file1));
    let reader2 = fastq::Reader::new(MultiGzDecoder::new(file2));

    let mut iter1 = reader1.records();
    let mut iter2 = reader2.records();

    let read_counts = Arc::new(Mutex::new(ReadCounts {
        basecall_failed: 0,
        constant_region_sequence_failed: 0,
        constant_region_position_failed: 0,
        barcode_1_failed: 0,
        barcode_2_failed: 0,
        discriminator_sequence_failed: 0,
        discriminator_position_failed: 0,
        valid: 0,
        total: 0,
    }));

    let rbs_counts: Arc<Mutex<RBSCounts<Sample, String>>> = Arc::new(Mutex::new(HashMap::new()));

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

                let mut read_counts_lock = read_counts.lock().unwrap();
                let mut rbs_counts_lock = rbs_counts.lock().unwrap();

                process_pairs(
                    &rec1.as_ref().unwrap(),
                    &rec2.as_ref().unwrap(),
                    &mut read_counts_lock,
                    &mut rbs_counts_lock,
                );

            });

        n += chunk1.len();
    }

    let read_counts = Arc::try_unwrap(read_counts).unwrap().into_inner().unwrap();
    let read_counts_csv = Path::new(output_dir).join("read_counts.csv");
    let mut wtr = Writer::from_path(read_counts_csv)
        .expect("Failed to create read counts CSV file");
    wtr.write_record(&["Read Type", "Count"])
        .expect("Failed to write CSV header");
    wtr.write_record(&["Basecall Fail", &read_counts.basecall_failed.to_string()])
        .expect("Failed to write basecall failed record");
    wtr.write_record(&["Constant Region Sequence Fail", &read_counts.constant_region_sequence_failed.to_string()])
        .expect("Failed to write constant region sequence failed record");
    wtr.write_record(&["Constant Region Position Fail", &read_counts.constant_region_position_failed.to_string()])
        .expect("Failed to write constant region position failed record");
    wtr.write_record(&["Barcode 1 Fail", &read_counts.barcode_1_failed.to_string()])
        .expect("Failed to write barcode 1 failed record");
    wtr.write_record(&["Barcode 2 Fail", &read_counts.barcode_2_failed.to_string()])
        .expect("Failed to write barcode 2 failed record");
    wtr.write_record(&["Discriminator Sequence Fail", &read_counts.discriminator_sequence_failed.to_string()])
        .expect("Failed to write discriminator sequence failed record");
    wtr.write_record(&["Discriminator Position Fail", &read_counts.discriminator_position_failed.to_string()])
        .expect("Failed to write discriminator position failed record");
    wtr.write_record(&["Valid Reads", &read_counts.valid.to_string()])
        .expect("Failed to write valid reads record");
    wtr.write_record(&["Total Reads", &read_counts.total.to_string()])
        .expect("Failed to write total reads record");
    wtr.flush().expect("Failed to flush CSV writer");

    for (sample, counts) in rbs_counts.lock().unwrap().iter() {
        let sample_csv = samples_dir.join(format!("{}_{}.csv", sample.barcode_1, sample.barcode_2));
        let mut wtr = Writer::from_path(sample_csv)
            .expect("Failed to create RBS counts CSV file");
        wtr.write_record(&["RBS", "Non-Flipped Count", "Flipped Count"])
            .expect("Failed to write RBS CSV header");

        for (rbs, flip_counts) in counts.iter() {
            wtr.write_record(&[
                rbs,
                &flip_counts.non_flipped_count.to_string(),
                &flip_counts.flipped_count.to_string(),
            ])
            .expect("Failed to write RBS record to CSV");
        }
        wtr.flush().expect("Failed to flush RBS CSV writer");
    }

    info!("Processing complete.");
}
