/// This module processes FASTQ files to count barcode pairs and RBS sequences.
use bio::io::fastq;
use dashmap::DashMap;
use flate2::read::MultiGzDecoder;
use polars::prelude::*;
use rayon::prelude::*;
use tracing::{error, info};

use strum::EnumCount;
use strum_macros::EnumCount;

use std::{
    fs::{self, File},
    hash::Hash,
    io::{self, BufReader},
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
    sync::Arc,
};

use crate::uaspire::constants;

// ---------- Sample table ----------

type SampleTable = DashMap<Sample, DashMap<String, [AtomicU64; 2]>>;

// ---------- Configuration ----------

#[derive(Clone, Debug)]
struct Config<'a> {
    barcodes1: &'a [&'a str],
    barcodes2: &'a [&'a str],
    const_region: &'a str,
    window: (usize, usize),
    rbs_len: usize,
    barcode_len: usize,
    max_n: usize,
    non_flipped: &'a str,
    flipped: &'a str,
    disc_offset: usize,
}

// ---------- Directory layout ----------

#[derive(Debug)]
pub struct DirLayout {
    pub root: PathBuf,
    pub data: PathBuf,
    pub counts: PathBuf,
    pub qc: PathBuf,
    pub tmp: PathBuf,
    pub parquet: PathBuf,
}

// ---------- Reads classication ---------

#[derive(Debug, EnumCount)]
enum FailReason {
    BaseCalls,
    ConstantSeq,
    ConstantPos,
    Barcode1,
    Barcode2,
    DiscSeq,
    DiscPos,
}

// ---------- Discriminator status ----------

#[derive(Debug)]
enum Flip {
    NonFlipped,
    Flipped,
}

// ---------- Sample is a barcode pair ----------

#[derive(Debug, Clone, Hash, Eq, PartialEq)]
struct Sample {
    barcode1: String,
    barcode2: String,
}

// ---------- Read types counter ----------

#[derive(Default)]
struct Counters {
    total: AtomicU64,
    valid: AtomicU64,
    fails: [AtomicU64; FailReason::COUNT],
}

impl Counters {
    fn inc_total(&self) {
        self.total.fetch_add(1, Ordering::Relaxed);
    }
    fn inc_valid(&self) {
        self.valid.fetch_add(1, Ordering::Relaxed);
    }
    fn inc_fail(&self, r: FailReason) {
        self.fails[r as usize].fetch_add(1, Ordering::Relaxed);
    }

    fn to_dataframe(&self) -> Result<DataFrame, polars::error::PolarsError> {
        let names = Series::new(
            "name".into(),
            &[
                "total",
                "valid",
                "base_calls",
                "constant_seq",
                "constant_pos",
                "barcode_1",
                "barcode_2",
                "disc_seq",
                "disc_pos",
            ],
        );

        let values = Series::new(
            "value".into(),
            &[
                self.total.load(Ordering::Relaxed),
                self.valid.load(Ordering::Relaxed),
                self.fails[FailReason::BaseCalls as usize]
                    .load(Ordering::Relaxed),
                self.fails[FailReason::ConstantSeq as usize]
                    .load(Ordering::Relaxed),
                self.fails[FailReason::ConstantPos as usize]
                    .load(Ordering::Relaxed),
                self.fails[FailReason::Barcode1 as usize]
                    .load(Ordering::Relaxed),
                self.fails[FailReason::Barcode2 as usize]
                    .load(Ordering::Relaxed),
                self.fails[FailReason::DiscSeq as usize]
                    .load(Ordering::Relaxed),
                self.fails[FailReason::DiscPos as usize]
                    .load(Ordering::Relaxed),
            ],
        );

        DataFrame::new(vec![names.into(), values.into()]).map_err(|e| {
            error!("Failed to create DataFrame: {}", e);
            panic!("Failed to create DataFrame");
        })
    }
}

// =========================================================
// Helper functions
// =========================================================

/// Validates that the IDs of two FASTQ records match.
fn validate_pairs(rec1: &fastq::Record, rec2: &fastq::Record) -> bool {
    let id1 = rec1.id();
    let id2 = rec2.id();

    if id1 != id2 {
        panic!("Record IDs do not match: {} vs {}", id1, id2);
    }

    true
}

/// Converts a `SampleTable` to a Polars `DataFrame`.
fn table_to_dataframe(
    table: &SampleTable,
) -> Result<DataFrame, polars::error::PolarsError> {
    let rows: Vec<_> = table
        .iter()
        .flat_map(|sample_map| {
            let sample = sample_map.key().clone();

            sample_map
                .value()
                .iter()
                .map(|rbs_counts| {
                    let rbs = rbs_counts.key();
                    let counts = rbs_counts.value();
                    let unflipped = counts[0].load(Ordering::Relaxed);
                    let flipped = counts[1].load(Ordering::Relaxed);

                    (
                        sample.barcode1.clone(),
                        sample.barcode2.clone(),
                        rbs.clone(),
                        unflipped,
                        flipped,
                    )
                })
                .collect::<Vec<(String, String, String, u64, u64)>>()
        })
        .collect();

    let df = {
        let mut bc1: Series = rows.iter().map(|r| r.0.clone()).collect();
        let mut bc2: Series = rows.iter().map(|r| r.1.clone()).collect();
        let mut rbs: Series = rows.iter().map(|r| r.2.clone()).collect();
        let mut unflipped: Series = rows.iter().map(|r| r.3.clone()).collect();
        let mut flipped: Series = rows.iter().map(|r| r.4.clone()).collect();

        DataFrame::new(vec![
            bc1.rename("barcode1".into()).clone().into(),
            bc2.rename("barcode2".into()).clone().into(),
            rbs.rename("gre".into()).clone().into(),
            unflipped.rename("unflipped".into()).clone().into(),
            flipped.rename("flipped".into()).clone().into(),
        ])?
    };

    Ok(df)
}

/// Write a `DataFrame` to a Parquet file on disk.
fn write_parquet_chunk(df: &DataFrame, path: &str) -> Result<u64, PolarsError> {
    let mut df = df.clone();
    let file = File::create(path)?;
    ParquetWriter::new(file).finish(&mut df)
}

/// Create all required directories under `output_dir`.
fn prepare_dirs(output_dir: impl AsRef<Path>) -> io::Result<DirLayout> {
    let root = output_dir.as_ref().to_path_buf();
    let data = root.join("data");
    let counts = data.join("counts");
    let qc = data.join("qc");
    let tmp = root.join("tmp");
    let parquet = tmp.join("parquet");

    let all = [&root, &data, &counts, &qc, &tmp, &parquet];

    for dir in all {
        fs::create_dir_all(dir)?;
    }

    Ok(DirLayout {
        root,
        data,
        counts,
        qc,
        tmp,
        parquet,
    })
}

/// List all Parquet files in a directory.
fn list_parquet_files(dir: impl AsRef<Path>) -> std::io::Result<Vec<PathBuf>> {
    let mut files: Vec<PathBuf> = fs::read_dir(dir)?
        .filter_map(|entry| {
            let entry = entry.ok()?;
            let p = entry.path();
            (p.extension()?.eq_ignore_ascii_case("parquet")).then_some(p)
        })
        .collect();

    files.sort();
    Ok(files)
}

/// Concatenate all Parquet files in a directory into a single `DataFrame`.
fn concat_parquet_dir(dir: impl AsRef<Path>) -> DataFrame {
    let files = match list_parquet_files(&dir) {
        Ok(files) => files,
        Err(e) => {
            error!("Failed to list parquet files: {}", e);
            panic!("Failed to list parquet files");
        }
    };

    if files.is_empty() {
        panic!("No parquet files");
    }

    let mut dfs = Vec::with_capacity(files.len());
    for path in files {
        let file =
            std::fs::File::open(&path).expect("Cannot open parquet file");
        let df = ParquetReader::new(file)
            .finish()
            .expect("Cannot read parquet file");
        dfs.push(df.lazy());
    }

    concat(&dfs, UnionArgs::default())
        .expect("Cannot concatenate dataframes")
        .collect()
        .expect("Cannot finalise dataframe")
        .lazy()
        .group_by_stable([col("barcode1"), col("barcode2"), col("gre")])
        .agg([
            col("unflipped").sum().alias("unflipped"),
            col("flipped").sum().alias("flipped"),
        ])
        .collect()
        .expect("Cannot convert LazyFrame to DataFrame")
}

fn write_qc_parquet(
    df: &DataFrame,
    output_root: &Path,
    sample: &str,
) -> PolarsResult<()> {
    let output_dir = output_root.join(format!("sample={sample}"));
    fs::create_dir_all(&output_dir)?;
    let path = output_dir.join("part-0.parquet");
    let file = std::fs::File::create(&path)?;
    ParquetWriter::new(file)
        .with_compression(ParquetCompression::Zstd(None))
        .finish(&mut df.clone())?;

    Ok(())
}

/// Write a partitioned Parquet file for each unique barcode pair.
fn write_partitioned_parquet(
    df: &DataFrame,
    output_root: &Path,
    sample_name: &str,
    parquet_size: Option<usize>,
) -> PolarsResult<()> {

    let parquet_size = parquet_size.unwrap_or(10_000);

    let barcode1s = df.column("barcode1")?.str()?.unique()?;
    let barcode2s = df.column("barcode2")?.str()?.unique()?;

    for b1 in barcode1s.into_no_null_iter() {
        for b2 in barcode2s.into_no_null_iter() {
            let filtered = df
                .clone()
                .lazy()
                .filter(
                    col("barcode1")
                        .eq(lit(b1))
                        .and(col("barcode2").eq(lit(b2))),
                )
                .collect()?;

            // Not all combinations are supposed to be present
            if filtered.height() == 0 {
                continue;
            }

            let path = output_root
                .join(format!("sample={}", sample_name))
                .join(format!("barcode1={}", b1))
                .join(format!("barcode2={}", b2));

            save_by_chunks(&filtered, &path, Some(parquet_size))?;
        }
    }

    Ok(())
}

/// Save a DataFrame to Parquet files in chunks.
fn save_by_chunks(
    df: &DataFrame,
    output_dir: &Path,
    chunk_size: Option<usize>,
) -> PolarsResult<()> {
    let chunk_size = chunk_size.unwrap_or(10_000);

    let n_rows = df.height();
    let n_chunks = (n_rows + chunk_size - 1) / chunk_size;
    let n_digits = ((n_chunks as f64).log10().floor() as usize) + 1;

    fs::create_dir_all(output_dir)?;

    for i in 0..n_chunks {
        let start = i * chunk_size;
        let len = chunk_size.min(n_rows - start);

        let chunk = df.slice(start as i64, len);
        let filename = format!("part-{i:0n_digits$}.parquet");
        let path = output_dir.join(filename);
        let file = std::fs::File::create(&path)?;

        ParquetWriter::new(file)
            .with_compression(ParquetCompression::Zstd(None))
            .finish(&mut chunk.clone())?;
    }

    Ok(())
}

// =========================================================
// Core logic
// =========================================================

fn classify_pair<'a>(
    cfg: &Config<'a>,
    rec1: &'a fastq::Record,
    rec2: &'a fastq::Record,
) -> Result<Result<(Sample, &'a str, Flip), FailReason>, std::str::Utf8Error> {
    validate_pairs(rec1, rec2);

    let seq1 = std::str::from_utf8(rec1.seq())?;
    let seq2 = std::str::from_utf8(rec2.seq())?;

    // -----------------------------------------------------
    // 1. Fast rejection for base call
    // -----------------------------------------------------
    if seq1.bytes().filter(|&b| b == b'N').count()
        + seq2.bytes().filter(|&b| b == b'N').count()
        > cfg.max_n
    {
        return Ok(Err(FailReason::BaseCalls));
    }

    // -----------------------------------------------------
    // 2. Reject when missing constant region
    // -----------------------------------------------------
    let (win_lo, win_hi) = cfg.window;
    let window = &seq2[win_lo - 1..win_hi];
    let const_offset = match window.find(cfg.const_region) {
        Some(local) => local + win_lo - 1,
        None => return Ok(Err(FailReason::ConstantSeq)),
    };

    // -----------------------------------------------------
    // 3. Reject when constant region is too skewed
    // -----------------------------------------------------
    if const_offset < cfg.barcode_len
        || const_offset + cfg.barcode_len + cfg.rbs_len > seq2.len()
    {
        return Ok(Err(FailReason::ConstantPos));
    }

    // -----------------------------------------------------
    // 4. Extract RBS
    // -----------------------------------------------------
    let rbs_start = const_offset + cfg.const_region.len();
    let rbs = &seq2[rbs_start..rbs_start + cfg.rbs_len];

    // -----------------------------------------------------
    // 5. Extract barcode 2
    // -----------------------------------------------------
    let barcode2 = &seq2[const_offset - cfg.barcode_len..const_offset];
    if !cfg.barcodes2.contains(&barcode2) {
        return Ok(Err(FailReason::Barcode2));
    }

    // -----------------------------------------------------
    // 6. Extract discriminator
    // -----------------------------------------------------
    let (disc_pos, flipped) =
        match (seq1.find(cfg.non_flipped), seq1.find(cfg.flipped)) {
            (Some(p), _) => (p, Flip::NonFlipped),
            (None, Some(p)) => (p, Flip::Flipped),
            _ => return Ok(Err(FailReason::DiscSeq)),
        };
    if disc_pos < cfg.disc_offset + cfg.barcode_len {
        return Ok(Err(FailReason::DiscPos));
    }

    // -----------------------------------------------------
    // 6. Extract barcode 1
    // -----------------------------------------------------
    let barcode1_start = disc_pos - cfg.disc_offset - cfg.barcode_len;
    let barcode1 = &seq1[barcode1_start..barcode1_start + cfg.barcode_len];
    if !cfg.barcodes1.contains(&barcode1) {
        return Ok(Err(FailReason::Barcode1));
    }

    // -----------------------------------------------------
    // 7. End
    // -----------------------------------------------------

    Ok(Ok((
        Sample {
            barcode1: barcode1.to_owned(),
            barcode2: barcode2.to_owned(),
        },
        rbs,
        flipped,
    )))
}

// =========================================================
// Main processing function
// =========================================================

pub fn process_fastq(
    path1: &str,
    path2: &str,
    sample_name: &str,
    output_dir: &str,
    chunk_size: usize,
    parquet_size: usize,
) {
    info!("Creating output directories if they do not exist");

    // -----------------------------------------------------
    // Create output directories
    // -----------------------------------------------------

    let dirs = match prepare_dirs(output_dir) {
        Ok(dir) => dir,
        Err(e) => {
            error!("Failed to create output directories: {}", e);
            panic!("Failed to create output directories");
        }
    };

    // -----------------------------------------------------
    // Configuration
    // -----------------------------------------------------

    let cfg = Config {
        barcodes1: &constants::BARCODES_1,
        barcodes2: &constants::BARCODES_2,
        const_region: constants::CONSTANT_REGION,
        window: constants::CONSTANT_REGION_WINDOW,
        rbs_len: constants::RBS_LEN,
        barcode_len: constants::BARCODE_LEN,
        max_n: constants::MAX_N_COUNT,
        non_flipped: constants::NON_FLIPPED_SEQ,
        flipped: constants::FLIPPED_SEQ,
        disc_offset: constants::DISCRIMINATOR_OFFSET,
    };

    // -----------------------------------------------------
    // Load FASTQ files
    // -----------------------------------------------------

    info!("Processing FASTQ files: {} and {}", path1, path2);

    let file1 = BufReader::new(File::open(path1).unwrap());
    let file2 = BufReader::new(File::open(path2).unwrap());

    let reader1 = fastq::Reader::new(MultiGzDecoder::new(file1));
    let reader2 = fastq::Reader::new(MultiGzDecoder::new(file2));

    let mut iter1 = reader1.records();
    let mut iter2 = reader2.records();

    // -----------------------------------------------------
    // Initialise counters
    // -----------------------------------------------------

    let mut i = 0;
    let mut n = 0;
    let counters = Arc::new(Counters::default());

    // -----------------------------------------------------
    // Process FASTQ files in chunks
    // -----------------------------------------------------
    loop {
        info!("Processing {}", n);

        let chunk1: Vec<_> = iter1.by_ref().take(chunk_size).collect();
        let chunk2: Vec<_> = iter2.by_ref().take(chunk_size).collect();

        if chunk1.is_empty() || chunk2.is_empty() {
            info!("No more records to process.");
            break;
        }

        let table: SampleTable = DashMap::new();

        chunk1
            .par_iter()
            .zip(chunk2.par_iter())
            .for_each(|(rec1, rec2)| {
                counters.inc_total();

                match classify_pair(
                    &cfg,
                    &rec1.as_ref().unwrap(),
                    &rec2.as_ref().unwrap(),
                ) {
                    Ok(Ok((sample, rbs, flipped))) => {
                        counters.inc_valid();

                        let inner = table
                            .entry(sample)
                            .or_insert_with(|| DashMap::new());

                        let cell = inner
                            .entry(rbs.to_owned())
                            .or_insert_with(|| Default::default());

                        let idx = match flipped {
                            Flip::NonFlipped => 0,
                            Flip::Flipped => 1,
                        };

                        cell[idx].fetch_add(1, Ordering::Relaxed);
                    }
                    Ok(Err(reason)) => counters.inc_fail(reason),
                    Err(func_err) => panic!("Problem in pairs: {}", func_err),
                }
            });

        n += chunk1.len();
        i += 1;

        // -----------------------------------------------------
        // Write results to Parquet
        // -----------------------------------------------------
        let df = match table_to_dataframe(&table) {
            Ok(df) => df,
            Err(e) => {
                // TODO: Handle zero padding more gracefully
                panic!("chunk {i:06}: table to DataFrame failed: {e}");
            }
        };

        let path = dirs.parquet.join(format!("chunk_{i:09}.parquet"));

        match write_parquet_chunk(&df, &path.to_string_lossy()) {
            Ok(_) => {
                info!("Wrote {} ({} rows)", path.display(), df.height());
            }
            Err(err) => panic!("chunk {i:06}: failed to write parquet: {err}"),
        }
    }

    // -----------------------------------------------------
    // Save QC results

    info!("Write QC parquet file");
    let qc = counters.to_dataframe().unwrap();

    match write_qc_parquet(&qc, &dirs.qc, "sample_name") {
        Ok(_) => info!("Wrote QC parquet file"),
        Err(err) => panic!("Couldn't write QC parquet file: {err}"),
    }

    // -----------------------------------------------------
    // Write final results to Parquet

    info!("Merging Parquet files...");
    let counts = concat_parquet_dir(&dirs.parquet);

    match write_partitioned_parquet(
        &counts,
        &dirs.counts,
        sample_name,
        Some(parquet_size),
    ) {
        Ok(_) => info!("Wrote counts parquet files"),
        Err(err) => panic!("Couldn't write counts parquet files: {err}"),
    }

    info!("Processing complete.");
}
