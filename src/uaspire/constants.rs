// Constants for the uASPIre method
//
// Quality thresholds
pub const MAX_N_COUNT: usize = 6;

// Barcodes
pub const BARCODES_1: [&str; 6] =
    ["ATCACG", "CGATGT", "CTTGTA", "GCCAAT", "ACAGTG", "ACTTGA"];
pub const BARCODES_2: [&str; 6] =
    ["ATCACG", "CGATGT", "CTTGTA", "GCCAAT", "ACAGTG", "ACTTGA"];
pub const BARCODE_LEN: usize = 6;

// Constant region
pub const CONSTANT_REGION: &str = "GAGCTCGCAT";
pub const CONSTANT_REGION_WINDOW: (usize, usize) = (7, 24);

// Flipped and non-flipped sequences
pub const NON_FLIPPED_SEQ: &str = "GGGTTTGTACCGTACAC";
pub const FLIPPED_SEQ: &str = "GCCCGGATGATCCTGAC";
pub const DISCRIMINATOR_OFFSET: usize = 6;

// RBSs
pub const RBS_LEN: usize = 17;
