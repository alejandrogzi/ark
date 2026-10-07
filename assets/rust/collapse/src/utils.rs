// Copyright (c) 2026 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Core module for collapsing BED files with deduplication and indexing
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main functions for efficiently collapsing BED files
//! by identifying and deduplicating identical genomic intervals. The module provides
//! flexible output modes including collapsed BED files with queue annotations or
//! separate binary index files for space-efficient storage, extending the original
//! BED file with additional columns for read name and queue information and fast
//! lookups fo read names.
//!
//! In short, every BED entry is fingerprinted using byte-level hashing accounting
//! for specific columns from the original BED file to ensure that only truly identical
//! rows are grouped together. The deduplicated entries are held in memory alongside
//! their corresponding read identifier queues [maintaining original order for reconstruction].

use std::cmp::{max, min};
use std::collections::HashMap;
use std::collections::{BTreeMap, VecDeque};
use std::fs::File;
use std::io::{self, prelude::*, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::RwLock;
use std::{fmt::Debug, io::BufReader};

use dashmap::DashMap;
use log;
use once_cell::sync::Lazy;
use rayon::prelude::*;

use crate::cli::CollapseMode;
use crate::{
    cli::RunMode,
    record::{BinKey, Queue, Record},
};

pub const MAGIC: &[u8; 5] = b"C0IDX";
pub const IDX_VERSION: u8 = 1;

static CHROM_INTERNER: Lazy<RwLock<HashMap<String, u32>>> =
    Lazy::new(|| RwLock::new(HashMap::with_capacity(1024)));
static NEXT_CHROM_ID: AtomicU32 = AtomicU32::new(1);

/// Unpacks and processes multiple BED12 files in parallel.
///
/// This function reads multiple BED12 files concurrently, parses their contents,
/// and groups records by their genomic coordinates into queues for efficient processing.
///
/// # Arguments
///
/// * `files` - A vector of file paths that implement `AsRef<Path> + Debug + Sync + Send`
///
/// # Returns
///
/// A `Result` containing a `HashMap<BinKey, Queue>` where:
/// - Keys are `BinKey` instances representing unique genomic intervals
/// - Values are `Queue` instances containing grouped reads for each interval
///
/// # Errors
///
/// Returns an error if:
/// - Any file cannot be read
/// - File contents cannot be parsed as valid BED12 format
/// - I/O operations fail during file reading
///
/// # Example
///
/// ```rust, ignore
/// let files = vec!["sample1.bed", "sample2.bed", "sample3.bed"];
/// let tracks = unpack(files)?;
/// println!("Processed {} unique genomic intervals", tracks.len());
/// ```
pub fn unpack<P: AsRef<Path> + Debug + Sync + Send>(
    files: Vec<P>,
    mode: CollapseMode,
    max_five_utr_length: u32,
    max_three_utr_length: u32,
) -> Result<DashMap<BinKey, Vec<Queue>>, Box<dyn std::error::Error>> {
    let contents = par_reader(files)?;
    let tracks = par_parse_tracks(&contents, mode, max_five_utr_length, max_three_utr_length)?;

    Ok(tracks)
}

/// Parses BED12 content in parallel and groups records by genomic coordinates.
///
/// This function processes multi-line BED12 content using parallel processing,
/// filtering out comments and invalid lines, then grouping valid records by
/// their `BinKey` (genomic coordinates).
///
/// # Arguments
///
/// * `contents` - A string slice containing BED12 formatted data
///
/// # Returns
///
/// A `Result` containing a `HashMap<BinKey, Queue>` with grouped genomic records.
///
/// # Errors
///
/// Returns an error if the parsing operations fail catastrophically.
///
/// # Processing Pipeline
///
/// 1. Split content into parallel lines
/// 2. Filter out comment lines (starting with '#')
/// 3. Parse each line into a `Record` (invalid lines are silently dropped)
/// 4. Group records by `BinKey` using parallel fold/reduce
/// 5. Merge duplicate keys by combining their queues
///
/// # Example
///
/// ```rust, ignore
/// let bed_content = "chr1\t1000\t2000\tread1\t100\t+\t1100\t1900\t0\t1\t1000\t0\n";
/// let tracks = par_parse_track
/// ```
fn par_parse_tracks(
    contents: &str,
    mode: CollapseMode,
    max_five_utr_length: u32,
    max_three_utr_length: u32,
) -> Result<DashMap<BinKey, Vec<Queue>>, Box<dyn std::error::Error>> {
    log::info!("INFO: Parsing BED12 content");

    let mut inner_tracks = contents
        .par_lines()
        .filter(|row| !row.starts_with("#"))
        .filter_map(|line| Record::parse(line, &mode).ok())
        .fold(HashMap::new, |mut acc, record| {
            // INFO: if record not in tracks, create a new queue
            acc.entry(record.key.clone())
                .and_modify(|records: &mut Vec<Record>| {
                    records.push(record.clone());
                })
                .or_insert(vec![record]);

            acc
        })
        .reduce(HashMap::new, |mut left, right| {
            // INFO: if record not in tracks, create a new queue
            for (key, right_queue) in right {
                left.entry(key)
                    .and_modify(|left_queue| {
                        left_queue.extend(right_queue.clone());
                    })
                    .or_insert(right_queue);
            }

            left
        });

    // INFO: 1) sort by start/end in descending order, 2) decide on records
    let tracks = DashMap::with_capacity(inner_tracks.len());
    inner_tracks.par_iter_mut().for_each(|(bin_key, records)| {
        records.par_sort_unstable_by(|a, b| {
            let (a_start, a_end) = a.bounds;
            let (b_start, b_end) = b.bounds;

            a_start.cmp(&b_start).then(b_end.cmp(&a_end))
        });

        let queues = collapse_to_queues(records, max_five_utr_length, max_three_utr_length);
        tracks.insert(bin_key.clone(), queues);
    });

    Ok(tracks)
}

/// Union-Find (DSU) data structure.
///
/// This struct represents a Union-Find data structure for efficiently
/// finding connected components in a graph. It maintains a parent array
/// and a size array to track the connected components.
pub(crate) struct Dsu {
    parent: Vec<usize>,
    size: Vec<usize>,
}

impl Dsu {
    /// Creates a new `Dsu` instance.
    ///
    /// # Arguments
    ///
    /// * `n` - The number of elements in the data structure
    ///
    /// # Returns
    ///
    /// A new `Dsu` instance with the specified number of elements.
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// let dsu = Dsu::new(10);
    /// ```
    pub(crate) fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            size: vec![1; n],
        }
    }

    /// Finds the parent of an element.
    ///
    /// # Arguments
    ///
    /// * `x` - The element to find the parent of
    ///
    /// # Returns
    ///
    /// The parent of the element.
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// let mut dsu = Dsu::new(10);
    /// dsu.find(5); // Returns 0
    /// ```
    pub(crate) fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    /// Performs union operation on two elements.
    ///
    /// # Arguments
    ///
    /// * `a` - The first element
    /// * `b` - The second element
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// let mut dsu = Dsu::new(10);
    /// dsu.union(5, 7); // Union of 5 and 7 is 5
    /// ```
    pub(crate) fn union(&mut self, a: usize, b: usize) {
        let mut ra = self.find(a);
        let mut rb = self.find(b);
        if ra == rb {
            return;
        }
        if self.size[ra] < self.size[rb] {
            std::mem::swap(&mut ra, &mut rb);
        }
        self.parent[rb] = ra;
        self.size[ra] += self.size[rb];
    }
}

/// Converts strand byte to boolean.
///
/// This function converts a strand byte to a boolean value.
///
/// # Arguments
///
/// * `strand` - A byte representing the strand
///
/// # Returns
///
/// A boolean value indicating the strand.
///
/// # Example
///
/// ```rust, ignore
/// let strand = b'+';
/// let is_plus = is_plus(strand);
/// assert_eq!(is_plus, true);
/// ```
#[inline]
fn is_plus(strand: u8) -> bool {
    strand == b'+'
}

/// Returns the 5' end of a `Record`.
///
/// This function returns the 5' end of a `Record` based on its strand.
///
/// # Arguments
///
/// * `rec` - A reference to a `Record` instance
///
/// # Returns
///
/// A `u32` representing the 5' end of the `Record`.
///
/// # Example
///
/// ```rust, ignore
/// let record = Record::parse("chr1\t1000\t2000\tread1\t100\t+\t1100\t1900\t0\t1\t1000\t0\n").unwrap();
/// let five = five_of(&record);
/// assert_eq!(five, 1000);
/// ```
#[inline]
fn five_of(rec: &Record) -> u32 {
    if is_plus(rec.strand) {
        rec.bounds.0
    } else {
        rec.bounds.1
    }
}

/// Returns the 3' end of a `Record`.
///
/// This function returns the 3' end of a `Record` based on its strand.
///
/// # Arguments
///
/// * `rec` - A reference to a `Record` instance
///
/// # Returns
///
/// A `u32` representing the 3' end of the `Record`.
///
/// # Example
///
/// ```rust, ignore
/// let record = Record::parse("chr1\t1000\t2000\tread1\t100\t+\t1100\t1900\t0\t1\t1000\t0\n").unwrap();
/// let three = three_of(&record);
/// assert_eq!(three, 1900);
/// ```
#[inline]
fn three_of(rec: &Record) -> u32 {
    if is_plus(rec.strand) {
        rec.bounds.1
    } else {
        rec.bounds.0
    }
}

/// Collapse a homogeneous set of `Record`s (same CDS/introns; UTRs vary)
/// into `Queue`s, using the rule:
///
///   |Δ5′| ≤ T5  AND  |Δ3′| ≤ T3
///
/// where T5 = max_five_utr_len (∞ if 0), T3 = max_three_utr_len (∞ if 0).
///
///
/// Arguments:
///
/// * `records` - A slice of `Record` instances
/// * `max_five_utr_len` - The maximum length of the 5'UTR
/// * `max_three_utr_len` - The maximum length of the 3'UTR
///
/// Returns:
///
/// A `Vec<Queue>` containing the collapsed `Queue` instances
///
/// # Example
///
/// ```rust, ignore
/// let records = vec![
///     Record::parse("chr1\t1000\t2000\tread1\t100\t+\t1100\t1900\t0\t1\t1000\t0\n").unwrap(),
///     Record::parse("chr1\t1000\t2000\tread2\t100\t+\t1100\t1900\t0\t1\t1000\t0\n").unwrap(),
///     Record::parse("chr1\t1000\t2000\tread3\t100\t+\t1100\t1900\t0\t1\t1000\t0\n").unwrap(),
/// ];
/// let queues = collapse_to_queues(&records, 100, 100);
/// assert_eq!(queues.len(), 1);
/// ```
pub fn collapse_to_queues(
    records: &[Record],
    max_five_utr_len: u32,
    max_three_utr_len: u32,
) -> Vec<Queue> {
    let n = records.len();
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![Queue {
            reads: vec![],
            count: 1,
            rep_line: records[0].line.clone(),
            header: records[0].read.clone(),
            bounds: records[0].bounds,
        }];
    }

    // INFO: strand normalization
    let plus = is_plus(records[0].strand);
    let mut five = Vec::with_capacity(n);
    let mut three = Vec::with_capacity(n);
    for r in records {
        five.push(five_of(r));
        three.push(three_of(r));
    }

    // INFO: order by five asc, then three asc (deterministic)
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by(|&i, &j| five[i].cmp(&five[j]).then(three[i].cmp(&three[j])));

    // INFO: effective thresholds (0 => "no limit")
    let t5 = if max_five_utr_len == 0 {
        u32::MAX
    } else {
        max_five_utr_len
    };
    let t3 = if max_three_utr_len == 0 {
        u32::MAX
    } else {
        max_three_utr_len
    };

    let mut dsu = Dsu::new(n);

    // INFO: sliding window over 'five' to enforce |Δ5′| ≤ T5
    let mut window: VecDeque<usize> = VecDeque::new();
    let mut alive = vec![false; n]; // tracks membership in window

    // INFO: For 3′ constraint, keep previous window members indexed by their 'three'
    // INFO: We don't physically remove; we check 'alive' on use.
    let mut by_three: BTreeMap<u32, Vec<usize>> = BTreeMap::new();

    let mut left = 0usize; // pointer into 'order'

    for &i in &order {
        let fi = five[i];
        let ti = three[i];

        // INFO: Evict from window anything with five < fi - t5
        while left < order.len() {
            let j = order[left];
            if five[j] >= fi.saturating_sub(t5) {
                break;
            }
            // INFO: Mark dead; leave in by_three (we'll skip by 'alive' check)
            alive[j] = false;
            left += 1;
            if let Some(front) = window.front() {
                if *front == j {
                    window.pop_front();
                }
            }
        }

        // INFO: 3′ candidates: those in [ti - t3, ti + t3]
        let lo = ti.saturating_sub(t3);
        let hi = ti.saturating_add(t3);
        for (_k, idxs) in by_three.range(lo..=hi) {
            for &j in idxs {
                if !alive[j] {
                    continue;
                } // INFO: ensure still in 5′ window
                  // 5′ constraint is guaranteed by construction of the window,
                  // but we can assert defensively:
                debug_assert!(fi.abs_diff(five[j]) <= t5);
                dsu.union(i, j);
            }
        }

        // INFO: Insert current index into window and by_three
        window.push_back(i);
        alive[i] = true;
        by_three.entry(ti).or_default().push(i);
    }

    // INFO: aggregate DSU components into queues
    #[derive(Clone)]
    struct GroupAgg {
        // INFO: extremes in normalized space
        min_five: u32,
        max_five: u32,
        min_three: u32,
        max_three: u32,
        rep: usize,
        reads: Vec<Vec<u8>>,
        count: u32,
        seen: bool,
    }

    let mut agg: Vec<GroupAgg> = vec![
        GroupAgg {
            min_five: u32::MAX,
            max_five: 0,
            min_three: u32::MAX,
            max_three: 0,
            rep: 0,
            reads: Vec::new(),
            count: 0,
            seen: false
        };
        n
    ];

    for i in 0..n {
        let r = dsu.find(i);
        let gi = &mut agg[r];
        if !gi.seen {
            gi.seen = true;
            gi.rep = i;
            gi.min_five = five[i];
            gi.max_five = five[i];
            gi.min_three = three[i];
            gi.max_three = three[i];
        } else {
            gi.min_five = min(gi.min_five, five[i]);
            gi.max_five = max(gi.max_five, five[i]);
            gi.min_three = min(gi.min_three, three[i]);
            gi.max_three = max(gi.max_three, three[i]);
        }

        gi.reads.push(records[i].read.clone());
        gi.count += 1;
    }

    let mut out = Vec::new();
    out.reserve(n);
    for mut g in agg.into_iter().filter(|g| g.seen) {
        let rep = &records[g.rep];

        // Map component extremes back to genomic (start, end)
        let new_five = if plus { g.min_five } else { g.max_five }; // most 5′-extreme
        let new_three = if plus { g.max_three } else { g.min_three }; // most 3′-extreme
        let new_bounds = if plus {
            (new_five, new_three)
        } else {
            (new_three, new_five)
        };

        // INFO: exclude header from queue reads
        g.reads.remove(0);

        out.push(Queue {
            reads: g.reads,
            count: g.count,
            rep_line: rep.line.clone(),
            header: rep.read.clone(),
            bounds: new_bounds,
        });
    }

    // Optional: keep your preferred order (start asc, end desc)
    out.sort_unstable_by(|a, b| {
        let (a_start, a_end) = a.bounds;
        let (b_start, b_end) = b.bounds;
        a_start.cmp(&b_start).then(b_end.cmp(&a_end))
    });

    out
}

/// Reads the entire contents of a file into a string.
///
/// # Arguments
///
/// * `file` - A file path that implements `AsRef<Path> + Debug`
///
/// # Returns
///
/// A `Result` containing the file contents as a `String`, or an error if reading fails.
///
/// # Errors
///
/// Returns an error if:
/// - The file cannot be opened (doesn't exist, permission denied, etc.)
/// - I/O errors occur during reading
/// - The file contains invalid UTF-8 sequences
///
/// # Example
///
/// ```rust, ignore
/// let content = reader("data.bed")?;
/// println!("File contains {} characters", content.len());
/// ```
fn reader<P: AsRef<Path> + Debug>(file: P) -> Result<String, Box<dyn std::error::Error>> {
    let mut file = File::open(file)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}

/// Reads multiple files in parallel and concatenates their contents.
///
/// This function leverages Rayon's parallel iterator to read multiple files
/// concurrently, then concatenates all file contents into a single string.
///
/// # Arguments
///
/// * `files` - A vector of file paths that implement `AsRef<Path> + Debug + Sync + Send`
///
/// # Returns
///
/// A `Result` containing the concatenated contents of all files as a single `String`.
///
/// # Errors
///
/// Returns an error if any file cannot be read successfully.
///
/// # Panics
///
/// Panics if any individual file reading operation fails, with a descriptive error message
/// indicating which file caused the failure.
///
/// # Performance Notes
///
/// - Files are read concurrently using Rayon's parallel iterators
/// - Memory usage scales with the total size of all input files
/// - No attempt is made to add separators between file contents
///
/// # Example
///
/// ```rust, ignore
/// let files = vec!["file1.bed", "file2.bed", "file3.bed"];
/// let combined_content = par_reader(files)?;
/// println!("Combined {} files into {} characters", files.len(), combined_content.len());
/// ```
fn par_reader<P: AsRef<Path> + Debug + Sync + Send>(
    files: Vec<P>,
) -> Result<String, Box<dyn std::error::Error>> {
    let contents: Vec<String> = files
        .par_iter()
        .map(|path| reader(path).unwrap_or_else(|e| panic!("Error reading file: {:?}", e)))
        .collect();

    Ok(contents.concat())
}

/// Interns a chromosome name to a small u32 identifier in a thread-safe manner.
///
/// This function provides a thread-safe string interning system for chromosome names,
/// converting them to compact u32 identifiers for memory efficiency. Uses a two-phase
/// locking strategy for optimal read performance.
///
/// # Arguments
///
/// * `chrom` - A string slice containing the chromosome name (e.g., "chr1", "chrX", "scaffold_1")
///
/// # Returns
///
/// A `u32` identifier that uniquely represents the chromosome name.
///
/// # Thread Safety
///
/// This function is fully thread-safe and uses the following locking strategy:
/// 1. **Fast path**: Attempts to read with a shared read lock
/// 2. **Slow path**: Acquires exclusive write lock only if chromosome is not found
/// 3. **Double-checked locking**: Rechecks existence after acquiring write lock
///
/// # Panics
///
/// Panics if the internal `RwLock` is poisoned due to a panic in another thread.
///
/// # Example
///
/// ```rust, ignore
/// let chr1_id = intern_chromosome("chr1");
/// let chr2_id = intern_chromosome("chr2");
/// let chr1_again = intern_chromosome("chr1");
///
/// assert_eq!(chr1_id, chr1_again); // Same chromosome returns same ID
/// assert_ne!(chr1_id, chr2_id);    // Different chromosomes get different IDs
/// ```
///
/// # Global State
///
/// This function modifies global static variables:
/// - `CHROM_INTERNER`: Thread-safe map of chromosome names to IDs
/// - `NEXT_CHROM_ID`: Atomic counter for assigning new IDs
pub fn intern_chromosome(chrom: &str) -> u32 {
    // INFO: Fast path (read lock)
    if let Ok(read) = CHROM_INTERNER.read() {
        if let Some(&id) = read.get(chrom) {
            return id;
        }
    }

    // INFO: Slow path (write lock)
    let mut write = CHROM_INTERNER
        .write()
        .unwrap_or_else(|e| panic!("ERROR: poisoned lock: {e} -> {chrom} - {CHROM_INTERNER:?}"));

    if let Some(&id) = write.get(chrom) {
        return id;
    }

    let id = NEXT_CHROM_ID.fetch_add(1, Ordering::Relaxed);
    write.insert(chrom.to_string(), id);
    id
}

/// Writes collapsed genomic tracks to a file in different output modes.
///
/// This function outputs the processed genomic tracks to a file, with different
/// formatting options based on the specified run mode.
///
/// # Arguments
///
/// * `tracks` - A `HashMap<BinKey, Queue>` containing the processed genomic data
/// * `mode` - A `RunMode` enum specifying the output format:
///   - `RunMode::Extend`: Includes read names as an additional tab-separated column
///   - `RunMode::Index`: Outputs only the representative line without extra columns
/// * `output` - A file path that implements `AsRef<Path> + Debug` where output will be written
///
/// # Output Formats
///
/// ## RunMode::Extend
/// ```text
/// chr1\t1000\t2000\titem1\t100\t+\t1100\t1900\t0\t1\t1000\t0\tread1,read2,read3
/// ```
///
/// ## RunMode::Index
/// ```text
/// chr1\t1000\t2000\titem1\t100\t+\t1100\t1900\t0\t1\t1000\t0
/// ```
///
/// # Panics
///
/// Panics if:
/// - The output file cannot be created
/// - Write operations fail
/// - Byte data cannot be converted to valid UTF-8
///
/// # Example
///
/// ```rust, ignore
/// let tracks = unpack(vec!["input.bed"])?;
/// __write_collapsed(tracks, RunMode::Extend, "output.bed");
/// ```
pub fn __write_collapsed<P: AsRef<Path> + Debug>(
    tracks: DashMap<BinKey, Vec<Queue>>,
    mode: RunMode,
    output: P,
) {
    log::info!("INFO: Writing collapsed file");
    let mut writer = BufWriter::new(
        File::create(&output)
            .unwrap_or_else(|e| panic!("ERROR: Could not create file {:?} -> {e}", output)),
    );

    let mut lines = 0;
    for (_, queues) in tracks {
        for queue in queues {
            lines += 1;
            match mode {
                RunMode::Extend => {
                    let mut line = std::str::from_utf8(&queue.rep_line)
                        .unwrap_or_else(|e| {
                            panic!("ERROR: Could not convert rep_line to utf8 -> {e}")
                        })
                        .to_string();

                    // INFO: appending queue.reads as extra column
                    if !queue.reads.is_empty() {
                        let tail = queue
                            .reads
                            .iter()
                            .map(|r| {
                                std::str::from_utf8(r).unwrap_or_else(|e| {
                                    panic!("ERROR: Could not convert read to utf8 -> {e}")
                                })
                            })
                            .collect::<Vec<_>>()
                            .join(",");

                        line = format!("{}\t{}", line, tail);
                    } else {
                        // INFO: if empty, append empty string to make format compatible!
                        line = format!("{}\tUNIQUE", line);
                    }

                    __write_line(&mut writer, &line);
                }
                RunMode::Index => {
                    // INFO: no extra column
                    let line = std::str::from_utf8(&queue.rep_line).unwrap_or_else(|e| {
                        panic!("ERROR: Could not convert rep_line to utf8 -> {e}")
                    });

                    __write_line(&mut writer, line);
                }
            }
        }
    }

    log::info!("INFO: Wrote {} lines to collapsed file", lines);
}

/// Writes a line to a file.
///
/// This function writes a line to a file using a buffered writer.
///
/// # Arguments
///
/// * `writer` - A reference to a `BufWriter` instance
/// * `line` - A string slice containing the line to write
///
/// # Panics
///
/// Panics if there is an error writing to the file.
fn __write_line(writer: &mut BufWriter<File>, line: &str) {
    writeln!(writer, "{}", line)
        .unwrap_or_else(|e| panic!("ERROR: Could not write to collapsed file -> {e}"));
}

/// Writes a binary index file containing genomic tracks data.
///
/// This function creates a custom binary index format that stores all genomic
/// track information in a compact, efficiently readable format suitable for
/// fast random access and reconstruction.
///
/// # Arguments
///
/// * `tracks` - A reference to `HashMap<BinKey, Queue>` containing the genomic data
/// * `output` - A file path that implements `AsRef<Path> + Debug` for the output index file
///
/// # Binary Format Structure
///
/// The index file uses the following binary layout:
/// ```text
/// Header:
///   - Magic bytes (5 bytes): File format identifier
///   - Version (1 byte): Index format version
///   - Number of records (8 bytes, little-endian u64)
///
/// For each record:
///   - Key length (4 bytes, little-endian u32)
///   - Key data (variable length): Canonical BinKey bytes
///   - Header length (4 bytes, little-endian u32)
///   - Header data (variable length bytes)
///   - Count (4 bytes, little-endian u32)
///   - Rep line length (4 bytes, little-endian u32)
///   - Rep line data (variable length bytes)
///   - Number of reads (4 bytes, little-endian u32)
///   - For each read:
///     - Read length (4 bytes, little-endian u32)
///     - Read data (variable length bytes)
/// ```
///
/// # Example
///
/// ```rust, ignore
/// let tracks = unpack(vec!["input.bed"])?;
/// __write_index(&tracks, "genome.idx");
///
/// // Later, read the index back
/// let reconstructed = read_index("genome.idx")?;
/// assert_eq!(tracks.len(), reconstructed.len());
/// ```
///
/// # See Also
///
/// - [`read_index`] for reading index files created by this function
/// - [`BinKey::to_bytes_canonical`] for the key serialization format
pub fn __write_index<P: AsRef<Path> + Debug>(tracks: &DashMap<BinKey, Vec<Queue>>, output: P) {
    log::info!("INFO: Writing index file");

    let mut writer = BufWriter::new(
        File::create(&output)
            .unwrap_or_else(|e| panic!("ERROR: Could not create file {:?} -> {e}", output)),
    );

    writer.write_all(MAGIC).unwrap_or_else(|e| {
        panic!(
            "ERROR: Could not write magic to index file {:?} -> {e}",
            output
        )
    });

    writer.write_all(&[IDX_VERSION]).unwrap_or_else(|e| {
        panic!(
            "ERROR: Could not write version to index file {:?} -> {e}",
            output
        )
    });

    let num_records = tracks.len() as u64;
    writer
        .write_all(&num_records.to_le_bytes())
        .unwrap_or_else(|e| {
            panic!(
                "ERROR: Could not write number of records to index file {:?} -> {e}",
                output
            )
        });

    for group in tracks.into_iter() {
        let key = group.key();
        let queues = group.value();

        // INFO: write key
        let key_bytes = key.to_bytes_canonical();
        let key_len = key_bytes.len() as u32;
        writer
            .write_all(&key_len.to_le_bytes())
            .unwrap_or_else(|e| {
                panic!(
                    "ERROR: Could not write key length to index file {:?} -> {e}",
                    output
                )
            });
        writer.write_all(&key_bytes).unwrap_or_else(|e| {
            panic!(
                "ERROR: Could not write key to index file {:?} -> {e}",
                output
            )
        });

        for queue in queues {
            // INFO: write header
            let header_len = queue.header.len() as u32;
            writer
                .write_all(&header_len.to_le_bytes())
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: Could not write header length to index file {:?} -> {e}",
                        output
                    )
                });
            writer.write_all(&queue.header).unwrap_or_else(|e| {
                panic!(
                    "ERROR: Could not write header to index file {:?} -> {e}",
                    output
                )
            });

            // INFO: write count
            writer
                .write_all(&queue.count.to_le_bytes())
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: Could not write count to index file {:?} -> {e}",
                        output
                    )
                });

            // INFO: write rep_line
            let rep_line_len = queue.rep_line.len() as u32;
            writer
                .write_all(&rep_line_len.to_le_bytes())
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: Could not write rep_line length to index file {:?} -> {e}",
                        output
                    )
                });
            writer.write_all(&queue.rep_line).unwrap_or_else(|e| {
                panic!(
                    "ERROR: Could not write rep_line to index file {:?} -> {e}",
                    output
                )
            });

            // INFO: write reads
            let n_reads = queue.reads.len() as u32;
            writer
                .write_all(&n_reads.to_le_bytes())
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: Could not write number of reads to index file {:?} -> {e}",
                        output
                    )
                });
            for read in &queue.reads {
                let rl = read.len() as u32;
                writer.write_all(&rl.to_le_bytes()).unwrap_or_else(|e| {
                    panic!(
                        "ERROR: Could not write read length to index file {:?} -> {e}",
                        output
                    )
                });
                writer.write_all(read).unwrap_or_else(|e| {
                    panic!(
                        "ERROR: Could not write read to index file {:?} -> {e}",
                        output
                    )
                });
            }
        }
    }

    log::info!("INFO: Wrote {} lines to index file", num_records);
}

/// Reads and reconstructs genomic tracks data from a binary index file.
///
/// This function reads a binary index file created by [`__write_index`] and
/// reconstructs the original `HashMap<BinKey, Queue>` data structure with
/// full fidelity to the original data.
///
/// # Arguments
///
/// * `input` - A file path that implements `AsRef<Path>` pointing to the index file
///
/// # Returns
///
/// An `io::Result<HashMap<BinKey, Queue>>` containing:
/// - `Ok(HashMap)`: Successfully reconstructed genomic tracks data
/// - `Err(io::Error)`: File I/O errors or format validation failures
///
/// # Errors
///
/// Returns `io::Error` for the following conditions:
///
/// ## File Access Errors
/// - File does not exist or cannot be opened
/// - Insufficient permissions to read the file
/// - I/O errors during reading operations
///
/// ## Format Validation Errors
/// - **Invalid Magic Bytes**: File is not a valid index format
/// - **Unsupported Version**: Index was created with an incompatible version
/// - **Truncated Data**: File ends unexpectedly during reading
/// - **Corrupted Length Fields**: Invalid length values that exceed remaining data
///
/// # Binary Format Compatibility
///
/// The function expects the exact binary format produced by [`__write_index`]:
/// - Magic bytes must match `MAGIC` constant exactly
/// - Version byte must match `IDX_VERSION` exactly
/// - All multi-byte values must be in little-endian format
/// - Length-prefixed data must be consistent
///
/// # Memory Usage
///
/// - Pre-allocates HashMap with capacity based on record count
/// - Memory usage scales with the total size of all genomic data
/// - Each `BinKey` is reconstructed from canonical bytes with fingerprint computation
///
/// # Example
///
/// ```rust, ignore
/// // Write an index file
/// let original_tracks = unpack(vec!["input.bed"])?;
/// __write_index(&original_tracks, "genome.idx");
///
/// // Read it back
/// match read_index("genome.idx") {
///     Ok(reconstructed_tracks) => {
///         assert_eq!(original_tracks.len(), reconstructed_tracks.len());
///         println!("Successfully loaded {} genomic intervals", reconstructed_tracks.len());
///     }
///     Err(e) => eprintln!("Failed to read index: {}", e),
/// }
/// ```
///
/// # Data Integrity
///
/// - BinKey fingerprints are recomputed during reconstruction
/// - All byte data is preserved exactly as originally stored
/// - HashMap iteration order may differ from original due to hash recomputation
///
/// # See Also
///
/// - [`__write_index`] for creating compatible index files
/// - [`BinKey::from_bytes_canonical`] for key reconstruction details
pub fn read_index<P: AsRef<Path>>(input: P) -> io::Result<HashMap<BinKey, Queue>> {
    let file = File::open(input)?;
    let mut reader = BufReader::new(file);

    // INFO: read magic and version
    let mut magic = [0u8; 5];
    reader.read_exact(&mut magic)?;
    if &magic != MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid index magic: expected {:?}, got {:?}", MAGIC, magic),
        ));
    }
    let mut version_buf = [0u8; 1];
    reader.read_exact(&mut version_buf)?;
    let version = version_buf[0];
    if version != IDX_VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Unsupported index version: {}", version),
        ));
    }

    // INFO: number of records
    let mut num_records_buf = [0u8; 8];
    reader.read_exact(&mut num_records_buf)?;
    let num_records = u64::from_le_bytes(num_records_buf);

    let mut map = HashMap::with_capacity(num_records as usize);

    for _ in 0..num_records {
        // INFO: key
        let mut u32_buf = [0u8; 4];

        reader.read_exact(&mut u32_buf)?;
        let key_len = u32::from_le_bytes(u32_buf) as usize;
        let mut key_bytes = vec![0u8; key_len];
        reader.read_exact(&mut key_bytes)?;

        // INFO: header
        reader.read_exact(&mut u32_buf)?;
        let header_len = u32::from_le_bytes(u32_buf) as usize;
        let mut header = vec![0u8; header_len];
        reader.read_exact(&mut header)?;

        // INFO: count
        reader.read_exact(&mut u32_buf)?;
        let count = u32::from_le_bytes(u32_buf);

        // INFO: rep_line
        reader.read_exact(&mut u32_buf)?;
        let rep_line_len = u32::from_le_bytes(u32_buf) as usize;
        let mut rep_line = vec![0u8; rep_line_len];
        reader.read_exact(&mut rep_line)?;

        // INFO: reads
        reader.read_exact(&mut u32_buf)?;
        let n_reads = u32::from_le_bytes(u32_buf) as usize;
        let mut reads: Vec<Vec<u8>> = Vec::with_capacity(n_reads);
        for _ in 0..n_reads {
            reader.read_exact(&mut u32_buf)?;
            let rl = u32::from_le_bytes(u32_buf) as usize;
            let mut r = vec![0u8; rl];
            reader.read_exact(&mut r)?;
            reads.push(r);
        }

        // INFO: reconstruct BinKey from canonical bytes
        let key = BinKey::from_bytes_canonical(&key_bytes);

        let queue = Queue {
            header,
            reads,
            count,
            rep_line,
            bounds: (0, 0), // INFO: placeholder bounds for now
        };

        map.insert(key, queue);
    }

    Ok(map)
}

/// Reads and swaps the index file for a `HashMap<Vec<u8>, Queue>`.
///
/// This function reads the index file and swaps the keys and values
/// in the `HashMap` to maintain the original order of the genomic intervals.
///
/// # Arguments
///
/// * `index` - A file path that implements `AsRef<Path>` pointing to the index file
///
/// # Returns
///
/// A `HashMap<Vec<u8>, Queue>` containing the swapped index data.
///
/// # Errors
///
/// Returns an error if the index file cannot be read.
///
/// # Example
///
/// ```rust, ignore
/// let original_tracks = unpack(vec!["input.bed"])?;
/// __write_index(&original_tracks, "genome.idx");
///
/// let reconstructed_tracks = read_and_swap_index("genome.idx");
/// assert_eq!(original_tracks.len(), reconstructed_tracks.len());
/// ```
pub fn read_and_swap_index<P: AsRef<Path> + Debug>(index: P) -> HashMap<Vec<u8>, Queue> {
    log::info!("INFO: Reading index file from {:?}", index.as_ref());
    let mut map = read_index(&index).unwrap_or_else(|e| {
        panic!(
            "ERROR: Could not read index file {:?} -> {e}",
            index.as_ref()
        )
    });

    log::info!("INFO: Swapping index file");
    let mut swapped = HashMap::new();
    for (_, queue) in map.drain() {
        swapped.insert(queue.header.clone(), queue);
    }

    log::info!("INFO: Swapped {} lines", swapped.len());
    swapped
}
