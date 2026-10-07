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

use std::fmt::Debug;
use std::hash::{Hash, Hasher};

use crate::{cli::CollapseMode, utils::intern_chromosome};

use xxhash_rust::xxh3::xxh3_64;

/// Represents a genomic record parsed from BED12 format.
///
/// A `Record` contains a binary key for efficient comparison and hashing,
/// along with the original read name and line data as byte vectors.
#[derive(Clone, Debug)]
pub struct Record {
    /// Binary key containing parsed genomic coordinates and metadata
    pub key: BinKey,
    /// Read name as bytes
    pub read: Vec<u8>,
    /// Original line as bytes
    pub line: Vec<u8>,
    /// Record bounds (start, end)
    pub bounds: (u32, u32),
    /// Record strand
    pub strand: u8,
}

impl Record {
    /// Parses a BED12 format line into a `Record`.
    ///
    /// # Arguments
    ///
    /// * `line` - A string slice containing a tab-separated BED12 format line
    ///
    /// # Returns
    ///
    /// A `Result` containing the parsed `Record` on success, or a boxed error on failure.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The input line is empty
    /// - Required fields are missing or cannot be parsed
    /// - Numeric fields contain invalid values
    ///
    /// # BED12 Format
    ///
    /// Expected fields (tab-separated):
    /// 1. chrom - Chromosome name
    /// 2. start - Start position (0-based)
    /// 3. end - End position (exclusive)
    /// 4. name - Name of the item
    /// 5. score - Score (0-1000)
    /// 6. strand - Strand (+ or -)
    /// 7. thickStart - Start of thick drawing
    /// 8. thickEnd - End of thick drawing
    /// 9. itemRgb - RGB color value
    /// 10. blockCount - Number of blocks
    /// 11. blockSizes - Comma-separated block sizes
    /// 12. blockStarts - Comma-separated block starts
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// let bed_line = "chr1\t1000\t2000\titem1\t100\t+\t1100\t1900\t0\t2\t500,400\t0,600";
    /// let record = Record::parse(bed_line)?;
    /// assert_eq!(record.key.start, 1000);
    /// assert_eq!(record.key.end, 2000);
    /// ```
    pub fn parse(line: &str, mode: &CollapseMode) -> Result<Self, Box<dyn std::error::Error>> {
        if line.is_empty() {
            return Err("ERROR: Empty line".into());
        }

        let mut fields = line.split('\t');

        let (
            chrom,
            start,
            end,
            name,
            _,
            strand,
            thick_start,
            thick_end,
            _,
            _,
            block_sizes,
            block_starts,
        ) = (
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Could not pasrse chrom field from -> {}", line)),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Could not pasrse start field from -> {}", line))
                .parse::<u32>()
                .unwrap_or_else(|e| {
                    panic!("ERROR: Could not pasrse start field from -> {}. {e}", line)
                }),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Could not pasrse end field from -> {}", line))
                .parse::<u32>()
                .unwrap_or_else(|e| {
                    panic!("ERROR: Could not pasrse end field from -> {}. {e}", line)
                }),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Could not pasrse name field from -> {}", line)),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Could not pasrse score field from -> {}", line)),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Could not pasrse strand field from -> {}", line))
                .parse::<char>()
                .unwrap_or_else(|e| {
                    panic!("ERROR: Could not pasrse score field from -> {}. {e}", line)
                }) as u8,
            fields
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: Could not pasrse thick_start field from -> {}", line)
                })
                .parse::<u32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: Could not pasrse thick_start field from -> {}. {e}",
                        line
                    )
                }),
            fields
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: Could not pasrse thick_end field from -> {}", line)
                })
                .parse::<u32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: Could not pasrse thick_end field from -> {}. {e}",
                        line
                    )
                }),
            fields.next().unwrap_or_else(|| {
                panic!("ERROR: Could not pasrse item_rgb field from -> {}", line)
            }),
            fields.next().unwrap_or_else(|| {
                panic!("ERROR: Could not pasrse block_count field from -> {}", line)
            }),
            fields
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: Could not pasrse block_sizes field from -> {}", line)
                })
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| {
                    s.parse::<u32>().unwrap_or_else(|e| {
                        panic!(
                            "ERROR: Could not pasrse block_sizes field from -> {}. {e}",
                            line
                        )
                    })
                })
                .collect::<Vec<_>>(),
            fields
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: Could not pasrse block_starts field from -> {}",
                        line
                    )
                })
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| {
                    s.parse::<u32>().unwrap_or_else(|e| {
                        panic!(
                            "ERROR: Could not pasrse block_starts field from -> {}. {e}",
                            line
                        )
                    })
                })
                .collect::<Vec<_>>(),
        );

        let gaps = block_sizes
            .iter()
            .zip(block_starts.iter())
            .map(|(block_size, block_start)| {
                (start + block_start, start + block_size + block_start)
            })
            .collect::<Vec<(u32, u32)>>()
            .windows(2)
            .map(|w| (w[0].1 + 1, w[1].0 - 1))
            .collect::<Vec<(u32, u32)>>();

        Ok(Self {
            key: BinKey::from_parts(
                chrom,
                start,
                end,
                strand,
                thick_start,
                thick_end,
                gaps,
                mode,
            ),
            read: name.as_bytes().to_vec(),
            line: line.as_bytes().to_vec(),
            bounds: (start, end),
            strand,
        })
    }
}

/// A binary key for efficient genomic coordinate comparison and hashing.
///
/// `BinKey` represents genomic intervals with chromosome, strand, coordinates,
/// CDS regions, and gaps (introns). It uses a cached fingerprint for fast
/// equality comparisons and hashing operations.
#[derive(Clone, Debug)]
pub struct BinKey {
    /// Interned chromosome identifier for memory efficiency
    pub chrom_id: u32,
    /// Strand as byte: b'+' for forward, b'-' for reverse
    pub strand: u8,
    /// Transcript start position (0-based, half-open interval)
    pub start: u32,
    /// Transcript end position (exclusive)
    pub end: u32,
    /// CDS thick start position (BED12 semantics)
    pub thick_start: u32,
    /// CDS thick end position (BED12 semantics)
    pub thick_end: u32,
    /// Gaps represented as vector of (intron_start, intron_end) tuples
    pub gaps: Vec<(u32, u32)>,
    /// Cached XXH3-64 hash fingerprint for fast comparisons
    pub fingerprint: u64,
}

impl BinKey {
    /// Creates a new `BinKey` from genomic coordinate components.
    ///
    /// # Arguments
    ///
    /// * `chrom` - Chromosome name (will be interned to ID)
    /// * `start` - Start position (0-based, inclusive)
    /// * `end` - End position (exclusive)
    /// * `strand` - Strand as byte (b'+' or b'-')
    /// * `thick_start` - CDS start position
    /// * `thick_end` - CDS end position
    /// * `gaps` - Vector of (start, end) tuples representing introns
    ///
    /// # Returns
    ///
    /// A new `BinKey` with computed fingerprint.
    ///
    /// # Notes
    ///
    /// - Invalid gaps (start >= end or outside transcript bounds) are filtered out
    /// - Gaps are clamped to the [start, end] interval
    /// - Chromosome names are interned for memory efficiency
    /// - A deterministic fingerprint is computed using XXH3-64
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// let gaps = vec![(1100, 1200), (1800, 1900)];
    /// let key = BinKey::from_parts("chr1", 1000, 2000, b'+', 1050, 1950, gaps);
    /// assert_eq!(key.start, 1000);
    /// assert_eq!(key.end, 2000);
    /// assert_eq!(key.gaps.len(), 2);
    /// ```
    #[allow(clippy::too_many_arguments)]
    pub fn from_parts(
        chrom: &str,
        mut start: u32,
        mut end: u32,
        strand: u8,
        mut thick_start: u32,
        mut thick_end: u32,
        mut gaps: Vec<(u32, u32)>,
        mode: &CollapseMode,
    ) -> Self {
        match mode {
            CollapseMode::Exon => {
                // INFO: ignores CDS and only performs exon coord matching
                //
                // read1: xxxxXXXX---XXXX----XXXXxxxx
                //          ^^||||           |^^^
                // read2: xxXXXXXX---XXXX----Xxxxxxxx
                //
                // Here, indenpendently of how thick blocks are definde, exon
                // coords will be used for matching. Note that gaps (introns) are
                // also taking into account.

                thick_start = 0;
                thick_end = 0;
            }
            CollapseMode::Coding => {
                // INFO: ignores transcription start and end, performs CDS coord matching
                //
                // read1: xxx---xxxXXXX---XXXX---XXXXx
                //                |^^^^
                // read2:         xXXXX---XXXX---XXXXx
                //
                // its negative variant:
                //
                // read1: xxxxXXXX---XXXX----XXXXxxxx
                //          ^^||||           ||||
                // read2: xxXXXXXX---XXXX----XXXXxxxx <- even gaps match, CDS does not
                //
                // Here, independently of exonic start and end, CDS and gaps will be
                // used for matching. Note that UTR structure will be completely ignored
                // in this category. Also, note that gaps are reduced within CDS bounds.

                start = thick_start;
                end = thick_end;

                // INFO: clamp gaps into [start, end]
                // INFO: also drop any intron that is invalid (start >= end or outside transcript)
                gaps.retain(|&(s, e)| s < e && s >= start && e <= end);
            }
            CollapseMode::GappedCoding => {
                // INFO: ignores exonic start and end but preserves all gaps
                //
                // read1: xxxxxxx--xxXXXX---XXX---XXXxx--xxx
                //            ||||||||||||||||||||||||||||||
                // read2:     xxx--xxXXXX---XXX---XXXxx--xxx
                // read3:      xx--xxXXXX---XXX---XXXxx--xxx
                // read4:       x--xxXXXX---XXX---XXXxx--xxx
                //
                // but excludes the following cases:
                //
                // read1: xx--xxxxXXXX---XXXX----XXXXxxxx
                //          **  ^^||||***||||****|^^^
                // read2: xx--xxXXXXXX---XXXX----Xxxxxxxx <- gaps match but not CDS
                //
                // read1: xxxxxxXXXXX----XXXX---XXXXXxxxx
                //          ^^^||||||||||||||||||||||||||
                // read2: xx---xXXXXX----XXXX---XXXXXxxxx <- CDS match but not gaps

                start = thick_start;
                end = thick_end;
            }
            CollapseMode::Transcript => {} // INFO: do nothing, all info is already in
        }

        // INFO: static indexing of chromosomes -> {chr: 1, chr2: 2, ...}
        let chrom_id = intern_chromosome(chrom);

        let mut key = BinKey {
            chrom_id,
            strand,
            start,
            end,
            thick_start,
            thick_end,
            gaps,
            fingerprint: 0,
        };

        // INFO: deterministic fingerprint and bytes
        let bytes = key.to_bytes_canonical();
        key.fingerprint = xxh3_64(&bytes);

        key
    }

    /// Converts the `BinKey` to a canonical byte representation.
    ///
    /// # Returns
    ///
    /// A `Vec<u8>` containing the canonical byte representation of the key.
    ///
    /// # Format
    ///
    /// The byte layout is:
    /// - chrom_id: 4 bytes (little-endian u32)
    /// - strand: 1 byte
    /// - start: 4 bytes (little-endian u32)
    /// - end: 4 bytes (little-endian u32)
    /// - thick_start: 4 bytes (little-endian u32)
    /// - thick_end: 4 bytes (little-endian u32)
    /// - n_gaps: 2 bytes (little-endian u16)
    /// - For each gap:
    ///   - gap_start: 4 bytes (little-endian u32)
    ///   - gap_length: 4 bytes (little-endian u32)
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// let key = BinKey::from_parts("chr1", 1000, 2000, b'+', 1100, 1900, vec![]);
    /// let bytes = key.to_bytes_canonical();
    /// assert!(bytes.len() >= 23); // Minimum size without gaps
    /// ```
    pub fn to_bytes_canonical(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(64 + self.gaps.len() * 8);

        bytes.extend_from_slice(&self.chrom_id.to_le_bytes());
        bytes.push(self.strand);
        bytes.extend_from_slice(&self.start.to_le_bytes());
        bytes.extend_from_slice(&self.end.to_le_bytes());
        bytes.extend_from_slice(&self.thick_start.to_le_bytes());
        bytes.extend_from_slice(&self.thick_end.to_le_bytes());

        let n_gaps = self.gaps.len() as u16;
        bytes.extend_from_slice(&n_gaps.to_le_bytes());

        for (s, e) in &self.gaps {
            bytes.extend_from_slice(&s.to_le_bytes());

            let len = e - s;
            bytes.extend_from_slice(&len.to_le_bytes());
        }

        bytes
    }

    /// Reconstructs a `BinKey` from its canonical byte representation.
    ///
    /// # Arguments
    ///
    /// * `bytes` - Byte slice containing the canonical representation
    ///
    /// # Returns
    ///
    /// A `BinKey` reconstructed from the byte data with computed fingerprint.
    ///
    /// # Panics
    ///
    /// Panics if the byte slice is malformed or too short for the expected format.
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// let original = BinKey::from_parts("chr1", 1000, 2000, b'+', 1100, 1900, vec![]);
    /// let bytes = original.to_bytes_canonical();
    /// let reconstructed = BinKey::from_bytes_canonical(&bytes);
    /// assert_eq!(original, reconstructed);
    /// ```
    pub fn from_bytes_canonical(bytes: &[u8]) -> Self {
        let buf = bytes;
        let mut bytes = bytes;

        let chrom_id = u32::from_le_bytes(bytes[..4].try_into().unwrap());
        bytes = &bytes[4..];

        let strand = bytes[0];
        bytes = &bytes[1..];

        let start = u32::from_le_bytes(bytes[..4].try_into().unwrap());
        bytes = &bytes[4..];

        let end = u32::from_le_bytes(bytes[..4].try_into().unwrap());
        bytes = &bytes[4..];

        let thick_start = u32::from_le_bytes(bytes[..4].try_into().unwrap());
        bytes = &bytes[4..];

        let thick_end = u32::from_le_bytes(bytes[..4].try_into().unwrap());
        bytes = &bytes[4..];

        let n_gaps = u16::from_le_bytes(bytes[..2].try_into().unwrap());
        bytes = &bytes[2..];

        let mut gaps = Vec::with_capacity(n_gaps as usize);
        for _ in 0..n_gaps {
            let s = u32::from_le_bytes(bytes[..4].try_into().unwrap());
            bytes = &bytes[4..];

            let len = u32::from_le_bytes(bytes[..4].try_into().unwrap());
            bytes = &bytes[4..];

            gaps.push((s, s + len));
        }

        Self {
            chrom_id,
            strand,
            start,
            end,
            thick_start,
            thick_end,
            gaps,
            fingerprint: xxh3_64(buf),
        }
    }
}

impl PartialEq for BinKey {
    /// Compares two `BinKey` instances for equality.
    ///
    /// # Arguments
    ///
    /// * `other` - Another `BinKey` to compare against
    ///
    /// # Returns
    ///
    /// `true` if the keys represent the same genomic interval, `false` otherwise.
    ///
    /// # Implementation Notes
    ///
    /// Uses a two-stage comparison:
    /// 1. Fast path: Compare fingerprints (different fingerprints = not equal)
    /// 2. Slow path: Compare canonical byte representations (guards against hash collisions)
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// let key1 = BinKey::from_parts("chr1", 1000, 2000, b'+', 1100, 1900, vec![]);
    /// let key2 = BinKey::from_parts("chr1", 1000, 2000, b'+', 1100, 1900, vec![]);
    /// assert_eq!(key1, key2);
    /// ```
    fn eq(&self, other: &Self) -> bool {
        // INFO: fast path, different fingerprints -> not equal
        if self.fingerprint != other.fingerprint {
            return false;
        }

        // INFO: verify by comparing canonical bytes (defends against hash collisions)
        self.to_bytes_canonical() == other.to_bytes_canonical()
    }
}

impl Eq for BinKey {}

impl Hash for BinKey {
    /// Computes the hash value for the `BinKey`.
    ///
    /// # Arguments
    ///
    /// * `state` - The hasher state to write to
    ///
    /// # Implementation Notes
    ///
    /// Uses the precomputed fingerprint for hashing to ensure:
    /// - Consistent hashing with equality comparison
    /// - Fast hash computation (no need to rehash all data)
    /// - Hash stability across different hasher implementations
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// use std::collections::HashMap;
    ///
    /// let mut map = HashMap::new();
    /// let key = BinKey::from_parts("chr1", 1000, 2000, b'+', 1100, 1900, vec![]);
    /// map.insert(key, "some_value");
    /// ```
    fn hash<H: Hasher>(&self, state: &mut H) {
        // INFO: we hash the fingerprint (u64). This keeps hashing cheap and consistent with Eq
        state.write_u64(self.fingerprint);
    }
}

/// Represents a collection of genomic reads with associated metadata.
///
/// `Queue` groups related reads together with their count, a representative line,
/// and header information for batch processing or output formatting.
#[derive(Debug, Clone)]
pub struct Queue {
    /// Collection of read names as byte vectors
    pub reads: Vec<Vec<u8>>,
    /// Total count of reads in this queue
    pub count: u32,
    /// Representative line for this group of reads
    pub rep_line: Vec<u8>,
    /// Header information for this queue
    pub header: Vec<u8>,
    /// Queue bounds (start, end)
    pub bounds: (u32, u32),
}

impl std::fmt::Display for Queue {
    /// Formats the `Queue` for display purposes.
    ///
    /// # Arguments
    ///
    /// * `f` - The formatter to write to
    ///
    /// # Returns
    ///
    /// A `Result` indicating success or failure of the formatting operation.
    ///
    /// # Output Format
    ///
    /// The display format includes:
    /// - Header: UTF-8 decoded header bytes
    /// - Count: Number of reads
    /// - Reads: List of UTF-8 decoded read names
    /// - Rep line: UTF-8 decoded representative line
    ///
    /// # Panics
    ///
    /// Panics if any byte vectors cannot be converted to valid UTF-8 strings.
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// let queue = Queue {
    ///     reads: vec![b"read1".to_vec(), b"read2".to_vec()],
    ///     count: 2,
    ///     rep_line: b"chr1\t1000\t2000\t...".to_vec(),
    ///     header: b"track name=example".to_vec(),
    /// };
    /// println!("{}", queue);
    /// ```
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "Header: {:?}",
            std::str::from_utf8(&self.header)
                .unwrap_or_else(|e| panic!("ERROR: Could not convert header to utf8 -> {e}"))
        )?;
        writeln!(f, "Count: {:?}", self.count)?;
        writeln!(
            f,
            "Reads: {:?}",
            self.reads
                .iter()
                .map(|r| std::str::from_utf8(r)
                    .unwrap_or_else(|e| panic!("ERROR: Could not convert read to utf8 -> {e}")))
                .collect::<Vec<_>>()
        )?;
        writeln!(
            f,
            "Rep line: {:?}",
            std::str::from_utf8(&self.rep_line)
                .unwrap_or_else(|e| panic!("ERROR: Could not convert rep_line to utf8 -> {e}"))
        )?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const MODE: CollapseMode = CollapseMode::Transcript;

    #[test]
    fn test_binkey_equality() {
        let key = BinKey::from_parts(
            "chr1",
            50,
            500,
            b'+',
            150,
            250,
            vec![(100, 200), (300, 400)],
            &MODE,
        );

        assert_eq!(key.chrom_id, 1);
        assert_eq!(key.strand, b'+');
        assert_eq!(key.start, 50);
        assert_eq!(key.end, 500);
        assert_eq!(key.thick_start, 150);
        assert_eq!(key.thick_end, 250);
        assert_eq!(key.gaps, vec![(100, 200), (300, 400)]);

        dbg!(&key.fingerprint);

        let key2 = BinKey::from_parts(
            "chr1",
            50,
            500,
            b'+',
            150,
            250,
            vec![(100, 200), (300, 400)],
            &MODE,
        );

        assert_eq!(key.fingerprint, key2.fingerprint);
        assert_eq!(key.to_bytes_canonical(), key2.to_bytes_canonical());
        assert_eq!(key, key2);
    }

    #[test]
    fn test_binkey_to_bytes() {
        let key = BinKey::from_parts(
            "chr1",
            100,
            200,
            b'+',
            150,
            250,
            vec![(100, 200), (300, 400)],
            &MODE,
        );

        let bytes = key.to_bytes_canonical();
        assert_eq!(bytes.len(), 23 + 2 * 8); // fixed header + 8 bytes per gap
    }

    #[test]
    fn test_binkey_from_bytes() {
        let key = BinKey::from_parts(
            "chr1",
            50,
            500,
            b'+',
            150,
            250,
            vec![(100, 200), (300, 400)],
            &MODE,
        );

        let bytes = key.to_bytes_canonical();
        let key2 = BinKey::from_bytes_canonical(&bytes);

        assert_eq!(key, key2);
    }
}
