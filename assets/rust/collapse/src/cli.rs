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

use clap::{self, ArgAction, ArgGroup, Parser, ValueEnum};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[clap(name="collapse", author = env!("CARGO_PKG_AUTHORS"), version = env!("CARGO_PKG_VERSION"), about = "shrink your .beds", long_about = None)]
pub struct Args {
    #[clap(subcommand)]
    pub command: Command,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[arg(
        short = 'L',
        long = "level",
        help = "Logging level",
        value_name = "LEVEL",
        default_value_t = log::Level::Info,
    )]
    pub level: log::Level,
}

#[derive(Parser, Debug)]
pub enum Command {
    Run(RunArgs),
    Read(ReadArgs),
    /// Call transcript models from segmented reads by intron chain
    Chain(ChainArgs),
}

#[derive(Parser, Debug)]
pub struct ChainArgs {
    #[arg(
        long = "bed",
        required = true,
        value_name = "PATHS",
        value_delimiter = ',',
        num_args = 1..,
        help = "Segmented BED12 files of one [sample, chr] delimited by comma"
    )]
    pub bed: Vec<PathBuf>,

    #[arg(
        long = "ref",
        value_name = "PATH",
        help = "Reference transcripts (BED12)"
    )]
    pub reference: PathBuf,

    #[arg(long = "prefix", value_name = "PREFIX", help = "Output prefix")]
    pub prefix: String,

    #[arg(long = "preset", value_enum, default_value = "balanced")]
    pub preset: Preset,

    #[arg(
        long = "min-support-novel",
        value_name = "INT",
        help = "Molecules for a novel multi-exon chain [preset]"
    )]
    pub min_support_novel: Option<usize>,

    #[arg(
        long = "min-support-mono",
        value_name = "INT",
        help = "Molecules for a novel mono-exonic cluster [preset]"
    )]
    pub min_support_mono: Option<usize>,

    #[arg(
        long = "min-read-fraction",
        value_name = "FLOAT",
        help = "Share of locus reads kept novel models must cover [preset]"
    )]
    pub min_read_fraction: Option<f64>,

    #[arg(long = "junction-wobble", value_name = "INT", default_value_t = 5)]
    pub junction_wobble: u32,

    #[arg(
        long = "junction-support",
        help = "Drop a novel chain resting on one singleton molecule unless each of its novel introns is in another read; dropped reads go to <prefix>.rejected.bed",
        action = ArgAction::SetTrue
    )]
    pub junction_support: bool,

    // ponytail: accepted for the Nextflow contract and ignored; one [sample, chr] is
    // single-threaded work, rayon over (chrom, strand) buckets if it ever matters
    #[arg(long = "threads", value_name = "THREADS", default_value_t = 1)]
    pub threads: usize,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum Preset {
    Sensitive,
    Balanced,
    Strict,
}

#[derive(Parser, Debug)]
#[command(group(ArgGroup::new("indexing").required(false).multiple(true).args(&["index", "write"])))]
#[command(group(ArgGroup::new("mode").required(true).multiple(true).args(&["index", "extend"])))]
pub struct RunArgs {
    #[arg(
        short = 'b',
        long = "bed",
        required = true,
        value_name = "PATHS",
        value_delimiter = ',',
        num_args = 1..,
        help = "Paths to BED12 files delimited by comma"
    )]
    pub bed: Vec<PathBuf>,

    #[arg(
        short = 'I',
        long = "index",
        help = "Create an index for the given bed file(s)",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
    )]
    pub index: bool,

    #[arg(
        short = 'W',
        long = "write",
        help = "Writes the collapse counterpart of each bed file using the index (only if --index is set)",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
    )]
    pub write: bool,

    #[arg(
        short = 'E',
        long = "extend",
        help = "Add the index (queue of collapse reads) as an additional column to the bed file(s)",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
    )]
    pub extend: bool,

    #[arg(
        short = 'c',
        long = "collapse-mode",
        help = "Mode to collapse reads. Affects the way the binnary key is processed and thus any equivalence",
        value_name = "COLLAPSE_MODE",
        default_value = "transcript"
    )]
    pub collapse_mode: CollapseMode,

    #[arg(
        short = 'o',
        long = "outdir",
        help = "Path to output directory (/collapse)",
        value_name = "PATH",
        required = false,
        default_value = "."
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 'n',
        long = "name",
        help = "Name of the output file",
        value_name = "NAME",
        required = false,
        default_value = "collapsed.bed"
    )]
    pub name: String,

    #[arg(
        short = 'U',
        long = "max-five-utr-length",
        help = "Maximum 5'UTR distance between queues to allow collapsing reads in bp",
        value_name = "INT",
        default_value = "50"
    )]
    pub max_five_utr_length: u32,

    #[arg(
        short = 'u',
        long = "max-three-utr-length",
        help = "Maximum 3'UTR distance between queues to allow collapsing reads in bp",
        value_name = "INT",
        default_value = "500"
    )]
    pub max_three_utr_length: u32,
}

#[derive(Parser, Debug)]
#[command(group(ArgGroup::new("lookup").required(true).multiple(false).args(&["read", "write", "expand"])))]
pub struct ReadArgs {
    #[arg(
        short = 'i',
        long = "index",
        help = "Path to indext file produced by collapse",
        value_name = "PATH",
        required = true
    )]
    pub index: PathBuf,

    #[arg(
        short = 'r',
        long = "read",
        help = "Read name to be looked up in index",
        value_name = "NAME",
        required = false
    )]
    pub read: Option<String>,

    #[arg(
        short = 'W',
        long = "write",
        help = "Writes a verbose uncompressed index file",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
    )]
    pub write: bool,

    #[arg(
        short = 'E',
        long = "expand",
        help = "Writes the collapse counterpart of each bed file using the index (only if --index is set)",
        value_name = "FLAG",
        action = ArgAction::SetTrue,
    )]
    pub expand: bool,

    #[arg(
        short = 'o',
        long = "outdir",
        help = "Path to output directory (/collapse)",
        value_name = "PATH",
        required = false,
        default_value = ".",
        conflicts_with = "read"
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 'n',
        long = "name",
        help = "Name of the output file",
        value_name = "NAME",
        required = false,
        default_value = "expanded",
        conflicts_with = "read"
    )]
    pub name: String,
}

#[derive(Debug, Clone, Copy)]
pub enum RunMode {
    Index,
    Extend,
}

impl RunMode {
    pub fn from(extend: &bool) -> Self {
        if *extend {
            RunMode::Extend
        } else {
            RunMode::Index
        }
    }
}

#[derive(Parser, Debug)]
pub enum ReadMode {
    Lookup,
    Write,
    Expand,
}

impl ReadMode {
    pub fn from(args: &ReadArgs) -> Self {
        if args.write {
            ReadMode::Write
        } else if args.expand {
            ReadMode::Expand
        } else {
            ReadMode::Lookup
        }
    }
}

#[derive(Parser, Debug, Clone, ValueEnum, Copy)]
pub enum CollapseMode {
    #[clap(name = "exon")]
    Exon, // INFO: cds_start, cds_end are 0
    #[clap(name = "cds")]
    Coding, // INFO: tx_start, tx_end are 0
    #[clap(name = "transcript")]
    Transcript, // INFO: cannonical
    #[clap(name = "gapped-cds")]
    GappedCoding, // INFO: tx_start, tx_end are 0, all gaps preserved
}

impl std::str::FromStr for CollapseMode {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "exon" => Ok(CollapseMode::Exon),
            "cds" => Ok(CollapseMode::Coding),
            "transcript" => Ok(CollapseMode::Transcript),
            "gapped-cds" => Ok(CollapseMode::GappedCoding),
            _ => Err(format!("Invalid collapse mode: {s}")),
        }
    }
}
