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

use collapse::{
    cli::{Args, Command, ReadMode, RunMode},
    read::{__expand, __lookup},
    utils::unpack,
    utils::{__write_collapsed, __write_index},
};

use clap::Parser;
use simple_logger::init_with_level;

fn main() {
    let start = std::time::Instant::now();

    let args = Args::parse();
    init_with_level(args.level).unwrap_or_else(|e| log::error!("ERROR: Logger init error -> {e}"));

    match args.command {
        Command::Run(args) => {
            let mode = RunMode::from(&args.extend);

            let target = args.outdir.join("collapsed");
            std::fs::create_dir_all(&target).unwrap_or_else(|e| {
                panic!("ERROR: Could not create directory {:?} -> {e}", &target)
            });

            let tracks = unpack(
                args.bed.clone(),
                args.collapse_mode,
                args.max_five_utr_length,
                args.max_three_utr_length,
            )
            .unwrap_or_else(|e| {
                panic!(
                    "ERROR: There was an error while unpacking {:?} -> {e}",
                    &args.bed
                )
            });

            match mode {
                RunMode::Extend => {
                    let output = target.join(args.name);
                    log::info!("INFO: Writing to collapsed file to {:?}", &output);

                    // INFO: extend writes collapsed bed file [merged], --write will write index
                    if args.write {
                        let index = target.join("index");
                        log::info!("INFO: Writing index to {:?}", &index);

                        __write_index(&tracks, &index);
                    }

                    __write_collapsed(tracks, mode, &output); // WARN: extra column appended!
                }
                RunMode::Index => {
                    // INFO: index write by default, --write will write collapsed file [merged]
                    let index = target.join("index");
                    log::info!("INFO: Writing index to {:?}", &index);
                    __write_index(&tracks, &index);

                    if args.write {
                        let output = target.join(args.name);
                        log::info!("INFO: Writing to collapsed file to {:?}", &output);

                        __write_collapsed(tracks, mode, output); // WARN: no extra column appended!
                    }
                }
            }
        }
        Command::Read(args) => {
            let mode = ReadMode::from(&args);

            match mode {
                ReadMode::Lookup => {
                    __lookup(&args.index, args.read);
                }
                ReadMode::Write => todo!(),
                ReadMode::Expand => {
                    let output = args
                        .outdir
                        .join("expanded")
                        .join(args.name)
                        .with_extension("bed");
                    log::info!("INFO: Writing to expanded file to {:?}", &output);

                    __expand(args.index, output);
                }
            };
        }
        Command::Chain(args) => {
            collapse::chain::run(args).unwrap_or_else(|e| panic!("ERROR: {e}"));
        }
    }

    let elapsed = start.elapsed();
    log::info!("INFO: Elapsed time: {elapsed:?}");
}
