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
//!
//! # Benchmark (`chain`, 0.1.1)
//!
//! LRGASP WTC11 rep1 SIRV-Set 4, de novo (empty annotation). Truth is the 60 multi-exon SIRVs. A model is a
//! hit when its intron chain matches a truth chain exactly; precision is over multi-exon models. Every tool
//! gets the same 20,465 reads that reach ARK's `CHAIN_COLLAPSE`, as the reads BED for collapse and as
//! the matching subset of ARK's minimap2 BAM for the others. Coordinates and chains are identical.
//!
//! - collapse: `chain --preset P [--junction-support] --ref empty.bed`;
//! - IsoQuant 4.0.0: `--data_type pacbio_ccs`, each `--model_construction_strategy`;
//! - isocall 1.3.0: `prep-isoforms` (empty GTF), `profile`, `merge`, `call` with its presets and one filter changed per config;
//! - cluster2 (isoseq 4.0.0, `--singletons`) on an unaligned BAM rebuilt from the same reads: one record per HQ
//!   cluster plus `#SG` singletons, which then go through collapse.
//!
//! | Caller | Multi-exon models | Recall | Precision | F1 |
//! |---|---|---|---|---|
//! | none (every read) | 11,728 | 0.983 | 0.935 (per read) | — |
//! | collapse sensitive | 190 | 0.983 | 0.311 | 0.472 |
//! | collapse sensitive + junction-support | 109 | 0.983 | 0.541 | 0.698 |
//! | collapse balanced (± junction-support) | 73 | 0.967 | 0.795 | 0.872 |
//! | collapse strict (± junction-support) | 54 | 0.850 | 0.944 | 0.895 |
//! | IsoQuant default / default_pacbio | 46 | 0.700 | 0.913 | 0.792 |
//! | IsoQuant sensitive_pacbio | 47 | 0.717 | 0.915 | 0.804 |
//! | IsoQuant fl_pacbio | 43 | 0.683 | 0.953 | 0.796 |
//! | IsoQuant reliable | 36 | 0.567 | 0.944 | 0.708 |
//! | isocall default | 40 | 0.633 | 0.950 | 0.760 |
//! | isocall min reads 2 / no internal-priming filter | 45 | 0.700 | 0.933 | 0.800 |
//! | isocall no relative-abundance filter | 44 | 0.650 | 0.886 | 0.750 |
//! | isocall yolo | 61 | 0.717 | 0.705 | 0.711 |
//! | isocall min reads 1 | 72 | 0.750 | 0.625 | 0.682 |
//! | isocall yolo + min reads 1 | 259 | 0.867 | 0.201 | 0.326 |
//! | cluster2 output | 1,562 | 0.983 | 0.782 | 0.871 |
//! | cluster2 → collapse sensitive | 159 | 0.967 | 0.365 | 0.530 |
//! | cluster2 → collapse sensitive + junction-support | 87 | 0.967 | 0.667 | 0.789 |
//! | cluster2 → collapse balanced (± junction-support) | 69 | 0.967 | 0.841 | 0.899 |
//!
//! `--junction-support` on giraffe muscle: two ARK runs from the same 642,806 subreads with cluster2, identical
//! except for the option. Models go from 5,250 to 3,578 and POLISH passes from 3,192 to 2,309. All 1,628
//! known models are kept, as are 829/831 annotated-intron and 373/375 independently supported extras. 911/924
//! flagged and 122/129 non-canonical extras are dropped, along with 624/642 unconfirmed canonical extras.

pub mod chain;
pub mod cli;
pub mod record;
pub mod utils;

pub mod read;
