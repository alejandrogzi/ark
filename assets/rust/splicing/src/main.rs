// Copyright (c) 2026 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use clap::Parser;
use log::error;
use simple_logger::init_with_level;

use splicing::{
    annotation::get_ss_from_annotation,
    calibration::{build_null_profile, derive_spliceai_score, make_bins, RegionFloor},
    cli::Args,
    genome::get_dinucleotide_count,
    spliceai::{get_splice_scores, include_region_splice_sites},
};

fn main() {
    let args = Args::parse();
    init_with_level(args.level).unwrap_or_else(|_| panic!("ERROR: Cannot initialize logger!"));

    if args.floor > args.ceiling {
        error!(
            "ERROR: floor ({}) cannot be greater than ceiling ({})",
            args.floor, args.ceiling
        );
        std::process::exit(1);
    }
    if let Some(region_floor) = args.region_floor {
        if region_floor < args.floor || region_floor > args.ceiling {
            error!(
                "ERROR: region floor ({}) must lie within [floor, ceiling] = [{}, {}]",
                region_floor, args.floor, args.ceiling
            );
            std::process::exit(1);
        }
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap_or_else(|e| {
            error!("ERROR: failed to initialize thread pool: {}", e);
            std::process::exit(1);
        });

    let (genome, dinucleotide_count) = get_dinucleotide_count(&args.sequence).unwrap_or_else(|e| {
        error!("ERROR: failed to get dinucleotide counts: {}", e);
        std::process::exit(1);
    });
    let chroms = genome.keys().cloned().collect::<Vec<Vec<u8>>>();

    let region_position = args
        .include_ss_from_regions
        .then_some(args.position_for_ss_regions);
    let (donors, acceptors, ss_coord_donors, ss_coord_acceptors, region_donors, region_acceptors) =
        get_ss_from_annotation(&genome, args.regions, region_position);

    let donor_profile = build_null_profile("donor", &donors, &dinucleotide_count);
    let acceptor_profile = build_null_profile("acceptor", &acceptors, &dinucleotide_count);

    let (spliceai_donor_scores, spliceai_acceptor_scores) =
        get_splice_scores(args.bw_dir, chroms, genome).unwrap_or_else(|e| {
            error!("ERROR: failed to load SpliceAI BigWigs: {}", e);
            std::process::exit(1);
        });

    // INFO: with --region-floor the region sites are applied after calibration (in derive_spliceai_score);
    // INFO: otherwise they enter the raw maps here, before calibration
    if args.include_ss_from_regions && args.region_floor.is_none() {
        include_region_splice_sites(
            &spliceai_donor_scores,
            &spliceai_acceptor_scores,
            &region_donors,
            &region_acceptors,
            args.score_for_ss_regions,
            args.bonus,
            args.bonus_score,
        );
    }

    let donor_bins = make_bins(
        &spliceai_donor_scores,
        &ss_coord_donors,
        &donors,
        &dinucleotide_count,
    );
    let acceptor_bins = make_bins(
        &spliceai_acceptor_scores,
        &ss_coord_acceptors,
        &acceptors,
        &dinucleotide_count,
    );

    derive_spliceai_score(
        &spliceai_donor_scores,
        &spliceai_acceptor_scores,
        &donor_bins,
        &acceptor_bins,
        &donor_profile,
        &acceptor_profile,
        args.floor,
        args.ceiling,
        args.region_floor.map(|score| RegionFloor {
            donors: &region_donors,
            acceptors: &region_acceptors,
            score,
        }),
        &args.prefix,
        &args.output_dir,
    )
    .unwrap_or_else(|e| {
        error!("ERROR: failed to derive splice scores: {}", e);
        std::process::exit(1);
    });
}
