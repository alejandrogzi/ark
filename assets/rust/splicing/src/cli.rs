// Copyright (c) 2026 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use clap::{Parser, ValueEnum};
use log::Level;
use std::path::PathBuf;

/// Minimum SpliceAI score to recover from BigWig files.
pub const SPLICE_AI_SCORE_RECOVERY_THRESHOLD: f32 = 0.001;
/// Default minimum derived score value.
pub const DEFAULT_SCORE_FLOOR: i32 = -4;
/// Default maximum derived score value.
pub const DEFAULT_SCORE_CEILING: i32 = 13;
/// Default synthetic SpliceAI score for splice sites included from regions.
pub const DEFAULT_SS_REGION_SCORE: f32 = 0.5;
/// Default bonus for annotation splice sites already present in BigWig files.
pub const DEFAULT_BONUS_FOR_SS_REGIONS: f32 = 0.25;

/// Region feature class used to select splice sites for synthetic inclusion.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum SsRegionPosition {
    /// Use splice sites from all exon junctions.
    Exon,
    /// Use splice sites from coding exon junctions only.
    Cds,
}

impl std::fmt::Display for SsRegionPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SsRegionPosition::Exon => write!(f, "exon"),
            SsRegionPosition::Cds => write!(f, "cds"),
        }
    }
}

/// Parses and validates a SpliceAI score-like value.
pub fn parse_score_for_ss_regions(value: &str) -> Result<f32, String> {
    let score = value
        .parse::<f32>()
        .map_err(|e| format!("invalid floating-point score: {e}"))?;

    if !score.is_finite() {
        return Err("score must be finite".to_string());
    }

    if !(0.0..=1.0).contains(&score) {
        return Err("score must be between 0.0 and 1.0".to_string());
    }

    Ok(score)
}

#[derive(Debug, Parser)]
#[command(name = "splicing", about = "derive spliceAi scores", version = env!("CARGO_PKG_VERSION"))]
pub struct Args {
    #[arg(
        short = 'b',
        long = "bigwig-dir",
        required = true,
        value_name = "PATH",
        help = "Path to directory containing SpliceAI BigWig files; filenames must include donor/acceptor and plus/minus tokens"
    )]
    pub bw_dir: PathBuf,

    #[clap(
        short = 's',
        long = "sequence",
        help = "Path to sequence file (FASTA/2bit, use '-' or omit to read stdin)",
        value_name = "SEQUENCE",
        default_value = "-"
    )]
    pub sequence: PathBuf,

    #[clap(
        short = 'r',
        long = "regions",
        help = "Path to regions file (BED/GTF/GFF/GZ)",
        value_name = "REGIONS",
        required = true
    )]
    pub regions: PathBuf,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[arg(
        short = 'l',
        long = "level",
        help = "Log level",
        value_name = "LEVEL",
        default_value_t = log::Level::Info,
    )]
    pub level: Level,

    #[arg(
        short = 'p',
        long = "prefix",
        required = false,
        value_name = "PATH",
        help = "Prefix for output files",
        default_value_t = String::from("spliceai")
    )]
    pub prefix: String,

    #[arg(
        short = 'o',
        long = "outdir",
        required = false,
        value_name = "PATH",
        help = "Output directory",
        default_value = "."
    )]
    pub output_dir: PathBuf,

    #[arg(
        long = "floor",
        required = false,
        value_name = "INT",
        help = "Minimum final rounded derived score",
        default_value_t = DEFAULT_SCORE_FLOOR
    )]
    pub floor: i32,

    #[arg(
        long = "ceiling",
        required = false,
        value_name = "INT",
        help = "Maximum final rounded derived score",
        default_value_t = DEFAULT_SCORE_CEILING
    )]
    pub ceiling: i32,

    #[arg(
        short = 'S',
        long = "include-ss-from-regions",
        required = false,
        help = "Include missing splice sites from --regions in the scored splice-site maps"
    )]
    pub include_ss_from_regions: bool,

    #[arg(
        short = 'K',
        long = "score-for-ss-regions",
        required = false,
        value_name = "FLOAT",
        help = "Synthetic SpliceAI score for missing splice sites from --regions when --include-ss-from-regions is set",
        default_value_t = DEFAULT_SS_REGION_SCORE,
        value_parser = parse_score_for_ss_regions
    )]
    pub score_for_ss_regions: f32,

    #[arg(
        short = 'B',
        long = "bonus",
        required = false,
        help = "Add --bonus-score to annotated splice sites already present in BigWigs when --include-ss-from-regions is set"
    )]
    pub bonus: bool,

    #[arg(
        short = 'J',
        long = "bonus-score",
        required = false,
        value_name = "FLOAT",
        help = "Bonus added to annotated splice sites already present in BigWigs when --bonus is set",
        default_value_t = DEFAULT_BONUS_FOR_SS_REGIONS,
        value_parser = parse_score_for_ss_regions
    )]
    pub bonus_score: f32,

    #[arg(
        short = 'R',
        long = "region-floor",
        required = false,
        value_name = "INT",
        requires = "include_ss_from_regions",
        conflicts_with = "bonus",
        help = "With --include-ss-from-regions: after calibration, give every splice site from --regions a derived score of at least INT (sites missing from the BigWigs are added at INT). Unlike --bonus and --score-for-ss-regions, which change raw scores before calibration, this leaves every other site's calibration untouched"
    )]
    pub region_floor: Option<i32>,

    #[arg(
        short = 'P',
        long = "position-for-ss-regions",
        required = false,
        value_name = "STRING",
        help = "Use splice sites from all exon junctions or coding exon junctions only when --include-ss-from-regions is set",
        value_enum,
        default_value_t = SsRegionPosition::Exon
    )]
    pub position_for_ss_regions: SsRegionPosition,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_region_splice_site_options() {
        let args = Args::try_parse_from([
            "splicing",
            "-b",
            "bw",
            "-r",
            "regions.gtf",
            "-S",
            "-K",
            "0.7",
            "-B",
            "-J",
            "0.3",
            "-P",
            "cds",
        ])
        .unwrap();

        assert!(args.include_ss_from_regions);
        assert_eq!(args.score_for_ss_regions, 0.7);
        assert!(args.bonus);
        assert_eq!(args.bonus_score, 0.3);
        assert_eq!(args.position_for_ss_regions, SsRegionPosition::Cds);
    }

    #[test]
    fn defaults_region_splice_site_bonus_options() {
        let args = Args::try_parse_from(["splicing", "-b", "bw", "-r", "regions.gtf"]).unwrap();

        assert!(!args.bonus);
        assert_eq!(args.bonus_score, DEFAULT_BONUS_FOR_SS_REGIONS);
    }

    #[test]
    fn rejects_invalid_region_splice_site_score() {
        for score in ["-0.1", "1.1", "NaN"] {
            let result =
                Args::try_parse_from(["splicing", "-b", "bw", "-r", "regions.gtf", "-K", score]);
            assert!(result.is_err(), "accepted invalid score {score}");
        }
    }

    #[test]
    fn rejects_invalid_region_splice_site_bonus_score() {
        for score in ["-0.1", "1.1", "NaN"] {
            let result =
                Args::try_parse_from(["splicing", "-b", "bw", "-r", "regions.gtf", "-J", score]);
            assert!(result.is_err(), "accepted invalid bonus score {score}");
        }
    }

    #[test]
    fn rejects_invalid_region_splice_site_position() {
        let result =
            Args::try_parse_from(["splicing", "-b", "bw", "-r", "regions.gtf", "-P", "intron"]);

        assert!(result.is_err());
    }
}
