#!/usr/bin/env python3
"""Run real Nextflow channel/process wiring with tiny stand-ins for sequencing tools.

Usage: python3 assets/ci/test_isoseq.py [path/to/nextflow]
"""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[2]
NEXTFLOW = sys.argv[1] if len(sys.argv) > 1 else "nextflow"
TOOLS = r'''#!/usr/bin/env python3
import json, os, shutil, sys
from pathlib import Path
name, args = Path(sys.argv[0]).name, sys.argv[1:]
if '--version' in args:
    print(name + ' ' + (' '.join(args[:-1]) + ' ' if name == 'isoseq' else '') + 'test')
    sys.exit(0)
pairs = json.loads(os.environ['TEST_PRIMER_PAIRS'])
sam = '@HD\tVN:1.6\n' + ''.join(f'r{i}\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\tis:i:{i}\n' for i in (1, 2))
def outputs(bam):
    Path(bam).write_text(sam)
    for suffix in ('.bam.pbi', '.consensusreadset.xml', '.filter_summary.report.json', '.report.csv'):
        Path(bam.removesuffix('.bam') + suffix).touch()
def option(flag):
    return args[args.index(flag) + 1]
if name == 'pbindex':
    Path(args[-1] + '.pbi').touch()
elif name == 'lima':
    assert Path(args[0] + '.pbi').is_file()
    stem = args[2].removesuffix('.bam')
    for pair in pairs:
        outputs(f'{stem}.{pair}.bam')
    for suffix in ('counts', 'report', 'summary'):
        Path(stem + '.lima.' + suffix).touch()
elif name == 'isoseq' and args[0] == 'refine':
    bam, primers, output = args[-3:]
    assert Path(bam).is_file() and Path(bam + '.pbi').is_file() and Path(primers).is_file()
    outputs(output)
elif name == 'isoseq' and args[0] == 'cluster2':
    bams = Path(args[1]).read_text().splitlines()
    assert bams and all(Path(bam).is_file() and bam.endswith('_flnc.bam') for bam in bams), bams
    assert len(bams) == (len(pairs) if args[2].startswith('pooled') else 1), bams
    outputs(args[2])
    Path(args[2].removesuffix('.bam') + '.cluster_report.csv').touch()
elif name == 'samtools':
    if args[0] == 'view' and '-bo' in args:
        Path(option('-bo')).write_text(sys.stdin.read())
    elif args[0] == 'view':
        print(sam, end='')
    elif args[0] == 'sort':
        Path(option('-o')).write_text(sys.stdin.read())
    elif args[0] == 'index':
        Path(args[-1] + '.bai').touch()
    elif args[0] == 'fasta':
        for line in Path(args[-1]).read_text().splitlines():
            if not line.startswith('@'):
                print('>' + line.split('\t')[0] + '\nACGT')
elif name == 'fxsplit':
    Path('chunks').mkdir()
    shutil.copyfile(option('-f'), 'chunks/0.fasta.gz')
elif name == 'minimap2':
    Path(option('-o')).write_text(sam)
elif name == 'iso-segment':
    Path('chr1@' + option('--prefix') + '.hq.bed').write_text('chr1\t0\t4\tr1\n')
elif name == 'iso-fusion':
    beds = option('--query').split(',')
    assert len(beds) == 2 and all(Path(bed).is_file() for bed in beds), beds
    directory = Path(option('--prefix'))
    directory.mkdir()
    (directory / 'fusions.free.bed').write_text('chr1\t0\t4\tr1\n')
else:
    raise AssertionError((name, args))
with open(os.environ['TEST_CALLS'], 'a') as log:
    log.write(json.dumps([name, args]) + '\n')
'''


with tempfile.TemporaryDirectory(prefix="ark-isoseq-") as temporary:
    temporary = Path(temporary)
    binary = temporary / "bin"
    binary.mkdir()
    for name in ("pbindex", "lima", "isoseq", "samtools", "fxsplit", "minimap2", "iso-segment", "iso-fusion"):
        executable = binary / name
        executable.write_text(TOOLS)
        executable.chmod(0o755)

    config = temporary / "test.config"
    config.write_text("""
process {
    withName: '.*' {
        cpus = 1
        memory = '256 MB'
        time = '1 min'
        errorStrategy = 'terminate'
    }
}
docker.enabled = false
apptainer.enabled = false
singularity.enabled = false
conda.enabled = false
params.aligner = 'mm2'
""")
    harness = temporary / "main.nf"
    harness.write_text("""
include { ISOSEQ } from 'REPO/src/subworkflows/isoseq/main.nf'
include { SPLIT_ALIGN_CLEAN_CHUNKS } from 'REPO/src/subworkflows/split_align/main.nf'
workflow {
    ISOSEQ(params.global_input_dir, params.global_primers, 1, params.isoseq_cluster2_mode,
        'pooled', params.entrypoint, params.entrypoint == 'refine')
    SPLIT_ALIGN_CLEAN_CHUNKS(ISOSEQ.out.reads, Channel.value(file(params.global_primers)),
        Channel.value([[:], file(params.global_primers)]),
        Channel.value([[:], file(params.global_primers)]), 'pooled', Channel.value([[:], []]),
        params.isoseq_cluster2_mode, params.entrypoint, 'mm2', false,
        false, false, false, false, Channel.empty())
    SPLIT_ALIGN_CLEAN_CHUNKS.out.reads.view { meta, bed -> 'RESULT\\t' + meta.id }
}
""".replace("REPO", str(ROOT)))
    environment = dict(os.environ, PATH=f"{binary}:{os.environ['PATH']}", NXF_OFFLINE="true", NXF_ANSI_LOG="false")
    base = [NEXTFLOW, "-C", f"{ROOT}/src/nextflow.config,{config}", "run"]

    isoseqx = [f"IsoSeqX_bc{i:02}_5p--IsoSeqX_3p" for i in (1, 2)]
    neb = ["NEB_5p--NEB_Clontech_3p", "NEB_5p--primer_3p"]
    for entrypoint, mode, pairs in (("refine", "per_sample", isoseqx), ("refine", "multi_sample", isoseqx),
                                    ("refine", "both", isoseqx), ("refine", "per_sample", isoseqx[:1]),
                                    ("ccs", "per_sample", isoseqx[:1]), ("ccs", "both", isoseqx),
                                    ("refine", "both", neb), ("ccs", "both", neb)):
        case = temporary / f"{entrypoint}-{mode}-{pairs[0]}-{len(pairs)}"
        inputs = case / "02_LIMA"
        inputs.mkdir(parents=True)
        stems = [f"movie.part.hifi_fl.{pair}" for pair in pairs]
        for i, stem in enumerate(stems if entrypoint == "refine" else ["movie.part.hifi"]):
            (inputs / f"{stem}.bam").touch()
            if i == 0:
                (inputs / f"{stem}.bam.pbi").touch()
        (inputs / "ignored.consensusreadset.xml").touch()
        (inputs / "ignored.lima.report").touch()
        primers = case / "primers.fasta"
        names = dict.fromkeys(primer for pair in pairs for primer in pair.split("--"))
        primers.write_text("".join(f">{name}\nACGT\n" for name in names))
        calls = case / "calls.jsonl"
        environment.update(TEST_PRIMER_PAIRS=json.dumps(pairs), TEST_CALLS=str(calls))
        command = base + [str(harness), "--global_input_dir", str(inputs), "--global_primers", str(primers),
                          "--global_output_dir", str(case / "results"), "--isoseq_cluster2_mode", mode,
                          "--entrypoint", entrypoint]
        result = subprocess.run(command, cwd=case, env=environment, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=120)
        assert result.returncode == 0, result.stdout
        actual = {line.split("\t")[1] for line in result.stdout.splitlines() if line.startswith("RESULT\t")}
        samples = {f"movie.part.hifi.{pair}" for pair in pairs}
        expected = (samples if mode != "multi_sample" else set()) | ({"pooled"} if mode != "per_sample" else set())
        assert actual == expected, (actual, expected, result.stdout)
        records = [json.loads(line) for line in calls.read_text().splitlines()]
        refine = [args for name, args in records if name == "isoseq" and args[0] == "refine"]
        cluster = [args for name, args in records if name == "isoseq" and args[0] == "cluster2"]
        assert len(refine) == len(pairs), refine
        assert len(cluster) == len(expected), cluster
        assert sum(name == "lima" for name, _ in records) == (entrypoint == "ccs")
        assert sum(name == "pbindex" for name, _ in records) == (entrypoint == "refine" and len(pairs) > 1)
        print(f"PASS {entrypoint}: {mode}, {', '.join(pairs)}", flush=True)

    # A checkpoint without a primer-pair suffix must fail before refinement.
    case = temporary / "invalid-lima-name"
    case.mkdir()
    (case / "movie.part.hifi_fl.NEB_5p.bam").touch()
    (case / "movie.part.hifi_fl.NEB_5p.bam.pbi").touch()
    calls = case / "calls.jsonl"
    environment.update(TEST_CALLS=str(calls))
    result = subprocess.run(base + [str(harness), "--entrypoint", "refine", "--global_input_dir", str(case),
                                   "--global_primers", str(primers)],
                            cwd=case, env=environment, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=120)
    assert result.returncode != 0 and "Unexpected LIMA BAM name: movie.part.hifi_fl.NEB_5p.bam" in result.stdout, result.stdout
    assert not calls.exists(), calls.read_text()
    print("PASS malformed LIMA filename rejected before refinement", flush=True)

    # Exercise the actual CLI validator for the refine checkpoint and invalid entrypoints.
    for start, diagnostic in (("refine", "missing required --global_primers"),
                              ("unknown", "Unknown entrypoint option")):
        result = subprocess.run(base + [str(ROOT / "src/main.nf"), "--entrypoint", start],
                                cwd=temporary, env=environment, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=120)
        assert result.returncode != 0 and diagnostic in result.stdout, result.stdout
    print("PASS --entrypoint validation and missing primers", flush=True)
