#!/usr/bin/env python3
"""Run real Nextflow channel/process wiring with tiny stand-ins for sequencing tools.

Usage: python3 assets/ci/test_isoseq.py [path/to/nextflow]
"""

import json
import gzip
import os
import re
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
    # Real chunk names carry the --suffix (meta.id), so every chunk FASTA has a distinct staged name.
    Path('chunks').mkdir()
    shutil.copyfile(option('-f'), f"chunks/tmp_chunk_0_{option('--suffix')}.fasta.gz")
elif name == 'minimap2':
    Path(option('-o')).write_text(sam)
elif name == 'iso-cigar':
    stem = option('--bam').removesuffix('.bam')
    Path(stem + '.extended.bam').touch()
    Path(stem + '.extended.bam.bai').touch()
elif name == 'iso-align':
    assert all(Path(path).is_file() for path in [option('--bam')] + option('--reads').split(',')), args
    Path(option('--report')).touch()
    Path(option('--output')).touch()  # empty: the module deletes it, so no fragment re-alignment
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
    for name in ("pbindex", "lima", "isoseq", "samtools", "fxsplit", "minimap2", "iso-cigar", "iso-align",
                 "iso-segment", "iso-fusion"):
        executable = binary / name
        executable.write_text(TOOLS)
        executable.chmod(0o755)
    # Cigar extension stages genome and annotation side by side, so they need distinct names.
    for name in ("genome.fa", "genome.mmi", "annotation.bed"):
        (temporary / name).touch()

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
params.cigar = false
""")
    # The second pass runs only for the ark cases (--aligner ark); --cigar toggles cigar extension there.
    harness = temporary / "main.nf"
    harness.write_text("""
include { ISOSEQ } from 'REPO/src/subworkflows/isoseq/main.nf'
include { POOL_READS } from 'REPO/src/subworkflows/pool_reads/main.nf'
include { SPLIT_ALIGN_CLEAN_CHUNKS } from 'REPO/src/subworkflows/split_align/main.nf'
workflow {
    if (params.entrypoint == 'flnc') {
        // Same read channel as src/subworkflows/preprocessing/main.nf builds for flnc,
        // routed through POOL_READS like the real pipeline so modes are honored here too.
        reads = Channel.fromPath("${params.global_input_dir}/*.fast*", checkIfExists: true).map { fastx ->
            [[id: fastx.baseName, sample_id: fastx.name.replaceFirst(/(?:\\.(?:hq|singletons))?\\.fast[aq](?:\\.gz)?$/, ''),
              single_end: true, singleton: fastx.baseName.contains('singleton')], fastx]
        }
        POOL_READS(reads, params.isoseq_cluster2_mode, 'pooled')
        reads = POOL_READS.out.reads
        POOL_READS.out.reads.view { meta, f -> 'POOL\\t' + meta.id + '\\t' + meta.sample_id + '\\t' + meta.singleton + '\\t' + f.size() }
    } else {
        ISOSEQ(params.global_input_dir, params.global_primers, 1, params.isoseq_cluster2_mode,
            'pooled', params.entrypoint, params.entrypoint in ['refine', 'cluster'])
        reads = ISOSEQ.out.reads
    }
    SPLIT_ALIGN_CLEAN_CHUNKS(reads, Channel.value(file("${projectDir}/genome.fa")),
        Channel.value([[:], file("${projectDir}/genome.mmi")]), Channel.value([[:], file("${projectDir}/annotation.bed")]),
        Channel.value([[:], []]), params.aligner, false, false, false, params.cigar, params.aligner == 'ark',
        Channel.empty())
    SPLIT_ALIGN_CLEAN_CHUNKS.out.reads.view { meta, bed -> 'RESULT\\t' + meta.id }
}
""".replace("REPO", str(ROOT)))
    environment = dict(os.environ, PATH=f"{binary}:{os.environ['PATH']}", NXF_OFFLINE="true", NXF_ANSI_LOG="false")
    base = [NEXTFLOW, "-C", f"{ROOT}/src/nextflow.config,{config}", "run"]

    isoseqx = [f"IsoSeqX_bc{i:02}_5p--IsoSeqX_3p" for i in (1, 2)]
    neb = ["NEB_5p--NEB_Clontech_3p", "NEB_5p--primer_3p"]
    # cigar None: mm2, second pass off. True/False: ark, second pass on, cigar extension on/off.
    # flnc honors the cluster mode through POOL_READS; "per_sample" only makes `expected` below the per-sample ids.
    for entrypoint, mode, pairs, cigar in (
            ("refine", "per_sample", isoseqx, None), ("refine", "multi_sample", isoseqx, None),
            ("refine", "both", isoseqx, None), ("refine", "per_sample", isoseqx[:1], None),
            ("ccs", "per_sample", isoseqx[:1], None), ("ccs", "both", isoseqx, None),
            ("refine", "both", neb, None), ("ccs", "both", neb, None), ("cluster", "both", isoseqx, None),
            ("refine", "multi_sample", isoseqx, True), ("refine", "per_sample", isoseqx, True),
            ("flnc", "per_sample", isoseqx, True), ("flnc", "per_sample", isoseqx, False),
            ("flnc", "multi_sample", isoseqx, True), ("flnc", "both", isoseqx, True)):
        case = temporary / f"{entrypoint}-{mode}-{pairs[0]}-{len(pairs)}-{cigar}"
        inputs = case / "02_LIMA"
        inputs.mkdir(parents=True)
        samples = {f"movie.part.hifi.{pair}" for pair in pairs}
        files = {"refine": [f"movie.part.hifi_fl.{pair}.bam" for pair in pairs], "ccs": ["movie.part.hifi.bam"],
                 "cluster": [f"movie.part.hifi.{pair}_flnc.bam" for pair in pairs],
                 "flnc": [f"{sample}.{kind}.fasta.gz" for sample in samples for kind in ("hq", "singletons")]}
        file_bytes = {}
        for i, name in enumerate(files[entrypoint]):
            if entrypoint == "flnc":
                # Real gzipped records so pooled concatenation is byte-assertable below.
                file_bytes[name] = gzip.compress(f">{name}\nACGT\n".encode())
                (inputs / name).write_bytes(file_bytes[name])
            else:
                (inputs / name).touch()
                if i == 0:
                    (inputs / f"{name}.pbi").touch()
        (inputs / "ignored.consensusreadset.xml").touch()
        (inputs / "ignored.lima.report").touch()
        primers = case / "primers.fasta"
        names = dict.fromkeys(primer for pair in pairs for primer in pair.split("--"))
        primers.write_text("".join(f">{name}\nACGT\n" for name in names))
        calls = case / "calls.jsonl"
        environment.update(TEST_PRIMER_PAIRS=json.dumps(pairs), TEST_CALLS=str(calls))
        command = base + [str(harness), "--global_input_dir", str(inputs),
                          "--global_output_dir", str(case / "results"), "--isoseq_cluster2_mode", mode,
                          "--entrypoint", entrypoint]
        command += ["--global_primers", str(primers)] if entrypoint != "cluster" else []  # cluster needs none
        command += ["--aligner", "ark", "--cigar", str(cigar).lower()] if cigar is not None else []
        result = subprocess.run(command, cwd=case, env=environment, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=120)
        assert result.returncode == 0, result.stdout
        actual = {line.split("\t")[1] for line in result.stdout.splitlines() if line.startswith("RESULT\t")}
        expected = (samples if mode != "multi_sample" else set()) | ({"pooled"} if mode != "per_sample" else set())
        assert actual == expected, (actual, expected, result.stdout)
        if entrypoint == "flnc":
            # Pooled items carry sample_id 'pooled' with ids distinct from per-sample ones;
            # pooled bytes are exactly the member files concatenated in sorted-name order.
            pools = [line.split("\t")[1:] for line in result.stdout.splitlines() if line.startswith("POOL\t")]
            by_id = {pid: (sample, singleton == "true", int(size)) for pid, sample, singleton, size in pools}
            members = {}
            for name in files["flnc"]:
                kind = "singletons" if ".singletons." in name else "hq"
                members.setdefault(kind, []).append(name)
                if mode != "multi_sample":
                    per_sample_id = name.removesuffix(".gz")
                    sample = re.sub(r"(?:\.(?:hq|singletons))?\.fast[aq](?:\.gz)?$", "", name)
                    assert by_id.get(per_sample_id) == (sample, kind == "singletons", len(file_bytes[name])), (by_id, name)
            pooled_ids = {"pooled.pooled", "pooled.pooled.singletons"}
            if mode == "per_sample":
                assert not (set(by_id) & pooled_ids), by_id
            else:
                for kind, cls in (("hq", ""), ("singletons", ".singletons")):
                    pid = f"pooled.pooled{cls}"
                    want = b"".join(file_bytes[n] for n in sorted(members[kind]))
                    assert by_id.get(pid) == ("pooled", kind == "singletons", len(want)), (by_id, pid)
                    hits = [p for p in case.rglob(f"{pid}.fasta.gz")]
                    assert hits, (pid, mode)
                    for path in hits:
                        assert gzip.decompress(path.read_bytes()) == gzip.decompress(want), (path, mode)
        records = [json.loads(line) for line in calls.read_text().splitlines()]
        refine = [args for name, args in records if name == "isoseq" and args[0] == "refine"]
        cluster = [args for name, args in records if name == "isoseq" and args[0] == "cluster2"]
        assert len(refine) == (len(pairs) if entrypoint in ("refine", "ccs") else 0), refine
        assert len(cluster) == (len(expected) if entrypoint != "flnc" else 0), cluster
        assert sum(name == "lima" for name, _ in records) == (entrypoint == "ccs")
        assert sum(name == "pbindex" for name, _ in records) == (entrypoint in ("refine", "cluster") and len(pairs) > 1)
        assert not any(name == "samtools" and args[0] == "merge" for name, args in records), records
        # Every result id has one hq and one singleton chunk. With the second pass on, iso-align runs once per
        # chunk BAM with exactly the chunk FASTA minimap2 aligned it from (reads name -> SAM name).
        chunks = {args[-3]: args[-1].removesuffix(".sam") for name, args in records if name == "minimap2"}
        found = sorted((args[args.index("--bam") + 1], args[args.index("--reads") + 1])
                       for name, args in records if name == "iso-align")
        suffix = ".extended.bam" if cigar else ".bam"
        assert len(chunks) == 2 * len(expected), chunks
        assert found == (sorted((sam + suffix, reads) for reads, sam in chunks.items()) if cigar is not None else []), (found, chunks)
        second_pass = "" if cigar is None else f", ark second pass, cigar extension {'on' if cigar else 'off'}"
        print(f"PASS {entrypoint}: {mode}, {', '.join(pairs)}{second_pass}", flush=True)

    # A checkpoint without a primer-pair (refine) or _flnc (cluster) suffix must fail before any tool runs.
    for start, bam, diagnostic in (("refine", "movie.part.hifi_fl.NEB_5p.bam", "Unexpected LIMA BAM name"),
                                   ("cluster", "movie.part.hifi.NEB_5p--NEB_Clontech_3p.bam",
                                    "Unexpected refined BAM name")):
        case = temporary / f"invalid-{start}-name"
        case.mkdir()
        (case / bam).touch()
        (case / f"{bam}.pbi").touch()
        calls = case / "calls.jsonl"
        environment.update(TEST_CALLS=str(calls))
        result = subprocess.run(base + [str(harness), "--entrypoint", start, "--global_input_dir", str(case),
                                       "--global_primers", str(primers), "--global_output_dir", str(case / "results")],
                                cwd=case, env=environment, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=120)
        assert result.returncode != 0 and f"{diagnostic}: {bam}" in result.stdout, result.stdout
        assert not calls.exists(), calls.read_text()
    print("PASS malformed LIMA and refined filenames rejected before any tool runs", flush=True)

    # Exercise the actual CLI validator: refine needs primers, cluster does not, unknown entrypoints fail.
    for start, diagnostic in (("refine", "missing required --global_primers"),
                              ("cluster", "Parameter validation failed"),
                              ("unknown", "Unknown entrypoint option")):
        result = subprocess.run(base + [str(ROOT / "src/main.nf"), "--entrypoint", start],
                                cwd=temporary, env=environment, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=120)
        assert result.returncode != 0 and diagnostic in result.stdout, result.stdout
        assert ("missing required --global_primers" in result.stdout) == (start == "refine"), result.stdout
    print("PASS --entrypoint validation and missing primers", flush=True)
