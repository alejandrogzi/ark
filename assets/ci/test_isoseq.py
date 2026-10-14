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
import gzip, json, os, shutil, sys
from pathlib import Path
name, args = Path(sys.argv[0]).name, sys.argv[1:]
if '--version' in args:
    print(name + ' ' + (' '.join(args[:-1]) + ' ' if name == 'isoseq' else '') + 'test')
    sys.exit(0)
if name == 'iso-align' and args == ['--help']:
    print('--bam --realigned-bam --rules')
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
    fastx = args[0].endswith('.fastq.gz')  # LIMA_FASTX: FASTQ in, one FASTQ per primer pair out
    assert fastx or Path(args[0] + '.pbi').is_file()
    # lima 26.2.1 hangs on non-PacBio read names, so FASTX_RENAME must have given CCS names
    assert not fastx or gzip.open(args[0], 'rt').readline().split()[0].endswith('/ccs'), args[0]
    stem = args[2].removesuffix('.fastq.gz' if fastx else '.bam')
    for pair in pairs:
        if fastx:
            shutil.copyfile(args[0], f'{stem}.{pair}.fastq.gz')
        else:
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
    if args[0] == 'view' and '-N' in args:
        # SAMTOOLS_SELECT_READS: keep (-o) or drop (-U) the named reads; the stub copies the input
        assert Path(option('-N')).is_file() and Path(args[-1]).is_file(), args
        shutil.copyfile(args[-1], option('-U') if '-U' in args else option('-o'))
    elif args[0] == 'view' and '-c' in args:
        print(1)
    elif args[0] == 'view' and '-bo' in args:
        Path(option('-bo')).write_text(sys.stdin.read())
    elif args[0] == 'view':
        print(sam, end='')
    elif args[0] == 'sort':
        # minimap2 pipes into sort (stdin); SAMTOOLS_SELECT_READS sorts a file argument
        source = args[-1] if args[-1].endswith('.bam') and args[-1] != option('-o') else None
        Path(option('-o')).write_text(Path(source).read_text() if source else sys.stdin.read())
    elif args[0] == 'index':
        Path(args[-1] + '.bai').touch()
    elif args[0] == 'fasta' and '-0' in args:
        # SAMTOOLS_FASTA: three reads per BAM, named after it, so pooled inputs are traceable.
        stem = Path(args[-1]).name
        Path(option('-0')).write_bytes(gzip.compress(''.join(f'>{stem}/{i}\nACGT\n' for i in range(3)).encode()))
    elif args[0] == 'fasta':
        for line in Path(args[-1]).read_text().splitlines():
            if not line.startswith('@'):
                print('>' + line.split('\t')[0] + '\nACGT')
elif name == 'iso-fastx' and args[0] == 'inspect':
    # State from the file name (.ccs., .mixed., .clustered., .subreads.), fl otherwise.
    fastx = Path(option('--fastx')).name
    state = next((s for s in ('ccs', 'mixed', 'clustered', 'subreads') if f'.{s}.' in fastx), 'fl')
    Path(option('--prefix') + '.inspect.tsv').write_text(f'file\treads\tstate\tkit\n{fastx}\t3\t{state}\tisoseqx\n')
elif name == 'iso-fastx' and args[0] == 'orient':
    shutil.copyfile(option('--fastx'), option('--output'))
    Path(option('--prefix') + '.orient.tsv').touch()
elif name == 'collapse' and args[0] == 'chain':
    shutil.copyfile(option('--bed').split(',')[0], option('--prefix') + '.models.bed')
    for suffix in ('.support.tsv', '.counts.tsv', '.members.tsv', '.excluded.bed'):
        Path(option('--prefix') + suffix).touch()
elif name == 'fxsplit':
    # Real chunk names carry the --suffix (meta.id), so every chunk FASTA has a distinct staged name.
    Path('chunks').mkdir()
    shutil.copyfile(option('-f'), f"chunks/tmp_chunk_0_{option('--suffix')}.fasta.gz")
    if os.environ.get('TEST_TWO_CHUNKS'):  # pooled fragment re-alignment needs several chunks per sample
        shutil.copyfile(option('-f'), f"chunks/tmp_chunk_1_{option('--suffix')}.fasta.gz")
elif name == 'minimap2':
    print(sam, end='')  # INFO: piped into samtools sort
elif name == 'iso-cigar':
    stem = option('--bam').removesuffix('.bam')
    Path(stem + '.extended.bam').touch()
    Path(stem + '.extended.bam.bai').touch()
elif name == 'iso-align' and '--realigned-bam' in args:
    assert Path(option('--bam')).is_file() and Path(option('--realigned-bam')).is_file(), args
    Path(option('--report')).touch()
    # singleton chunks accept nothing: their base BAM must skip the drop step and they emit no fragments BAM
    Path(option('--output')).write_text('' if '.singleton.' in option('--output') else 'r1\n')
elif name == 'iso-align':
    assert all(Path(path).is_file() for path in [option('--bam')] + option('--reads').split(',')), args
    Path(option('--report')).touch()
    # empty unless TEST_FRAGMENTS: the module deletes an empty FASTA, so no fragment re-alignment
    Path(option('--output')).write_text('>r1\nACGT\n' if os.environ.get('TEST_FRAGMENTS') else '')
elif name == 'iso-segment':
    Path('chr1@' + option('--prefix') + '.hq.bed').write_text('chr1\t0\t4\tr1\n')
elif name == 'iso-fusion':
    beds = option('--query').split(',')
    assert beds and all(Path(bed).is_file() for bed in beds), beds
    directory = Path(option('--prefix'))
    directory.mkdir()
    (directory / 'fusions.free.bed').write_text('chr1\t0\t4\tr1\n')
elif name in ('splicetransformer', 'spliceai') and args[0] == 'chunk':
    # one chunk per contig of genome.fa, chunk numbers restarting per contig (real chunker layout)
    Path('chunks').mkdir()
    gz = '--gzip' in args
    for contig in ('chrA', 'chrB'):
        Path(f"chunks/tmp.{contig}.chunk.0.fa{'.gz' if gz else ''}").touch()
elif name in ('splicetransformer', 'spliceai') and args[0] == 'predict':
    outdir, stem = Path(option('--outdir')), Path(option('--sequence')).name.split('.fa')[0]
    outdir.mkdir()
    for track in ('donor_plus', 'donor_minus', 'acceptor_plus', 'acceptor_minus'):
        (outdir / f'{stem}.{track}.wig').write_text(f'{stem}\n')
elif name == 'wigToBigWig':
    shutil.copyfile(args[-3], args[-1])
elif name == 'bigtools' and args[0] == 'bigwigmerge':
    Path(args[-1]).write_text(''.join(Path(bw).read_text() for bw in args[1:-1]))
elif name == 'splicing':
    tracks = sorted(path.name for path in Path(option('--bigwig-dir')).iterdir())
    assert len(tracks) == 4 and all(any(site in t for site in ('donor', 'acceptor')) and any(strand in t for strand in ('plus', 'minus'))
                                    for t in tracks), tracks
    Path(option('--prefix') + '.derived.tsv').touch()
else:
    raise AssertionError((name, args))
with open(os.environ['TEST_CALLS'], 'a') as log:
    log.write(json.dumps([name, args, os.getcwd()]) + '\n')
'''


with tempfile.TemporaryDirectory(prefix="ark-isoseq-") as temporary:
    temporary = Path(temporary)
    binary = temporary / "bin"
    binary.mkdir()
    for name in ("pbindex", "lima", "isoseq", "samtools", "iso-fastx", "collapse", "fxsplit", "minimap2", "iso-cigar",
                 "iso-align", "iso-segment", "iso-fusion", "splicetransformer", "spliceai", "wigToBigWig", "bigtools",
                 "splicing"):
        executable = binary / name
        executable.write_text(TOOLS)
        executable.chmod(0o755)
    # Cigar extension stages genome and annotation side by side, so they need distinct names.
    for name in ("genome.fa", "genome.mmi", "annotation.bed", "chrom.sizes", "SpTransformer_pytorch.ckpt"):
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
params.accept = false
params.pool = false
params.cluster_engine = 'isoseq'
params.test_bigwigs = null
params.test_algorithm = 'spliceai'
params.test_scores = null
params.test_weights = null
""")
    # The second pass runs only for the ark cases (--aligner ark); --cigar toggles cigar extension there.
    harness = temporary / "main.nf"
    harness.write_text("""
include { ISOSEQ } from 'REPO/src/subworkflows/isoseq/main.nf'
include { FASTX_PREPARE } from 'REPO/src/subworkflows/fastx/main.nf'
include { SPLIT_ALIGN_CLEAN_CHUNKS } from 'REPO/src/subworkflows/split_align/main.nf'
workflow {
    if (params.entrypoint == 'flnc') {
        // Same call as src/subworkflows/preprocessing/main.nf: state per file, lima/orient, pooling.
        FASTX_PREPARE(params.global_input_dir, null, params.flnc_input_state, params.cluster_mode, 'pooled')
        reads = FASTX_PREPARE.out.reads
        reads.view { meta, f -> 'POOL\\t' + meta.id + '\\t' + meta.sample_id + '\\t' + meta.singleton + '\\t' + f.size() }
    } else {
        ISOSEQ(params.global_input_dir, params.global_primers, 1, params.cluster_mode, params.cluster_engine,
            'pooled', params.entrypoint, params.entrypoint in ['refine', 'cluster'])
        reads = ISOSEQ.out.reads
    }
    SPLIT_ALIGN_CLEAN_CHUNKS(reads, Channel.value(file("${projectDir}/genome.fa")),
        Channel.value([[:], file("${projectDir}/genome.mmi")]), Channel.value([[:], file("${projectDir}/annotation.bed")]),
        Channel.value([[:], []]), params.aligner, false, false, false, params.cigar, params.aligner == 'ark', false,
        params.accept, params.pool, params.reconstruct_engine, Channel.empty())
    SPLIT_ALIGN_CLEAN_CHUNKS.out.reads.view { meta, bed -> 'RESULT\\t' + meta.id }
}
""".replace("REPO", str(ROOT)))
    environment = dict(os.environ, PATH=f"{binary}:{os.environ['PATH']}", NXF_OFFLINE="true", NXF_ANSI_LOG="false")
    base = [NEXTFLOW, "-C", f"{ROOT}/src/nextflow.config,{config}", "run"]

    isoseqx = [f"IsoSeqX_bc{i:02}_5p--IsoSeqX_3p" for i in (1, 2)]
    neb = ["NEB_5p--NEB_Clontech_3p", "NEB_5p--primer_3p"]
    # cigar None: mm2, second pass off. True/False: ark, second pass on, cigar extension on/off.
    # flnc honors the cluster mode through POOL_READS; "per_sample" only makes `expected` below the per-sample ids.
    # Engine (5th field; default isoseq): none aligns refined reads directly (SAMTOOLS_FASTA + POOL_READS).
    for entrypoint, mode, pairs, cigar, *engine in (
            ("refine", "per_sample", isoseqx, None), ("refine", "multi_sample", isoseqx, None),
            ("refine", "both", isoseqx, None), ("refine", "per_sample", isoseqx[:1], None),
            ("ccs", "per_sample", isoseqx[:1], None), ("ccs", "both", isoseqx, None),
            ("refine", "both", neb, None), ("ccs", "both", neb, None), ("cluster", "both", isoseqx, None),
            ("refine", "multi_sample", isoseqx, True), ("refine", "per_sample", isoseqx, True),
            ("flnc", "per_sample", isoseqx, True), ("flnc", "per_sample", isoseqx, False),
            ("flnc", "multi_sample", isoseqx, True), ("flnc", "both", isoseqx, True),
            ("cluster", "both", isoseqx, None, "none"), ("refine", "per_sample", isoseqx, None, "none"),
            ("ccs", "multi_sample", isoseqx, None, "none"), ("flnc", "per_sample", isoseqx, True, "isoseq", "accept"),
            ("flnc", "per_sample", isoseqx, True, "isoseq", "accept", "pool")):
        # 6th field "accept": iso-align selects one read per chunk and fragment acceptance is on;
        # 7th field "pool": two chunks per sample, fragment reads realigned once per sample
        accept = "accept" in engine[1:]
        pool = "pool" in engine[1:]
        engine = engine[0] if engine else "isoseq"
        case = temporary / f"{entrypoint}-{mode}-{pairs[0]}-{len(pairs)}-{cigar}-{engine}{'-accept' if accept else ''}{'-pool' if pool else ''}"
        inputs = case / "02_LIMA"
        inputs.mkdir(parents=True)
        samples = {f"movie.part.hifi.{pair}" for pair in pairs}
        files = {"refine": [f"movie.part.hifi_fl.{pair}.bam" for pair in pairs], "ccs": ["movie.part.hifi.bam"],
                 "cluster": [f"movie.part.hifi.{pair}_flnc.bam" for pair in pairs],
                 "flnc": [f"{sample}.{kind}.fasta.gz" for sample in samples for kind in ("hq", "singletons")]}
        file_bytes = {}
        for i, name in enumerate(files[entrypoint]):
            if name.endswith(".fastq.gz"):
                stem = name.removesuffix(".fastq.gz")
                reads = "" if "empty" in name else "".join(f"@{stem}/{n} desc\nACGT\n+\nIIII\n" for n in range(3))
                (inputs / name).write_bytes(gzip.compress(reads.encode()))
                if i == 0:
                    (inputs / f"{name}.pbi").touch()
            elif entrypoint == "flnc":
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
        environment.pop("TEST_FRAGMENTS", None)
        environment.pop("TEST_TWO_CHUNKS", None)
        if accept:
            environment["TEST_FRAGMENTS"] = "1"
        if pool:
            environment["TEST_TWO_CHUNKS"] = "1"
        command = base + [str(harness), "--global_input_dir", str(inputs),
                          "--global_output_dir", str(case / "results"), "--cluster_mode", mode,
                          "--entrypoint", entrypoint]
        command += ["--global_primers", str(primers)] if entrypoint != "cluster" else []  # cluster needs none
        command += ["--aligner", "ark", "--cigar", str(cigar).lower()] if cigar is not None else []
        command += ["--cluster_engine", engine]
        command += ["--accept", "true"] if accept else []
        command += ["--pool", "true"] if pool else []
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
        refine = [args for name, args, *_ in records if name == "isoseq" and args[0] == "refine"]
        cluster = [args for name, args, *_ in records if name == "isoseq" and args[0] == "cluster2"]
        assert len(refine) == (len(pairs) if entrypoint in ("refine", "ccs") else 0), refine
        assert len(cluster) == (len(expected) if entrypoint != "flnc" and engine == "isoseq" else 0), cluster
        assert sum(name == "lima" for name, *_ in records) == (entrypoint == "ccs")
        assert sum(name == "pbindex" for name, *_ in records) == (
            (entrypoint == "refine" or (entrypoint == "cluster" and engine == "isoseq")) and len(pairs) > 1)
        to_fasta = [args for name, args, *_ in records if name == "samtools" and args[0] == "fasta" and "-0" in args]
        # Engine none converts every refined BAM to FASTA once, also in 'both'.
        assert len(to_fasta) == (len(pairs) if engine == "none" else 0), to_fasta
        # collapse chain runs once per result group (the stub fusion detector emits free reads only);
        # cluster2 records are clustered, so they keep every chain (sensitive preset).
        chains = [args for name, args, *_ in records if name == "collapse" and args[0] == "chain"]
        assert len(chains) == len(expected), chains
        presets = {args[args.index("--preset") + 1] for args in chains}
        assert presets == ({"sensitive"} if entrypoint != "flnc" and engine == "isoseq" else {"balanced"}), presets
        assert not any(name == "samtools" and args[0] == "merge" for name, args, *_ in records), records
        # Every result id has one hq and one singleton chunk (only hq with engine none). With the second pass on, iso-align runs once per
        # chunk BAM with exactly the chunk FASTA minimap2 aligned it from (reads name -> BAM name).
        # INFO: minimap2 pipes into samtools sort in the same task directory: pair them by it
        sorted_bams = {cwd: args[args.index("-o") + 1] for name, args, cwd in records if name == "samtools" and args[0] == "sort"}
        is_fragments = lambda reads: ".fragments." in Path(reads).name
        chunks = {args[-1]: sorted_bams[cwd].removesuffix(".bam") for name, args, cwd in records
                  if name == "minimap2" and not is_fragments(args[-1])}
        found = sorted((args[args.index("--bam") + 1], args[args.index("--reads") + 1])
                       for name, args, *_ in records if name == "iso-align" and "--realigned-bam" not in args)
        suffix = ".extended.bam" if cigar else ".bam"
        assert len(chunks) == (1 if engine == "none" else 2) * len(expected) * (2 if pool else 1), chunks  # none: no singleton class
        assert found == (sorted((sam + suffix, reads) for reads, sam in chunks.items()) if cigar is not None else []), (found, chunks)
        # Acceptance: per chunk, iso-align compares the pass-1 BAM (not the iso-cigar BAM) with the fragments BAM;
        # samtools drops the accepted reads from the base (iso-cigar) BAM and keeps only them in the fragments BAM;
        # iso-segment reads exactly those two filtered BAMs per chunk.
        accepts = sorted((Path(args[args.index("--bam") + 1]).name, Path(args[args.index("--realigned-bam") + 1]).name,
                          ".singleton." not in args[args.index("--output") + 1])
                         for name, args, *_ in records if name == "iso-align" and "--realigned-bam" in args)
        fragment_bams = sorted(Path(sorted_bams[cwd]).name for name, args, cwd in records
                               if name == "minimap2" and is_fragments(args[-1]))
        drops = sorted(Path(args[-1]).name for name, args, *_ in records if name == "samtools" and "-U" in args)
        # each selection is re-sorted (iso-cigar BAMs are not strictly sorted) before indexing
        resorted = [args for name, args, *_ in records if name == "samtools" and args[0] == "sort" and args[-1] == "selected.unsorted.bam"]
        assert len(resorted) == (2 * len(drops) if accept else 0), resorted
        keeps = sorted(Path(args[-1]).name for name, args, *_ in records
                       if name == "samtools" and "-N" in args and "-U" not in args)
        selected = sorted(Path(args[args.index("-o") + 1]).name for name, args, *_ in records
                          if name == "samtools" and args[0] == "sort" and args[-1] == "selected.unsorted.bam"
                          and args[args.index("-o") + 1].endswith(".selected.bam"))
        segmented = sorted(Path(args[args.index("--bam") + 1]).name for name, args, *_ in records if name == "iso-segment")
        if accept:
            pass1 = sorted(Path(sam).name + ".bam" for sam in chunks.values())
            if pool:
                # one re-alignment per sample and read class, shared by both chunks of that sample
                assert len(fragment_bams) == len(chunks) // 2, (fragment_bams, chunks)
                assert all(".pooled" in f for f in fragment_bams), fragment_bams
                assert [p for p, *_ in accepts] == pass1 and sorted(set(f for _, f, _ in accepts)) == fragment_bams, (accepts, pass1)
            else:
                assert len(fragment_bams) == len(chunks), (fragment_bams, chunks)
                assert [p for p, *_ in accepts] == pass1 and sorted(f for _, f, _ in accepts) == fragment_bams, (accepts, pass1)
            # chunks with accepted reads: base BAM filtered, fragments BAM kept; chunks without: base BAM untouched, no fragments
            hit = [(p, f) for p, f, any_accepted in accepts if any_accepted]
            assert hit and len(hit) < len(accepts), accepts
            assert drops == sorted(p.removesuffix(".bam") + ".extended.bam" for p, _ in hit), drops
            assert keeps == sorted(f for _, f in hit), keeps
            # one kept fragments BAM per accepting chunk, named after the chunk (also from a pooled input)
            assert len(selected) == len(set(selected)) == len(hit), (selected, hit)
            assert segmented == sorted([d.removesuffix(".bam") + ".filtered.bam" for d in drops] +
                                       [p.removesuffix(".bam") + ".extended.bam" for p, _, a in accepts if not a] +
                                       selected), segmented
        else:
            assert not accepts and not fragment_bams and not drops and not keeps, (accepts, fragment_bams, drops, keeps)
        second_pass = "" if cigar is None else f", ark second pass, cigar extension {'on' if cigar else 'off'}"
        second_pass += ", fragment acceptance" if accept else ""
        second_pass += ", pooled per sample" if pool else ""
        print(f"PASS {entrypoint}: {mode}, {engine}, {', '.join(pairs)}{second_pass}", flush=True)

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

    # flnc routing: ccs -> lima, mixed -> orient, clustered -> sensitive preset, fl as is; subreads stop the run.
    case = temporary / "flnc-routing"
    case.mkdir()
    for name in ("s1.ccs.fastq.gz", "s2.mixed.fastq.gz", "s3.clustered.fastq.gz", "s4.fastq.gz"):
        (case / name).write_bytes(gzip.compress(f"@{name}/0\nACGT\n+\nIIII\n".encode()))
    calls = case / "calls.jsonl"
    environment.update(TEST_PRIMER_PAIRS=json.dumps(isoseqx[:1]), TEST_CALLS=str(calls))
    result = subprocess.run(base + [str(harness), "--entrypoint", "flnc", "--global_input_dir", str(case),
                                   "--global_output_dir", str(case / "results"), "--cluster_mode", "per_sample"],
                            cwd=case, env=environment, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=120)
    assert result.returncode == 0, result.stdout
    records = [json.loads(line) for line in calls.read_text().splitlines()]
    assert [Path(a[0]).name for n, a, *_ in records if n == "lima"] == ["s1.ccs.fastq.ccs_names.fastq.gz"], records
    assert [Path(a[a.index("--fastx") + 1]).name for n, a, *_ in records if n == "iso-fastx" and a[0] == "orient"] == \
        ["s2.mixed.fastq.gz"], records
    presets = {a[a.index("--prefix") + 1].split(".chr1")[0]: a[a.index("--preset") + 1]
               for n, a, *_ in records if n == "collapse"}
    assert presets == {"s1.ccs": "balanced", "s2.mixed": "balanced", "s3.clustered": "sensitive", "s4": "balanced"}, presets
    (case / "s5.subreads.fastq.gz").write_bytes(gzip.compress(b"@s5/0\nACGT\n+\nIIII\n"))
    result = subprocess.run(base + [str(harness), "--entrypoint", "flnc", "--global_input_dir", str(case),
                                   "--global_output_dir", str(case / "results")],
                            cwd=case, env=environment, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=120)
    assert result.returncode != 0 and "looks like 'subreads'" in result.stdout, result.stdout
    print("PASS flnc routing: ccs -> lima, mixed -> orient, clustered -> sensitive, subreads rejected", flush=True)

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

    for flag, value, diagnostic in (("--cluster_engine", "foo", "Unknown cluster_engine option"),
                                    ("--cluster_engine", "cdhit", "removed in v2.1.0"),
                                    ("--reconstruct_engine", "isoquant", "not implemented in v2.1.0"),
                                    ("--isoseq_cluster2_mode", "both", "isoseq_cluster2_mode was renamed to cluster_mode")):
        result = subprocess.run(base + [str(ROOT / "src/main.nf"), "--entrypoint", "cluster", flag, value],
                                cwd=temporary, env=environment, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=120)
        assert result.returncode != 0 and diagnostic in result.stdout, result.stdout
    print("PASS --cluster_engine/--reconstruct_engine validation and isoseq_cluster2_mode rename", flush=True)

    # Splice scores: SPLICING runs splice_model (SpliceAI or SpliceTransformer) or reads its bigWig dir, and the
    # same splicing crate (SPLICEAI_DERIVE, region floor on) calibrates either model's four bigWigs.
    splicing = temporary / "splicing.nf"
    splicing.write_text("""
include { SPLICING } from 'REPO/src/subworkflows/splicing/main.nf'
workflow {
    SPLICING(Channel.value(file("${projectDir}/genome.fa")), Channel.value([[:], file("${projectDir}/annotation.bed")]),
        params.test_bigwigs, params.test_algorithm, params.splice_model, params.test_scores, null,
        Channel.value(file("${projectDir}/chrom.sizes")), params.spliceai_chunk_compression, params.test_weights, Channel.empty())
    SPLICING.out.scores.view { meta, scores -> 'SCORES\\t' + scores.name }
    SPLICING.out.bigwigs.view { meta, dir -> 'BIGWIGS\\t' + (meta.id ?: '') + '\\t' + (dir ? dir.name : '') }
}
""".replace("REPO", str(ROOT)))
    tracks = ["acceptor_minus.bw", "acceptor_plus.bw", "donor_minus.bw", "donor_plus.bw"]
    given = temporary / "given_bigwigs"
    given.mkdir()
    for track in tracks:  # benchmark naming of the SpliceTransformer bigWigs
        (given / track).touch()
    minisplice_scores = temporary / "minisplice.tsv"
    minisplice_scores.touch()
    # model, extra flags, model tool run, bigWig dir id/name, derived scores (None: no derive)
    for model, extra, tool, bigwigs, scores in (
            ("splicetransformer", [], "splicetransformer", ("splicetransformer", "splicetransformer"), "genome.derived.tsv"),
            ("splicetransformer", ["--spliceai_chunk_compression", "true", "--test_weights", str(temporary / "SpTransformer_pytorch.ckpt")],
             "splicetransformer", ("splicetransformer", "splicetransformer"), "genome.derived.tsv"),
            ("splicetransformer", ["--test_bigwigs", str(given)], None, ("splicetransformer", "given_bigwigs"), "genome.derived.tsv"),
            ("spliceai", [], "spliceai", ("spliceai", "spliceai"), "genome.derived.tsv"),
            ("spliceai", ["--test_algorithm", "minisplice", "--test_scores", str(minisplice_scores)], None, ("", ""), None)):
        case = temporary / f"splicing-{model}-{len(extra)}"
        case.mkdir()
        calls = case / "calls.jsonl"
        environment.update(TEST_CALLS=str(calls))
        result = subprocess.run(base + [str(splicing), "--splice_model", model, "--minimap2_align_use_splice_scores", "true",
                                       "--global_output_dir", str(case / "results")] + extra,
                                cwd=case, env=environment, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=120)
        assert result.returncode == 0, result.stdout
        lines = result.stdout.splitlines()
        assert [l.split("\t")[1] for l in lines if l.startswith("SCORES\t")] == [scores or "minisplice.tsv"], result.stdout
        assert [tuple(l.split("\t")[1:]) for l in lines if l.startswith("BIGWIGS\t")] == [bigwigs], result.stdout
        records = [json.loads(line) for line in calls.read_text().splitlines()] if calls.exists() else []
        runs = [(name, args) for name, args, *_ in records if name in ("splicetransformer", "spliceai")]
        assert {name for name, _ in runs} == ({tool} if tool else set()), runs
        if tool:
            chunk = [args for _, args in runs if args[0] == "chunk"]
            predict = [args for _, args in runs if args[0] == "predict"]
            flank = "4000" if tool == "splicetransformer" else "50000"
            assert len(chunk) == 1 and chunk[0][chunk[0].index("--flank-size") + 1] == flank, chunk
            # one predict per contig chunk, flank trimmed with the chunker's flank size
            assert sorted(Path(a[a.index("--sequence") + 1]).name for a in predict) == \
                sorted(f"tmp.{c}.chunk.0.fa{'.gz' if '--gzip' in chunk[0] else ''}" for c in ("chrA", "chrB")), predict
            assert all(a[a.index("--offset") + 1] == flank for a in predict), predict
            assert ("--gzip" in chunk[0]) == ("--spliceai_chunk_compression" in extra and tool == "splicetransformer"), chunk
            assert all(("--model" in a) == ("--test_weights" in extra) for a in predict), predict
            assert sum(name == "wigToBigWig" for name, *_ in records) == 8, records
            merged = sorted(Path(a[-1]).name for name, a, *_ in records if name == "bigtools")
            assert merged == tracks, merged
            published = sorted(p.name for p in (case / "results" / f"00_{tool.upper()}" / tool).iterdir())
            assert published == sorted(tracks), published
        derive = [args for name, args, *_ in records if name == "splicing"]
        assert len(derive) == (1 if scores else 0), derive
        if scores:
            # TOGA2 CDS region floor, unchanged for either model
            assert " ".join(derive[0]).startswith("--include-ss-from-regions --position-for-ss-regions cds --region-floor 11"), derive
            assert Path(derive[0][derive[0].index("--bigwig-dir") + 1]).name == bigwigs[1], derive
        print(f"PASS splice scores: splice_model {model}{' ' + ' '.join(Path(e).name for e in extra) if extra else ''}", flush=True)

    result = subprocess.run(base + [str(ROOT / "src/main.nf"), "--entrypoint", "cluster", "--splice_model", "foo"],
                            cwd=temporary, env=environment, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=120)
    assert result.returncode != 0 and "Unknown splice_model option" in result.stdout, result.stdout
    print("PASS --splice_model validation", flush=True)
