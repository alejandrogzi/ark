#!/usr/bin/env python3
"""Exercise the FXSPLIT errorStrategy from src/nextflow.config with stub fxsplit.

Usage: python3 assets/ci/test_fxsplit.py [path/to/nextflow]

Probed real-tool behavior (ghcr.io/alejandrogzi/fxsplit:latest):
- empty input -> exit 1 ("No FASTA records found", no chunks): benign, run must continue
- OOM-kill -> exit 137: must retry, then fail the run loudly (never silently starve downstream)

The test config below deliberately sets NO errorStrategy, so the repo's
'.*:FXSPLIT' closure is what actually runs in both cases.
"""

import gzip
import os
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[2]
NEXTFLOW = sys.argv[1] if len(sys.argv) > 1 else "nextflow"
TOOLS = r'''#!/usr/bin/env python3
import os, shutil, sys
from pathlib import Path
args = sys.argv[1:]
if '--version' in args:
    print('fxsplit test')
    sys.exit(0)
def option(flag):
    return args[args.index(flag) + 1]
code = int(os.environ.get('TEST_FXSPLIT_EXIT', '0'))
if code == 1:
    print('ERROR [fxsplit] ERROR: No FASTA records found')
    sys.exit(1)
if code == 137:
    sys.exit(137)
Path('chunks').mkdir()
shutil.copyfile(option('-f'), f"chunks/tmp_chunk_0_{option('--suffix')}.fasta.gz")
'''


with tempfile.TemporaryDirectory(prefix="ark-fxsplit-") as temporary:
    temporary = Path(temporary)
    binary = temporary / "bin"
    binary.mkdir()
    executable = binary / "fxsplit"
    executable.write_text(TOOLS)
    executable.chmod(0o755)

    reads = temporary / "reads.fasta.gz"
    reads.write_bytes(gzip.compress(b">r1\nACGT\n"))

    config = temporary / "test.config"
    config.write_text("""
process {
    withName: '.*' {
        cpus = 1
        memory = '256 MB'
        time = '1 min'
    }
}
docker.enabled = false
apptainer.enabled = false
singularity.enabled = false
conda.enabled = false
""")
    harness = temporary / "main.nf"
    harness.write_text("""
include { FXSPLIT } from 'REPO/src/modules/custom/fxsplit/main.nf'
// INFO: processes run inside a named workflow so the task is INNER:FXSPLIT and the
// INFO: repo's '.*:FXSPLIT' config selectors (ext.*, errorStrategy) apply as in the pipeline.
workflow INNER {
    take:
    input
    main:
    FXSPLIT(input)
    emit:
    out = FXSPLIT.out.fastx_gz
}
workflow {
    reads = Channel.of([ [id: 'sample', single_end: true, singleton: false], file(params.test_fasta) ])
    INNER(reads)
    INNER.out.view { meta, f -> 'CHUNK\\t' + meta.id }
}
""".replace("REPO", str(ROOT)))
    environment = dict(os.environ, PATH=f"{binary}:{os.environ['PATH']}", NXF_OFFLINE="true", NXF_ANSI_LOG="false")
    base = [NEXTFLOW, "-C", f"{ROOT}/src/nextflow.config,{config}", "run", str(harness),
            "--test_fasta", str(reads)]

    # Benign: empty-input exit 1 is ignored, the run completes with no chunks.
    environment.update(TEST_FXSPLIT_EXIT="1")
    result = subprocess.run(base, cwd=temporary, env=environment, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=120)
    assert result.returncode == 0, result.stdout
    assert not [line for line in result.stdout.splitlines() if line.startswith("CHUNK\t")], result.stdout
    print("PASS fxsplit exit 1 (empty input) ignored, run completes", flush=True)

    # Fatal: OOM-kill exit 137 retries, then fails the run instead of starving downstream.
    environment.update(TEST_FXSPLIT_EXIT="137")
    result = subprocess.run(base, cwd=temporary, env=environment, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=300)
    assert result.returncode != 0, result.stdout
    print("PASS fxsplit exit 137 (OOM) fails the run loudly", flush=True)
