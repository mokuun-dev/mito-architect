# Mito-Architect

Mito-Architect is a research workbench for exploring human mitochondrial DNA
(mtDNA) in long-read sequencing data. It reads aligned SAM, BAM, or CRAM files
and presents coverage, sequence variants, structural changes, haplogroup
candidates, and supporting reads in a web interface. Results are for research
use; they are not clinical diagnoses.

## Install

On Debian/Ubuntu or Arch Linux, check the required tools and install missing
dependencies:

```bash
bash scripts/bootstrap.sh --check
bash scripts/bootstrap.sh --install
npm ci
```

The build needs C++20, CMake, Rust, Node.js, and htslib for BAM/CRAM input.
The bootstrap script prints any missing system tools.

## Run the web interface

Start the API and web application in separate terminals:

```bash
cargo run --release -p mito-server
```

```bash
npm --workspace web run dev
```

Open <http://127.0.0.1:5173> and upload an aligned SAM/BAM/CRAM file. The
default analysis uses result schema 0.5. The extended schema 0.6 evidence
graph is optional and may reach its resource limit on noisy long reads.
The local API listens on <http://127.0.0.1:8080>.

Raw FASTQ reads have no reference coordinates. Align them before upload if
you need variant, coverage, or structural-change results. For NUMT
specificity assessment, align against a reference that contains both nuclear
and mitochondrial sequences. Alignment only to the bundled mtDNA reference
leaves that assessment unavailable.

## Try public data

This command downloads a bounded selection of the public ONT run
`SRR18110025` and aligns it to the bundled mtDNA reference:

```bash
bash scripts/fetch_public_fixture.sh
```

Upload `.data/public/SRR18110025/alignment/competitive.bam` in the web
interface, or save a command-line result:

```bash
cargo run --release -p mito-cli -- analyze \
  --input .data/public/SRR18110025/alignment/competitive.bam \
  --json > .data/public/SRR18110025/result.json
```

The default download selects only the first 200 complete reads. It is a
software demonstration, not a representative biological sample or a
validated variant callset. The source and selection settings are recorded in
`fixtures/public/SRR18110025.json` and the generated files.

## Resource limits and input safety

Compressed BAM and CRAM can decode to far more data than their file size
suggests, so the engine caps decoded input and fails closed with
`MITO-E1601` when a cap is exceeded:

| Setting | CLI flag | Server variable | Default (CLI / server) |
| --- | --- | --- | --- |
| Decoded records | `--max-input-records` | `MITO_MAX_INPUT_RECORDS` | 50,000,000 / 2,000,000 |
| Decoded bases | `--max-input-bases` | `MITO_MAX_INPUT_BASES` | 10 Gbp / 1 Gbp |

CRAM is decoded only against the analysis reference (the bundled rCRS
unless `--reference` is given). A contig matches it by MD5, so `chrM`, `MT`
and `NC_012920.1` encodings of rCRS all decode. The engine ignores the
reference locations named inside a CRAM header (`UR` paths or URLs, `M5`
lookups through `REF_PATH`/`REF_CACHE`), because an uploaded file could
otherwise make the server read arbitrary local paths or contact remote
hosts. Reads on other contigs then fail to decode. For trusted local files,
`mito-cli analyze --allow-external-cram-references` restores htslib's
standard reference lookup; alternatively convert the CRAM to BAM.

The server keeps finished results in memory for at most 24 hours and at most
`MITO_MAX_STORED_JOBS` jobs (default 64) holding `MITO_MAX_STORED_BYTES`
(default 2 GiB). The oldest finished jobs are evicted first. Uploaded files
are deleted as soon as their analysis ends.

Without `MITO_API_KEY` the API only listens on loopback. It then accepts
requests only for loopback host names (add others with
`MITO_ALLOWED_HOSTS`), which blocks DNS-rebinding attacks. Browser requests
from any origin outside `MITO_CORS_ORIGINS` are refused, so other web pages
cannot submit uploads.

## Docker

```bash
docker compose up --build
```

The `server` service builds the API image and publishes it on
`127.0.0.1:8080`; the `web` service runs the development web server on
<http://127.0.0.1:5173>. A container has to bind `0.0.0.0`, so the compose
file sets `MITO_ALLOW_UNAUTHENTICATED=1` to run without an API key. Keep the
port on host loopback in that mode, or set `MITO_API_KEY` instead.

Installed or relocated binaries find the bundled reference and annotation
files through `MITO_DATA_DIR` (the image sets it to `/app/core/data`). The
native core needs htslib 1.18 or newer.

## Releases

CI (`.github/workflows/ci.yml`) runs `scripts/verify.sh`, an MSRV check,
ASan/UBSan and TSan test builds, a container smoke test
(`scripts/verify_container.sh`), and a dependency audit. Pushing a tag that
matches the crate version, for example `v0.5.0`, runs the same checks and
then publishes `ghcr.io/<owner>/mito-architect` and a GitHub release with
Linux binaries, the static web build, and `SHA256SUMS`.

## Check the installation

```bash
cargo run -p mito-cli -- doctor
bash scripts/verify.sh
```

The verification script runs the native, Rust, web, and fixture checks.
For source navigation: `core/` contains the C++ analysis engine, `ffi/`
connects it to Rust, `cli/` and `server/` expose it, and `web/` with
`visualization-lib/` renders the results.
