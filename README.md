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

## Check the installation

```bash
cargo run -p mito-cli -- doctor
bash scripts/verify.sh
```

The verification script runs the native, Rust, web, and fixture checks.
For source navigation: `core/` contains the C++ analysis engine, `ffi/`
connects it to Rust, `cli/` and `server/` expose it, and `web/` with
`visualization-lib/` renders the results.
