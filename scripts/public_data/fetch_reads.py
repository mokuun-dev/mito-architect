#!/usr/bin/env python3
"""Download a bounded prefix of the pinned ENA run; no SRA prefetch/cache."""
import argparse
import hashlib
import json
import shutil
import tempfile
import time
import urllib.request
import zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def fastq_prefix(chunks, output, max_reads, max_decoded_bytes):
    """Validate complete records from a possibly truncated gzip prefix.

    The prefix is intentionally not described as CRC/whole-file MD5 verified.
    The output SHA256 identifies exactly the selected raw FASTQ bytes.
    """
    decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
    pending = b""
    record = []
    decoded = count = 0
    digest = hashlib.sha256()
    for chunk in chunks:
        while chunk:
            if decoder.eof:
                decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
            block = decoder.decompress(chunk, max_decoded_bytes - decoded + 1)
            decoded += len(block)
            if decoded > max_decoded_bytes:
                raise ValueError("decompressed byte limit exceeded")
            chunk = decoder.unused_data
            pending += block
            while b"\n" in pending:
                line, pending = pending.split(b"\n", 1)
                record.append(line)
                if len(record) == 4:
                    header, sequence, plus, quality = record
                    if (not header.startswith(b"@") or not plus.startswith(b"+")
                            or not sequence or len(sequence) != len(quality)
                            or any(c < 33 or c > 126 for c in quality)):
                        raise ValueError("malformed upstream FASTQ record")
                    data = b"\n".join(record) + b"\n"
                    output.write(data)
                    digest.update(data)
                    count += 1
                    record = []
                    if count == max_reads:
                        return {"selected_reads": count, "fastq_sha256": digest.hexdigest()}
            if len(pending) > 4 * 1024 * 1024:
                raise ValueError("FASTQ line exceeds 4 MiB limit")
    raise ValueError(f"download ended before {max_reads} complete records (got {count})")


def fetch(args):
    manifest = json.loads((ROOT / "fixtures/public/SRR18110025.json").read_text())
    if args.output.exists():
        raise ValueError(f"output already exists: {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".public-reads-", dir=args.output.parent))
    started = time.monotonic()
    received = 0
    download_hash = hashlib.sha256()
    request = urllib.request.Request(manifest["fastq_https"], headers={
        "Range": f"bytes=0-{args.max_download_bytes - 1}",
        "User-Agent": "Mito-Architect/public-prefix-1.0",
        "Accept-Encoding": "identity",
    })
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            if response.status not in (200, 206):
                raise ValueError(f"unexpected HTTP status {response.status}")
            if response.status == 206 and not response.headers.get("Content-Range", "").startswith("bytes 0-"):
                raise ValueError("server did not return the requested zero-offset prefix")

            def chunks():
                nonlocal received
                while received < args.max_download_bytes:
                    if time.monotonic() - started > args.timeout:
                        raise ValueError("download time limit exceeded")
                    chunk = response.read(min(65536, args.max_download_bytes - received))
                    if not chunk:
                        break
                    received += len(chunk)
                    download_hash.update(chunk)
                    yield chunk

            with (stage / "SRR18110025.fastq").open("wb") as output:
                selection = fastq_prefix(chunks(), output, args.max_reads, args.max_decoded_bytes)
            selection.update({"http_status": response.status,
                              "content_range": response.headers.get("Content-Range")})
        selection.update({"source": manifest, "requested_reads": args.max_reads,
                          "downloaded_prefix_bytes": received,
                          "downloaded_prefix_sha256": download_hash.hexdigest(),
                          "max_download_bytes": args.max_download_bytes,
                          "max_decoded_bytes": args.max_decoded_bytes,
                          "full_upstream_md5_verified": False})
        (stage / "selection.json").write_text(json.dumps(selection, indent=2) + "\n")
        stage.rename(args.output)
        print(json.dumps(selection, indent=2))
    finally:
        if stage.exists():
            shutil.rmtree(stage)


def positive(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-reads", type=positive, default=200)
    parser.add_argument("--max-download-bytes", type=positive, default=8 * 1024 * 1024)
    parser.add_argument("--max-decoded-bytes", type=positive, default=64 * 1024 * 1024)
    parser.add_argument("--timeout", type=positive, default=120)
    args = parser.parse_args()
    try:
        fetch(args)
    except (OSError, ValueError, zlib.error) as error:
        parser.exit(1, f"public data acquisition failed: {error}\n")


if __name__ == "__main__":
    main()
