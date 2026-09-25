"""Strict streaming FASTQ IO and reproducible, pair-preserving sampling."""
import gzip
import hashlib
import json
import sqlite3
import tempfile
from itertools import zip_longest
from functools import lru_cache
from pathlib import Path
from .progress import Progress


def fragment_name(name):
    name = name.split()[0]
    return name[:-2] if name.endswith(("/1", "/2")) else name


@lru_cache(maxsize=100000)
def priority(name, seed, replicate=0):
    if isinstance(name, tuple):
        read_group, read_name = name
        normalized_name = fragment_name(read_name)
        identity = (json.dumps([read_group, normalized_name], separators=(",", ":"))
                    if read_group else normalized_name)
    else:
        identity = fragment_name(name)
    key = f"{seed}:{replicate}:{identity}".encode()
    return int.from_bytes(hashlib.blake2b(key, digest_size=8).digest(), "big") / 2**64


def fastq(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="ascii") as handle:
        index = 0
        while True:
            header = handle.readline()
            if not header:
                break
            index += 1
            seq, plus, qual = [handle.readline().rstrip("\r\n") for _ in range(3)]
            if not header.startswith("@") or not plus.startswith("+") or not seq or len(seq) != len(qual):
                raise ValueError(f"Malformed four-line FASTQ record {index} in {path}")
            if any(ord(c) < 33 or ord(c) > 126 for c in qual):
                raise ValueError(f"Invalid Phred+33 quality at record {index} in {path}")
            yield header[1:].strip(), seq.upper(), qual


def fragments(r1, r2=None):
    if r2 is None:
        for read in fastq(r1):
            yield fragment_name(read[0]), [read]
    else:
        for a, b in zip_longest(fastq(r1), fastq(r2)):
            if a is None or b is None or fragment_name(a[0]) != fragment_name(b[0]):
                raise ValueError("Paired FASTQs must have matching names and record counts")
            yield fragment_name(a[0]), [a, b]


def fingerprint(path):
    p = Path(path).resolve()
    stat = p.stat()
    return {"path": str(p), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def bam_record_total(bam):
    """Index-derived alignment count; no extra scan just to obtain progress."""
    try:
        # CRAI indexes do not provide reliable counts; BAM BAI/CSI indexes do.
        return bam.mapped + bam.unmapped if bam.is_bam and bam.has_index() else None
    except (ValueError, OSError):
        return None


def bam_fragments(path, tempdir, progress=None):
    """Group primary records on disk; recover original orientation, including unmapped reads."""
    import pysam
    progress = progress or Progress()
    with tempfile.TemporaryDirectory(prefix="bam-reads-", dir=tempdir) as work:
        conn = sqlite3.connect(str(Path(work) / "reads.sqlite"))
        try:
            conn.execute("CREATE TABLE reads (rg TEXT, name TEXT, mate INTEGER, seq TEXT, qual TEXT, "
                         "PRIMARY KEY(rg,name,mate))")
            with pysam.AlignmentFile(path, "rb") as bam:
                progress.phase("Preparing BAM reads for k-mers", bam_record_total(bam), "records")
                for read in bam.fetch(until_eof=True):
                    progress.advance()
                    if read.is_secondary or read.is_supplementary:
                        continue
                    seq, qual = read.get_forward_sequence(), read.get_forward_qualities()
                    if not read.query_name or seq is None or qual is None:
                        raise ValueError("BAM primary read lacks name, sequence or qualities; supply FASTQ with --kmer-source fastq, or use --kmer-source none")
                    if any(op == 5 for op, length in (read.cigartuples or [])):
                        raise ValueError("BAM primary read is hard-clipped; full sequence unavailable. Use --kmer-source fastq or none")
                    rg = str(read.get_tag("RG")) if read.has_tag("RG") else ""
                    name = fragment_name(read.query_name)
                    mate = 1 if read.is_read1 else 2 if read.is_read2 else 0
                    try:
                        conn.execute("INSERT INTO reads VALUES (?,?,?,?,?)", (rg, name, mate, seq.upper(), "".join(chr(q+33) for q in qual)))
                    except sqlite3.IntegrityError as exc:
                        raise ValueError(f"Repeated primary BAM read/mate: {name}") from exc
            conn.commit()
            progress.phase("Counting k-mer reads", conn.execute("SELECT count(*) FROM reads").fetchone()[0], "reads")
            previous, group = None, []
            for rg, name, seq, qual in conn.execute("SELECT rg,name,seq,qual FROM reads ORDER BY rg,name,mate"):
                identity = (rg, name)
                if previous is not None and identity != previous:
                    yield previous, group
                    group = []
                group.append((name, seq, qual))
                previous = identity
            if group:
                yield previous, group
        finally:
            conn.close()


def compare_bam_fastq(bam_path, r1, r2, progress=None):
    """Order-independent record fingerprints of names, original sequences and qualities."""
    import pysam
    progress = progress or Progress()
    def signature(records):
        count = digest = 0
        for name, seq, quality in records:
            payload = json.dumps([fragment_name(name), seq.upper(), quality], separators=(",", ":"))
            digest = (digest + int.from_bytes(hashlib.sha256(payload.encode()).digest(), "big")) % (1 << 256)
            count += 1
        return {"records": count, "fingerprint": f"{digest:064x}"}
    def bam_records():
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            progress.phase("Checking BAM identity", bam_record_total(bam), "records")
            for read in bam.fetch(until_eof=True):
                progress.advance()
                if read.is_secondary or read.is_supplementary:
                    continue
                quality = read.get_forward_qualities()
                yield read.query_name or "*", read.get_forward_sequence() or "*", None if quality is None else "".join(chr(q+33) for q in quality)
    a = signature(bam_records())
    progress.phase("Checking FASTQ identity", unit="reads")
    def fastq_records():
        for _, reads in fragments(r1, r2):
            for read in reads:
                progress.advance()
                yield read
    b = signature(fastq_records())
    return {"bam": a, "fastq": b, "match": a == b,
            "method": "Order-independent SHA256 multiset fingerprint of normalized names, forward sequences and qualities; excludes secondary/supplementary BAM records"}
