"""Small, reproducible synthetic isolate with a low-abundance second segment."""
import random
from pathlib import Path
import pysam


def make_demo(outdir):
    outdir = Path(outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Demo output directory must be empty")
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(2026)
    seqs = ["".join(rng.choices("ACGT", k=n)) for n in (2400, 1600)]
    with (outdir / "reference.fa").open("w") as f:
        for i,s in enumerate(seqs):
            f.write(f">segment_{i+1}\n{s}\n")
    header = {"HD": {"VN":"1.6"}, "SQ":[{"SN":f"segment_{i+1}","LN":len(s)} for i,s in enumerate(seqs)]}
    with pysam.AlignmentFile(outdir / "unsorted.bam", "wb", header=header) as bam, \
         (outdir / "reads.fastq").open("w") as fq:
        for i in range(640):
            target = 0 if i < 560 else 1
            start = rng.randrange(len(seqs[target])-100+1)
            sequence = seqs[target][start:start+100]
            read = pysam.AlignedSegment()
            read.query_name = f"fragment_{i}"
            read.query_sequence = sequence
            read.query_qualities = pysam.qualitystring_to_array("I"*100)
            read.reference_id = target
            read.reference_start = start
            read.mapping_quality = 60
            read.cigarstring = "100M"
            bam.write(read)
            fq.write(f"@{read.query_name}\n{sequence}\n+\n{'I'*100}\n")
    pysam.sort("-o", str(outdir / "sample.bam"), str(outdir / "unsorted.bam"))
    pysam.index(str(outdir / "sample.bam"))
    (outdir / "unsorted.bam").unlink()
    (outdir / "targets.tsv").write_text("contig\tgroup\nsegment_1\tvirus\nsegment_2\tvirus\n")
    print(outdir)
