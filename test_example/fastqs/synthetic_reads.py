"""Helpers for building the synthetic FASTQs of the test example.

The plates of the test example are built rather than sampled, so that every read parses and
the counts the pipeline finds are exactly the counts planted. These are the pieces both
generators need: the flanks a read is wrapped in, how a barcode is turned into a read, and
how the reads of a well are written to a reproducible gzip file.

The flanks are those of ``illumina_barcode_parser_params`` in ``config.yml``, with the
``upstream2`` the synthetic plates are given there.

"""

import gzip

# The flanks as they appear in a read, which is the reverse complement of barcode,
# upstream, upstream2, as `bc_orientation` is R2.
UPSTREAM = "CTCCCTACAATGTCGGATTTGTATTTAATAG"
UPSTREAM2 = "GCTACA"

# the quality string of a planted read, high enough everywhere that none is filtered
READ_QUALITY = "I"


def revcomp(seq):
    """Reverse complement of `seq`."""
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def hamming(s1, s2):
    """Hamming distance between equal-length strings."""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2, strict=True))


def mutate(barcode, position):
    """`barcode` with the nucleotide at `position` changed."""
    new = "C" if barcode[position] != "C" else "G"
    return barcode[:position] + new + barcode[position + 1 :]


def random_barcode(rng, length, far_from, min_distance, prefix=""):
    """A random barcode starting `prefix`, at least `min_distance` from all of `far_from`."""
    while True:
        barcode = prefix + "".join(
            rng.choice("ACGT") for _ in range(length - len(prefix))
        )
        if all(hamming(barcode, other) >= min_distance for other in far_from):
            return barcode


def read_of(barcode, outer=UPSTREAM2):
    """The sequence of a read of `barcode`, which parses unless `outer` is wrong."""
    return revcomp(barcode) + revcomp(UPSTREAM) + revcomp(outer)


def split_evenly(total, n):
    """Split `total` reads as evenly as possible into `n` parts."""
    return [total // n + (1 if i < total % n else 0) for i in range(n)]


def write_fastq(path, name, sequences):
    """Write `sequences` as the reads of one well, reproducibly.

    `mtime=0` so that the gzip header, and so the file, is the same on every run: the
    generators are meant to be rerunnable without changing their output.

    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with (
        open(path, "wb") as raw,
        gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as f,
    ):
        for i, seq in enumerate(sequences):
            f.write(f"@{name}:{i}\n{seq}\n+\n{READ_QUALITY * len(seq)}\n".encode())
