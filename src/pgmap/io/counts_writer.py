import csv
from typing import Counter, Optional, IO, Iterator
from contextlib import contextmanager
import sys


def write_counts(output_path: Optional[str], counts: Counter[tuple[str, str, str]], barcodes: dict[str, str], id_mappings: dict[str, str]):
    """
    Write counts to an output file tsv with a header [id, seq_1, seq2, sample_id_1, ..., sample_id_n]. If no
    output_path is specified, write to STDOUT.

    Args
        output_path (Optional[str]): The output path to write to. If none, then the counts will be written to STDOUT.
        counts (Counter[tuple[str, str, str]]): The counts for each (gRNA1, gRNA2, barcode).
        id_mappings (dict[str, str]): A mapping from (gRNA1, gRNA2) to the paired guide id.
    """
    # TODO delete file if it already exists? manually test to see what happens if you don't do this

    with _open_file_or_stdout(output_path) as f:
        writer = csv.writer(f, delimiter='\t')

        sample_ids = list(barcodes.values())
        sample_id_to_barcode = {v: k for k, v in barcodes.items()}

        writer.writerow(["id", "seq_1", "seq_2"] + sample_ids)

        for (gRNA1, gRNA2), pg_id in id_mappings.items():
            row = [pg_id, gRNA1, gRNA2] + \
                [counts[(gRNA1, gRNA2, sample_id_to_barcode[sample_id])]
                 for sample_id in sample_ids]

            writer.writerow(row)


@contextmanager
def _open_file_or_stdout(path) -> Iterator[IO]:
    if path is None:
        yield sys.stdout
    else:
        with open(path, 'w', newline='') as f:
            yield f
