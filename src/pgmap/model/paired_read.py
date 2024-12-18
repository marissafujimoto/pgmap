import typing


class PairedRead(typing.NamedTuple):
    """
    A single paired read containing a candidate for counting.

    Attributes:
        gRNA1_candidate (str): The candidate gRNA1.
        gRNA2_candidate (str): The candidate gRNA2.
        barcode_candidate (str): The candidate barcode.
    """
    gRNA1_candidate: str
    gRNA2_candidate: str
    barcode_candidate: str
