import typing


class PairedRead(typing.NamedTuple):
    gRNA1_candidate: str
    gRNA2_candidate: str
    barcode_candidate: str
