import typing


class PairedRead(typing.NamedTuple):
    R1_candidate: str
    R2_candidate: str
    barcode_candidate: str
