
import typing

from pgmap.model.trim_coordinate import TrimCoordinate


class TrimStrategy(typing.NamedTuple):
    """
    A container describing a trimming stategy for paired guide sequencing.

    Attributes:
        gRNA1 (TrimCoordinate): The trim coordinate for the gRNA1.
        gRNA2 (TrimCoordinate): The trim coordinate for the gRNA2.
        barcode (TrimCoordinate): The trim coordinate for the barcode.
    """
    gRNA1: TrimCoordinate
    gRNA2: TrimCoordinate
    barcode: TrimCoordinate


DEFAULT_TWO_READ_TRIM_STRATEGY = TrimStrategy(
    gRNA1=TrimCoordinate(file_index=0, start=0, end=20),
    gRNA2=TrimCoordinate(file_index=1, start=1, end=21),
    barcode=TrimCoordinate(file_index=1, start=160, end=166)
)

DEFAULT_THREE_READ_TRIM_STRATEGY = TrimStrategy(
    gRNA1=TrimCoordinate(file_index=0, start=0, end=20),
    gRNA2=TrimCoordinate(file_index=1, start=1, end=21),
    barcode=TrimCoordinate(file_index=2, start=0, end=6)
)
