
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
