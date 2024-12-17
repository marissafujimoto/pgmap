
import typing

from pgmap.model.trim_coordinate import TrimCoordinate


class TrimStrategy(typing.NamedTuple):
    """
    A container describing a trimming stategy for paired guide sequencing.

    Attributes:
        R1 (TrimCoordinate): The trim coordinate for the gRNA1.
        R2 (TrimCoordinate): The trim coordinate for the gRNA2.
        barcode (TrimCoordinate): The trim coordinate for the barcode.
    """
    # TODO rename to gRNA1, gRNA2 instead of R1 R2?
    R1: TrimCoordinate
    R2: TrimCoordinate
    barcode: TrimCoordinate
