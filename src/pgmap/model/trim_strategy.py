
import typing

from pgmap.model.trim_coordinate import TrimCoordinate


class TrimStrategy(typing.NamedTuple):
    R1: TrimCoordinate
    R2: TrimCoordinate
    barcode: TrimCoordinate
