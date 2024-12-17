import typing


class TrimCoordinate(typing.NamedTuple):
    """
    A container describing a trim location within a read.

    Attributes:
        file_index (int): The zero indexed file index to trim from.
        start (int): The zero indexed, inclusive start location for the trim.
        end (int): The zero indexed, exclusive end location for the trim.
    """
    file_index: int
    start: int
    end: int
