from typing import NamedTuple


class PointMutation(NamedTuple):
    index: int
    old_nucleotide: str
    new_nucleotide: str

