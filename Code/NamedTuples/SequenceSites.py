from typing import NamedTuple


class SequenceSites(NamedTuple):
    start: int
    end: int

    def get_mean(self):
        return (self.start + self.end)//2


