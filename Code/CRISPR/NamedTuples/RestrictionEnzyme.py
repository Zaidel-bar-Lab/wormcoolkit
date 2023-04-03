from typing import NamedTuple


class RestrictionEnzyme(NamedTuple):
    name: str
    site: str
    derivatives: tuple
    full_site: str


