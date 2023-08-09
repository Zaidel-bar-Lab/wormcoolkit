from typing import NamedTuple
from Code.CRISPR.NamedTuples.SequenceSites import SequenceSites


class CodonData(NamedTuple):
    codon: str
    codon_sites: SequenceSites

# CodonData = namedtuple("CodonData", ['codon', 'codon_sites'])

