from Code.CRISPR.NamedTuples.SequenceSites import SequenceSites
from Code.CRISPR.Enum.DNASection import DNASection
from typing import NamedTuple


class MutateSection(NamedTuple):
    number_of_mutations: int
    section_type: DNASection
    section_sites: SequenceSites



