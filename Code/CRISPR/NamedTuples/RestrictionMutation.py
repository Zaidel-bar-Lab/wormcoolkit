from typing import NamedTuple
from Code.CRISPR.NamedTuples.RestrictionSite import RestrictionSite


class RestrictionMutation(NamedTuple):
    restriction_site: RestrictionSite
    number_of_mutations: int
    mutated_sites: list
    mutated_strand: str
    codon_mutations: list
    reattachment_mutations: list
    repair_template: list
