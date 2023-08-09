from typing import NamedTuple

from Code.CRISPR.Enum.RestrictionSiteType import RestrictionSiteType
from Code.CRISPR.NamedTuples.RestrictionEnzyme import RestrictionEnzyme
from Code.CRISPR.NamedTuples.SequenceSites import SequenceSites


class RestrictionSite(NamedTuple):
    index: SequenceSites
    enzyme: RestrictionEnzyme
    rest_site_type: RestrictionSiteType
