import enum


class RestrictionSiteType(enum.Enum):
    INSERTED_WITHOUT = '+'
    INSERTED = '++'
    REMOVED_WITHOUT = '-'
    REMOVED = '--'

