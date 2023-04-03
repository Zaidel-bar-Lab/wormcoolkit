import enum


class FileType(enum.Enum):
    TSV = 1
    CSV = 2
    HTML = 3
    XLS = 4
    UNCLEAR = 5
    FILE = 6
    CONSOLE = 7