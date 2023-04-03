import enum


class FileMode(enum.Enum):
    READ = 'r'
    WRITE = 'w'
    APPEND = 'a'
    UPDATE = '+'
