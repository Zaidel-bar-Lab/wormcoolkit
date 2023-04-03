import enum


class BlastType(enum.Enum):
    BLASTN = 'blastn'
    BLASTP = 'blastp'
    BLASTX = 'blastx'
    TBLASTN = 'tblastn'
    TBLASTX = 'tblasrx'

    def fromNtoP(b):
        if b is BlastType.BLASTN.value:
            return BlastType.BLASTP.value
        elif b is BlastType.BLASTP.value:
            return BlastType.BLASTN.value
