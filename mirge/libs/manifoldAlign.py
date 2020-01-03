import subprocess
from pathlib import Path


def bwtAlign(args,pdDataFrame,workDir):
    """
    THIS FUNCTION LOOKS FOR DEPENDENCIES REQUIRED TO EXECUTE miRge3.0.
    """
    bwtCommand = Path(args.bowtie_path)/"bowtie " if args.bowtie_path else "bowtie "
    SequenceToAlign = pdDataFrame[pdDataFrame['SeqLength'] == 0].index.tolist()
    bwtInput = Path(workDir)/("bwtInput.fasta")
    with open(bwtInput, 'w') as wseq:
        for sequences in SequenceToAlign:
            wseq.write(">"+str(sequences))
            wseq.write(str(sequences))

        
