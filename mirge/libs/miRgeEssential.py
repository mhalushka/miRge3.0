import subprocess
from pathlib import Path


def check_dependencies(args):
    bwtCommand = Path(args.bowtieBinary)/"bowtie --version" if args.bowtieBinary else "bowtie --version"
    samtoolsCommand = Path(args.samtoolsBinary)/"samtools --version" if args.samtoolsBinary else "samtools --version"
    cutadaptCommand = "cutadapt --version"
    rnaFoldCommand = Path(args.rnafoldBinary)/"RNAfold --version" if args.rnafoldBinary else "RNAfold --version"


    # Checking bowtie version #
    bowtie = subprocess.run(str(bwtCommand), shell=True, capture_output=True, text=True)
    bwtver = ["1.2.1", "1.2.2", "1.2.3"]
    if bowtie.returncode==0:
        if not (bowtie.stdout.split('\n')[0].split(' ')[2]) in bwtver:
            print("bowtie error!: incorrect version. Require - bowtie (1.2.1, 1.2.2 or 1.2.3) \nUse argument -pb <name of the directory>")
            exit()
        else:
            print("bowtie version: "+ str(bowtie.stdout.split('\n')[0].split(' ')[2]))
    else:
        print("bowtie error!: bowtie, command not found \nUse argument -pb <name of the directory>")
        exit()


    # Checking cutadapt version #
    cutadapt = subprocess.run(str(cutadaptCommand), shell=True, capture_output=True, text=True)
    if not cutadapt.returncode==0 and float(cutadapt.stdout.strip()) >= 2.7:
        print("cutadapt error!. Required: cutadapt =2.7")
        exit()
    else:
        print("cutadapt version: "+ str(cutadapt.stdout.strip()))


    # Checking samtools version #
    samtools = subprocess.run(str(samtoolsCommand), shell=True, capture_output=True, text=True)
    if not samtools.returncode==0 and float(samtools.stdout.split('\n')[0].split(' ')[1]) >= 1.5:
        print("Samtools error!: Can't locate or version incorrect. Require - samtools >1.5\nUse argument -ps <name of the directory>")
        exit()
    else:
        print("Samtools version: "+ str(samtools.stdout.split('\n')[0].split(' ')[1]))
    
    
    # Checking RNAfold version #
    if args.novel_miRNA: 
        rnafold = subprocess.run(str(rnaFoldCommand), shell=True, capture_output=True, text=True)
        if rnafold.returncode==0:
            if float(rnafold.stdout.split('\n')[0].split(' ')[1]) != "2.4.14":
                print("RNAfold error!: Can't locate or version incorrect. Require - RNAfold = 2.4.14\nUse argument -pr <name of the directory>")
                exit()
            else:
                print("RNAfold version: "+ str(rnafold.stdout.split('\n')[0].split(' ')[1]))
        else:
            print("RNAfold error!: RNAfold, command not found \nUse argument -pr <name of the directory>")
            exit()
