import subprocess
from pathlib import Path


def check_dependencies(args):
    bwtCommand = Path(args.bowtie_path)/"bowtie --version" if args.bowtie_path else "bowtie --version"
    samtoolsCommand = Path(args.samtools_path)/"samtools --version" if args.samtools_path else "samtools --version"
    cutadaptCommand = "cutadapt --version"
    rnaFoldCommand = Path(args.RNAfold_path)/"RNAfold --version" if args.RNAfold_path else "RNAfold --version"
    

    # Checking bowtie version #
    bowtie = subprocess.run(str(bwtCommand), shell=True, capture_output=True, text=True)
    bwtver = ["1.2.1", "1.2.2", "1.2.3"]
    if bowtie.returncode==0:
        if not (bowtie.stdout.split('\n')[0].split(' ')[2]) in bwtver:
            print("bowtie error!: incorrect version. Require - bowtie (1.2.1, 1.2.2 or 1.2.3) \nUse argument -pbwt <name of the directory>")
            exit()
        else:
            print("bowtie version: "+ str(bowtie.stdout.split('\n')[0].split(' ')[2]))
    else:
        print("bowtie error!: bowtie, command not found \nUse argument -pbwt <name of the directory>")
        exit()

    # Checking cutadapt version #
    cutadapt = subprocess.run(str(cutadaptCommand), shell=True, capture_output=True, text=True)
    try:
        if not cutadapt.returncode==0 and float(cutadapt.stdout.strip()) >= 2.7:
            print("cutadapt error!. Required: cutadapt =2.7")
            exit()
        else:
            print("cutadapt version: "+ str(cutadapt.stdout.strip()))
    except ValueError:
        print("cutadapt error!: cutadapt not found\nPlease install cutadapt version = 2.7.")

    # Checking samtools version #
    samtools = subprocess.run(str(samtoolsCommand), shell=True, capture_output=True, text=True)
    if samtools.returncode==0:
        if not float(samtools.stdout.split('\n')[0].split(' ')[1]) >= 1.5:
    #if not samtools.returncode==0 and float(samtools.stdout.split('\n')[0].split(' ')[1]) >= 1.5:
            print("Samtools error!: incorrect version. Require - samtools >1.5\nUse argument -psam <name of the directory>")
            exit()
        else:
            print("Samtools version: "+ str(samtools.stdout.split('\n')[0].split(' ')[1]))
    else:
        print("Samtools error!: samtools, command not found\n Use argument -psam <name of the directory>")
        exit()
    
    
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

def validate_files(in_fileArray, fastq_fullPath):
    for files in in_fileArray:
       filetype = ''.join(Path(files).suffixes) if Path(files).suffix == ".gz" else Path(files).suffix
       if Path(files).exists() and (filetype == ".fastq" or filetype == ".fastq.gz"):
           fastq_fullPath.append(str(Path(files).resolve()))
       else:
           print(f"\nWARNING: File {files} does not exists!") if not Path(files).exists() else print(f"\nWARNING: File {files} is neither fastq or fastq.gz format!")
           print(f"Omitting file {files}")
    #Validating files in the list where the list should not be empty:
    if not fastq_fullPath:
        print("\nERROR!: No valid input files were available!\nPlease verify miRge -s arguments")
        exit()
    return fastq_fullPath
