#!/usr/bin/env python

#Built-in libraries 
from pathlib import Path
import time
import sys

# GitHub libraries
#import gzip
from cutadapt.modifiers import AdapterCutter, QualityTrimmer, UnconditionalCutter, QualityTrimmer
import cutadapt

#Custom miRge libraries 
from libs.parse import parseArg
from libs.miRgeEssential import check_dependencies, validate_files

#with gzip.open('/home/joe/file.txt.gz', 'rb') as f:

def main():
    #WORKING 
    args = parseArg()
    check_dependencies(args)
    
    #TESTING 
    samples = args.samples
    tStamp = time.strftime('%Y-%m-%d_%H-%M-%S',time.localtime(time.time()))
    ourDir = "miRge." + tStamp
    workDir = Path(args.outDir)/ourDir if args.outDir else Path.cwd()/ourDir
    Path(workDir).mkdir(exist_ok=True, parents=True)
    db_keys = {"mirbase":"miRBase", "mirgenedb":"MirGeneDB"}
    ref_db = db_keys.get(args.mir_DB.lower()) if args.mir_DB.lower() in db_keys else sys.exit("ERROR: Require valid database (-d miRBase or MirGeneDB)")
    

    file_exts = ['.txt', '.csv']
    file_list = samples[0].split(',')
    fastq_fullPath = [] #Initiating an empty list to collect all the fastq or fastq.gz files based on user arguments!
    if Path(file_list[0]).is_dir():
        file_list = [str(x) for x in Path(file_list[0]).iterdir() if x.is_file()]
        fastq_fullPath = validate_files(file_list, fastq_fullPath)
    elif Path(file_list[0]).exists() and Path(file_list[0]).suffix in file_exts: # READ TXT OR CSV FILE HERE
        with open(file_list[0]) as file:
            lines = [line.strip() for line in file]
            fastq_fullPath = validate_files(lines, fastq_fullPath)
    else:  # READ FASTQ OR FASTQ.gz FILES HERE
        fastq_fullPath = validate_files(file_list, fastq_fullPath)
    print(fastq_fullPath)
    base_names = [Path(bn).stem for bn in fastq_fullPath]
    print(base_names)

if __name__ == '__main__':
    main()
    
