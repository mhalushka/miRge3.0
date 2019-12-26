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
#from libs.trim_file import *
#from libs.trim_adapt import *
#from libs.feeding import *
from libs.digest import *
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
    if Path(file_list[0]).is_dir():
        file_list = [str(x) for x in Path(file_list[0]).iterdir() if x.is_file()]
        fastq_fullPath,base_names = validate_files(file_list)
    elif Path(file_list[0]).exists() and Path(file_list[0]).suffix in file_exts: # READ TXT OR CSV FILE HERE
        with open(file_list[0]) as file:
            lines = [line.strip() for line in file]
            fastq_fullPath, base_names = validate_files(lines)
    else:  # READ FASTQ OR FASTQ.gz FILES HERE
        fastq_fullPath, base_names = validate_files(file_list)
    print(f"\nmiRge will process {len(fastq_fullPath)} out of {len(file_list)} input files.")
    print(fastq_fullPath)
    print(base_names)
    print(args)
    #process_reads(args, fastq_fullPath)
    baking(args, fastq_fullPath, base_names)
    adapter = args.adapters
    numCPU = args.threads 
    cleanedReads = base_names[0]+'.trim.fastq'
    #phred, processed_count, kept_count = trim_file(fastq_fullPath[0], adapter, cleanedReads, int(numCPU))



if __name__ == '__main__':
    main()
    
