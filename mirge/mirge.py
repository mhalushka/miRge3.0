#!/usr/bin/env python

#Built-in libraries 
from pathlib import Path
import time
import sys
import os

# GitHub libraries
import pandas
from cutadapt.modifiers import AdapterCutter, QualityTrimmer, UnconditionalCutter, QualityTrimmer
import cutadapt

#Custom miRge libraries 
from libs.parse import parseArg
from libs.miRgeEssential import check_dependencies, validate_files
#from libs.feeding import *
from libs.digest import baking 


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
    print(f"\nmiRge3.0 will process {len(fastq_fullPath)} out of {len(file_list)} input file(s).\n")
    pdDataFrame = baking(args, fastq_fullPath, base_names, workDir)
    fileToCSV = Path(workDir)/"miRge3_collapsed.csv"
    pdDataFrame.to_csv(fileToCSV)
    


if __name__ == '__main__':
    main()
    
