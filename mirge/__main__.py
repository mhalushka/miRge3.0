#!/usr/bin/env python

#Built-in libraries 
from pathlib import Path
import time
import sys
import os
import multiprocessing

# GitHub libraries
import pandas
from cutadapt.modifiers import AdapterCutter, QualityTrimmer, UnconditionalCutter, QualityTrimmer
import cutadapt

#Custom miRge libraries 
from mirge.libs.parse import parseArg
from mirge.libs.miRgeEssential import check_dependencies, validate_files
from mirge.libs.digest import baking 
from mirge.libs.summary import summarize
from mirge.libs.manifoldAlign import bwtAlign
from mirge.libs.novel_mir import predict_nmir


def main():
    globalstart = time.perf_counter()     
    args = parseArg()
    check_dependencies(args)
    samples = args.samples
    tStamp = time.strftime('%Y-%m-%d_%H-%M-%S',time.localtime(time.time()))
    ourDir = "miRge." + tStamp
    workDir = Path(args.outDir)/ourDir if args.outDir else Path.cwd()/ourDir
    Path(workDir).mkdir(exist_ok=True, parents=True)
    db_keys = {"mirbase":"miRBase", "mirgenedb":"MirGeneDB"}
    if args.tRNA_frag and args.organism_name != "human":
        sys.exit("ERROR: Detection of tRF(tRNA fragments) is only supported for human.")

    if args.threads == 0:
        args.threads = multiprocessing.cpu_count()

    ref_db = db_keys.get(args.mir_DB.lower()) if args.mir_DB.lower() in db_keys else sys.exit("ERROR: Require valid database (-d miRBase or MirGeneDB)")
    if len(args.adapters) == 2:
        back = list(args.adapters[0])
        if back[1] == "illumina":
            back[1] = 'TGGAATTCTCGGGTGCCAAGGAACTCCAG'
        args.adapters[0] = tuple(back)

        front = list(args.adapters[1])
        if front[1] == "illumina":
            front[1] = 'GTTCAGAGTTCTACAGTCCGACGATC'
        args.adapters[1] = tuple(front)

    if len(args.adapters) == 1:
        somewhere = list(args.adapters[0])
        if somewhere[0] == "back" and somewhere[1] == "illumina":
            somewhere[1] = 'TGGAATTCTCGGGTGCCAAGGAACTCCAG'
            args.adapters[0] = tuple(somewhere)
        elif somewhere[0] == "front" and somewhere[1] == "illumina": 
            somewhere[1] = 'GTTCAGAGTTCTACAGTCCGACGATC'
            args.adapters[0] = tuple(somewhere)

    print("Collecting and validating input files...")
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
    pdDataFrame,sampleReadCounts,trimmedReadCounts,trimmedReadCountsUnique = baking(args, fastq_fullPath, base_names, workDir)
    pdDataFrame = bwtAlign(args,pdDataFrame,workDir,ref_db)
    print(f"Summarizing and tabulating results...")
    summary_Start_time = time.perf_counter()
    pdMapped = pdDataFrame[pdDataFrame.annotFlag.eq(1)]
    pdUnmapped = pdDataFrame[pdDataFrame.annotFlag.eq(0)]
    summarize(args, workDir, ref_db, base_names, pdMapped, sampleReadCounts, trimmedReadCounts, trimmedReadCountsUnique)

    #fileToCSV = Path(workDir)/"miRge3_collapsed.csv"
    mappedfileToCSV = Path(workDir)/"mapped.csv"
    unmappedfileToCSV = Path(workDir)/"unmapped.csv"
    #pdDataFrame.to_csv(fileToCSV)
    pdMapped.to_csv(mappedfileToCSV)
    pdUnmapped.to_csv(unmappedfileToCSV)
    summary_End_time = time.perf_counter()
    print(f'Summary completed in {round(summary_End_time-summary_Start_time, 4)} second(s)\n')     
    if args.novel_miRNA:
        print("Predicting novel miRNAs\n")
        predict_nmir(args, workDir, ref_db, base_names, pdUnmapped)
    globalend_time = time.perf_counter()
    print(f'\nThe analysis completed in {round(globalend_time-globalstart, 4)} second(s)\n')     



if __name__ == '__main__':
    main()
    
