#!/usr/bin/env python

#Built-in libraries 
from pathlib import Path
import time
import sys
import os
import multiprocessing
import pickle

# GitHub libraries
import pandas as pd

#Custom miRge libraries 
from mirge.libs.parse import parseArg
from mirge.libs.miRgeEssential import check_dependencies, validate_files
from mirge.libs.digest import baking 
from mirge.libs.digestEC import bakingEC
from mirge.libs.summary import summarize
from mirge.libs.manifoldAlign import bwtAlign
from mirge.libs.novel_mir import predict_nmir
from mirge.classes.exportHTML import FormatHTML

def main():
    globalstart = time.perf_counter()     
    args = parseArg()
    samples = args.samples
    if args.outDirName:
        ourDir_n = str(args.outDirName)
        workDir = Path(args.outDir)/ourDir_n if args.outDir else Path.cwd()/ourDir_n
    else:
        tStamp = time.strftime('%Y-%m-%d_%H-%M-%S',time.localtime(time.time()))
        ourDir = "miRge." + tStamp
        workDir = Path(args.outDir)/ourDir if args.outDir else Path.cwd()/ourDir
    Path(workDir).mkdir(exist_ok=True, parents=True)
    db_keys = {"mirbase":"miRBase", "mirgenedb":"MirGeneDB"}
    runlogFile = Path(workDir)/"run.log"
    outlog = open(str(runlogFile),"a+")
    outlog.write(" ".join(sys.argv))
    outlog.write("\n")
    outlog.close()
    check_dependencies(args, str(runlogFile))
    outlog = open(str(runlogFile),"a+")
    if args.tRNA_frag and args.organism_name != "human":
        outlog.write("ERROR: Detection of tRF(tRNA fragments) is only supported for human.\n")
        sys.exit("ERROR: Detection of tRF(tRNA fragments) is only supported for human.")

    indexFiles = Path(args.libraries_path)/args.organism_name/"index.Libs"
    isExist = os.path.exists(indexFiles)
    if not isExist:
        print("\n ERROR: The path to miRge libraries is incorrect or does not exist!\n")
        outlog = open(str(runlogFile),"a+")
        outlog.write("\n ERROR: The path to miRge libraries is incorrect or does not exist!\n")
        outlog.close()
        exit()

    if args.threads == 0:
        args.threads = multiprocessing.cpu_count()

    ref_db = db_keys.get(args.mir_DB.lower()) if args.mir_DB.lower() in db_keys else sys.exit("ERROR: Require valid database (-d miRBase or MirGeneDB)")
    if args.organism_name == "hamster":
        if "mirbase" in args.mir_DB.lower():
            print("Library for hamster is not developed for miRBase, therefore, MirGeneDB is used\n")
        ref_db = "MirGeneDB"
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
    #print(args.adapters)
    if not args.quiet:
        print("Collecting and validating input files...")
    outlog.write("Collecting and validating input files...\n")
    file_exts = ['.txt', '.csv']
    file_list = samples[0].split(',')
    if Path(file_list[0]).is_dir():
        if args.resume:
            rootToPKL = str(Path(file_list[0]).absolute())
            pdf = str(Path(rootToPKL)/'collapsed.pkl')
            pdfacc = str(Path(rootToPKL)/'collapsed_accessories.pkl')
            if not Path(pdf).exists() and not Path(pdfacc).exists():
                print("\nERROR: The provided path doesn't contain pickle files with .pkl extensions! Please refer to documentation.\n")
                outlog.write("\nERROR: The provided path doesn't contain pickle files with .pkl extensions! Please refer to documentation.\n")
                outlog.close()
                exit()
            else:
                pdDataFrame = pd.read_pickle(pdf)
                with open(pdfacc, 'rb') as pklin:
                    bpkl = pickle.load(pklin)
                sampleReadCounts = bpkl[0]
                trimmedReadCounts = bpkl[1]
                trimmedReadCountsUnique = bpkl[2]
                fastq_fullPath = bpkl[3]
                base_names = bpkl[4]
        else:
            file_list = [str(x) for x in Path(file_list[0]).iterdir() if x.is_file()]
            file_list = sorted(file_list)
            fastq_fullPath,base_names = validate_files(args, file_list, str(runlogFile))
    elif Path(file_list[0]).exists() and Path(file_list[0]).suffix in file_exts: # READ TXT OR CSV FILE HERE
        with open(file_list[0]) as file:
            lines = [line.strip() for line in file]
            fastq_fullPath, base_names = validate_files(args, lines, str(runlogFile))
    else:  # READ FASTQ OR FASTQ.gz FILES HERE
        fastq_fullPath, base_names = validate_files(args, file_list, str(runlogFile))
    if not args.quiet and not args.resume:
        print(f"\nmiRge3.0 will process {len(fastq_fullPath)} out of {len(file_list)} input file(s).\n")
        outlog.write(f"\nmiRge3.0 will process {len(fastq_fullPath)} out of {len(file_list)} input file(s).\n\n")
        outlog.close()
    elif args.resume:
        print(f"\nmiRge3.0 will process {len(fastq_fullPath)} saved run(s) from binary pickle file.\n")
        outlog.write(f"\nmiRge3.0 will process {len(fastq_fullPath)} saved run(s) from binary pickle file.\n\n")
        outlog.close()


    if args.miREC and not args.resume:
        if args.uniq_mol_ids:
        #if args.umi:
            outlog = open(str(runlogFile),"a+")
            print(f"\nCurrently, error correction does not apply for the data with UMI sequences\n")
            outlog.write(f"\nCurrently, error correction does not apply for the data with UMI sequences\n")
            outlog.close()
            exit()
        else:
            pdDataFrame,sampleReadCounts,trimmedReadCounts,trimmedReadCountsUnique = bakingEC(args, fastq_fullPath, base_names, workDir)
    elif not args.resume:
        pdDataFrame,sampleReadCounts,trimmedReadCounts,trimmedReadCountsUnique = baking(args, fastq_fullPath, base_names, workDir)

    if args.save_pkl:
        pickleFile = str(Path(workDir)/'collapsed.pkl')
        pickleFile_acc = str(Path(workDir)/'collapsed_accessories.pkl')
        pdDataFrame.to_pickle(pickleFile)
        a = [sampleReadCounts, trimmedReadCounts,trimmedReadCountsUnique,fastq_fullPath, base_names]
        with open(pickleFile_acc, 'wb') as pklac:
            pickle.dump(a, pklac, protocol=pickle.HIGHEST_PROTOCOL)
    
    if args.save_pkl == args.resume and not args.save_pkl != "False":
        outlog = open(str(runlogFile),"a+")
        print(f"\nERROR: The arguments -pkl and -rr are mutually exclusive! Please refer to documentation.\n" )
        outlog.write(f"\nERROR: The arguments -pkl and -rr are mutually exclusive! Please refer to documentation.\n")
        outlog.close()
        exit()

    pdDataFrame = bwtAlign(args,pdDataFrame,workDir,ref_db)
    outlog = open(str(runlogFile),"a+")
    if not args.quiet:
        print(f"Summarizing and tabulating results...")
    outlog.write("\nSummarizing and tabulating results...\n")
    outlog.close()
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
    """
    Enabling Visualization HTML format
    """
    html = FormatHTML(workDir)
    html.beginHTML()
    html.histReadLen(len(base_names))
    if args.gff_out:
        html.isomirsTab(len(base_names), True)
    else:
        html.isomirsTab(len(base_names), False)
    
    html.exprTab(len(base_names))
    
    if args.uniq_mol_ids:
        html.umiTab(len(base_names), True)
    else:
        html.umiTab(len(base_names), False)

    outlog = open(str(runlogFile),"a+")
    if not args.quiet:
        print(f'Summary completed in {round(summary_End_time-summary_Start_time, 4)} second(s)\n')     
    outlog.write(f"Summary completed in {round(summary_End_time-summary_Start_time, 4)} second(s)\n")
    if args.novel_miRNA:
        html.novelTab(1)
        if not args.quiet:
            print("Predicting novel miRNAs\n")
        outlog.write("Predicting novel miRNAs\n")
        outlog.close()
        predict_nmir(args, workDir, ref_db, base_names, pdUnmapped)
        outlog = open(str(runlogFile),"a+")
    else:
        html.novelTab(0)
    #    novelTab
    for fname in os.listdir(str(Path(workDir))):
        if fname.endswith('.sam'):
            try:
                allSamFiles = Path(workDir)/"*.sam"
                os.system('rm -r %s'%(allSamFiles))
                break
            except OSError:    
                pass
    html.closeHTML()
    if args.diffex:
        if args.metadata:
            miRCountsFile = str(Path(workDir)/"miR.Counts.csv")
            diffOutFile = str(Path(workDir)/"miR.DiffExpn.txt")
            diffRData = str(Path(workDir)/"miR.DESeq2.RData")
            metaDataFile = str(Path(args.metadata).resolve())
            if not args.quiet:
                print("Performing differential expression...\n")
            outlog.write("Performing differential expression...\n")
            outlog.close()
            outlog = open(str(runlogFile),"a+")
            RscriptDirTmp = Path(__file__).resolve().parents[0]
            RscriptDir = Path(RscriptDirTmp)/('rScripts')/('difxp_DESeq.R')
            outvpdf = Path(workDir)/('mirge_diffex.volcano.pdf')
            outpcapdf = Path(workDir)/('mirge_diffex.pca.pdf')
            os.system('Rscript %s %s %s %s %s %s %s'%(RscriptDir, miRCountsFile, metaDataFile, outvpdf, outpcapdf, diffOutFile, diffRData))
        else: 
            print("Metadata file is required\n")
            # Rscript  /home/arun/repositories/Project_120919/mirge/rScripts/difxp_DESeq.R miR.Counts.csv DESmetadata.csv 

    globalend_time = time.perf_counter()
    if not args.quiet:
        resultsDir = str(Path(workDir).absolute())
        print(f'\nThe path to ourput directory: {resultsDir}')
        print(f'\nThe analysis completed in {round(globalend_time-globalstart, 4)} second(s)\n')     
    outlog.write(f"\nThe path to ourput directory: {resultsDir}")
    outlog.write(f"\nThe analysis completed in {round(globalend_time-globalstart, 4)} second(s)\n")
    outlog.close()



if __name__ == '__main__':
    main()
    
