import subprocess
from pathlib import Path
import pandas as pd
import time
import os
import re
import concurrent.futures

from mirge.libs.miRgeEssential import UID


def alignPlusParse(bwtExec, iter_number, pdDataFrame, args, workDir):
    """
    ALIGN TO BOWTIE, PARSE SAM FILE AND UPDATE THE DATAFRAME
    """
    #indexNames = ['_mirna_', '_hairpin_', '_mature_trna', '_pre_trna', '_snorna', '_rrna', '_ncrna_others', '_mrna', '_mirna_', '_spike-in']
    colnames = list(pdDataFrame.columns)
    colToAct = 1 + int(iter_number)
    bowtie = subprocess.run(str(bwtExec), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
    if args.bam_out: 
        if iter_number == 0 or iter_number == 8:
            bwtoutput = Path(workDir)/"miRge3_miRNA.sam"
            bwto = open(bwtoutput, "a+")
        elif iter_number == 1:
            bwtoutput = Path(workDir)/"miRge3_hairpin_miRNA.sam"
            bwto = open(bwtoutput, "a+")
        elif iter_number == 4: 
            bwtoutput = Path(workDir)/"miRge3_snorna.sam"
            bwto = open(bwtoutput, "a+")
        elif iter_number == 5:
            bwtoutput = Path(workDir)/"miRge3_rrna.sam" 
            bwto = open(bwtoutput, "a+")
        elif iter_number == 6:
            bwtoutput = Path(workDir)/"miRge3_ncrna_others.sam" 
            bwto = open(bwtoutput, "a+")
        elif iter_number == 7:
            bwtoutput = Path(workDir)/"miRge3_mrna.sam" 
            bwto = open(bwtoutput, "a+")
    if args.tRNA_frag:
        if iter_number == 2:
            bwtoutput = Path(workDir)/"miRge3_tRNA.sam"
            bwto = open(bwtoutput, "a+")
        elif iter_number == 3:
            bwtoutput = Path(workDir)/"miRge3_pre_tRNA.sam"
            bwto = open(bwtoutput, "a+")

    if bowtie.returncode==0:
        bwtOut = bowtie.stdout
        bwtErr = bowtie.stderr
    for srow in bwtOut.split('\n'):
        if not srow.startswith('@'):
            sam_line = srow.split('\t')
            if sam_line != ['']:
                if sam_line[2] != "*":
                    pdDataFrame.at[sam_line[0], colnames[colToAct]] = sam_line[2]
                    pdDataFrame.at[sam_line[0], colnames[0]] = 1
                    if iter_number != 2 and iter_number != 3 and iter_number != 9: #tRNA and pre_tRNA and spike-ins are ignored 
                        if args.bam_out:
                            bwto.write(srow+"\n")
                    elif args.tRNA_frag:
                        if iter_number != 9:
                            bwto.write(srow+"\n")

    return pdDataFrame



def bwtAlign(args,pdDataFrame,workDir,ref_db):
    """
    THIS FUNCTION COLLECTS DATAFRAME AND USER ARGUMENTS TO MAP TO VARIOUS DATABASES USING BOWTIE. CALLED FIRST AND ONCE. 
    """
    global threads
    threads = args.threads
    begningTime = time.perf_counter()
    bwtCommand = str(Path(args.bowtie_path)/"bowtie ") if args.bowtie_path else "bowtie "
    if args.bowtieVersion == "False": # That is if version is v1.3.0 
        bwtCommand += " -x "
    bwtInput = Path(workDir)/"bwtInput.fasta"
    runlogFile = Path(workDir)/"run.log"
    outlog = open(str(runlogFile),"a+")
    if not args.quiet:
        print("Alignment in progress ...")
    outlog.write("Alignment in progress ...\n")
    indexNames = ['_mirna_', '_hairpin_', '_mature_trna', '_pre_trna', '_snorna', '_rrna', '_ncrna_others', '_mrna', '_mirna_', '_spike-in']
    parameters = [' -n 0 -f --norc -S --threads ', ' -n 1 -f --norc -S --threads ', ' -v 1 -f -a --best --strata --norc -S --threads ', ' -v 0 -f -a --best --strata --norc -S --threads ', ' -n 1 -f --norc -S --threads ', ' -n 1 -f --norc -S --threads ', ' -n 1 -f --norc -S --threads ', ' -n 0 -f --norc -S --threads ', ' -5 1 -3 2 -v 2 -f --norc --best -S --threads ', ' -n 0 -f --norc -S --threads ']
    if args.spikeIn:
        iterations = 10
    else:
        iterations = 9
    for bwt_iter in range(iterations):
        if bwt_iter == 0:
            with open(bwtInput, 'w') as wseq:
                for sequences in (pdDataFrame.index[pdDataFrame.index.str.len() < 26]):
                    wseq.write(">"+str(sequences)+"\n")
                    wseq.write(str(sequences)+"\n")

            indexName  = str(args.organism_name) + str(indexNames[bwt_iter]) + str(ref_db)
            indexFiles = Path(args.libraries_path)/args.organism_name/"index.Libs"/indexName
            bwtExec = str(bwtCommand) + str(indexFiles) + str(parameters[bwt_iter]) + str(args.threads) + " " + str(bwtInput) 
            alignPlusParse(bwtExec, bwt_iter, pdDataFrame, args, workDir)
        
        elif bwt_iter == 1:
            with open(bwtInput, 'w') as wseq: 
                for sequences in (pdDataFrame.index[pdDataFrame.index.str.len() > 25]):
                    wseq.write(">"+str(sequences)+"\n")
                    wseq.write(str(sequences)+"\n")

            indexName  = str(args.organism_name) + str(indexNames[bwt_iter]) + str(ref_db)
            indexFiles = Path(args.libraries_path)/args.organism_name/"index.Libs"/indexName
            bwtExec = str(bwtCommand) + str(indexFiles) + str(parameters[bwt_iter]) + str(args.threads) + " " + str(bwtInput) 
            alignPlusParse(bwtExec, bwt_iter, pdDataFrame, args, workDir)

        else:
            if bwt_iter == 8: 
                indexName  = str(args.organism_name) + str(indexNames[bwt_iter]) + str(ref_db)
            else:
                indexName  = str(args.organism_name) + str(indexNames[bwt_iter])
            if bwt_iter == 3:
                with open(bwtInput, 'w') as wseq:
                    for sequences in (pdDataFrame.index[pdDataFrame.annotFlag.eq(0)]):
                        try:
                            footer = sequences[:(re.search('T{3,}$', sequences).span(0)[0])]
                            wseq.write(">"+str(sequences)+"\n")
                            wseq.write(str(footer)+"\n")
                        except AttributeError:
                            pass
            else:
                with open(bwtInput, 'w') as wseq:
                    for sequences in (pdDataFrame.index[pdDataFrame.annotFlag.eq(0)]):
                        wseq.write(">"+str(sequences)+"\n")
                        wseq.write(str(sequences)+"\n")

            indexFiles = Path(args.libraries_path)/args.organism_name/"index.Libs"/indexName
            bwtExec = str(bwtCommand) + str(indexFiles) + str(parameters[bwt_iter]) + str(args.threads) + " " + str(bwtInput) 
            alignPlusParse(bwtExec, bwt_iter, pdDataFrame, args, workDir)
    finish = time.perf_counter()
    if not args.spikeIn:
        pdDataFrame = pdDataFrame.drop(columns=['spike-in'])
    
    os.remove(bwtInput)
    pdDataFrame = pdDataFrame.fillna('')
    if not args.quiet:
        print(f'Alignment completed in {round(finish-begningTime, 4)} second(s)\n')
    outlog.write(f'Alignment completed in {round(finish-begningTime, 4)} second(s)\n')
    outlog.close()
    return pdDataFrame
