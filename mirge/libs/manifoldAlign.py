import subprocess
from pathlib import Path
import pandas as pd
import time
import os
import re
import concurrent.futures

from libs.miRgeEssential import UID



def createFastaInput(SequenceToAlign, bwtInput, bwt_iter):
    """
    CREATE FASTA FOR EACH ITERATIONS FOR BOWTIE ALIGNMENT 
    """
    with open(bwtInput, 'w') as wseq:
        for sequences in SequenceToAlign:
            if bwt_iter == 3:
                try:
                    footer = sequences[:(re.search('T{3,}$', sequences).span(0)[0])]
                    wseq.write(">"+str(sequences)+"\n")
                    wseq.write(str(footer)+"\n")
                except AttributeError:
                    pass
            else:
                wseq.write(">"+str(sequences)+"\n")
                wseq.write(str(sequences)+"\n")
    return True



def parseSAM(sam_rows):
    """
    FUNCTION FOR PARALLEL PROCESSING, WHICH OPTIMIZES MEMORY UTILIZATION BY PROCESSING ONLY FEW CHUNKS OF DATA AT ONCE
    """
    collective_list=[]
    for srow in sam_rows:
        if not srow.startswith('@'):
            sam_line = srow.split('\t')
            if sam_line != ['']:
                if sam_line[2] != "*":
                    variable_item = [sam_line[0], sam_line[2]]
                    collective_list.append(variable_item)
    return collective_list



def alignPlusParse(bwtExec, iter_number, pdDataFrame):
    """
    ALIGN TO BOWTIE, PARSE SAM FILE AND UPDATE THE DATAFRAME
    """
    colnames = list(pdDataFrame.columns)
    colToAct = 2 + int(iter_number)
    bowtie = subprocess.run(str(bwtExec), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
    if bowtie.returncode==0:
        bwtOut = bowtie.stdout
        bwtErr = bowtie.stderr
    readobj=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for item in bwtOut.split('\n'):
            readobj.append(item)
            if len(readobj) == 100000:
                future = [executor.submit(parseSAM, readobj[i:i+10000]) for i in range(0, len(readobj), 10000)]
                for sam_match in concurrent.futures.as_completed(future):
                    for each_list in sam_match.result():
                        pdDataFrame.at[each_list[0], colnames[colToAct]] = each_list[1]
                        pdDataFrame.at[each_list[0], colnames[1]] = 1
                readobj=[]
        future = [executor.submit(parseSAM, readobj[i:i+10000]) for i in range(0, len(readobj), 10000)]
        for sam_match in concurrent.futures.as_completed(future):
            for each_list in sam_match.result():
                pdDataFrame.at[each_list[0], colnames[colToAct]] = each_list[1]
                pdDataFrame.at[each_list[0], colnames[1]] = 1
        
        readobj=[]
    return pdDataFrame



def bwtAlign(args,pdDataFrame,workDir,ref_db):
    """
    THIS FUNCTION COLLECTS DATAFRAME AND USER ARGUMENTS TO MAP TO VARIOUS DATABASES USING BOWTIE. CALLED FIRST AND ONCE. 
    """
    global threads
    threads = args.threads
    begningTime = time.perf_counter()
    bwtCommand = Path(args.bowtie_path)/"bowtie " if args.bowtie_path else "bowtie "
    bwtInput = Path(workDir)/"bwtInput.fasta"
    #outSam = Path(workDir)/"SeqToAnnot.sam"
    print("Alignment in progress ...")
    indexNames = ['_mirna_', '_hairpin_', '_mature_trna', '_pre_trna', '_snorna', '_rrna', '_ncrna_others', '_mrna', '_mirna_', '_spike-in']
    parameters = [' -n 0 -f --norc -S --threads ', ' -n 1 -f --norc -S --threads ', ' -v 1 -f -a --best --strata --norc -S --threads ', ' -v 0 -f -a --best --strata --norc -S --threads ', ' -n 1 -f --norc -S --threads ', ' -n 1 -f --norc -S --threads ', ' -n 1 -f --norc -S --threads ', ' -n 0 -f --norc -S --threads ', ' -5 1 -3 2 -v 2 -f --norc --best -S --threads ', ' -n 0 -f --norc -S --threads ']
    if args.spikeIn:
        iterations = 10
    else:
        iterations = 9
    for bwt_iter in range(iterations):
        if bwt_iter == 0:
            SequenceToAlign = pdDataFrame[pdDataFrame['SeqLength'] <= 25].index.tolist()
            createFastaInput(SequenceToAlign, bwtInput, bwt_iter)
            indexName  = str(args.organism_name) + str(indexNames[bwt_iter]) + str(ref_db)
            #indexName  = str(args.organism_name) + "_mirna_"+ str(ref_db)
            indexFiles = Path(args.libraries_path)/args.organism_name/"index.Libs"/indexName
            bwtExec = str(bwtCommand) + " " + str(indexFiles) + str(parameters[bwt_iter]) + str(args.threads) + " " + str(bwtInput) 
            #bwtExec = str(bwtCommand) + " " + str(indexFiles) + " -n 0 -f --norc -S --threads " + str(args.threads) + " " + str(bwtInput) 
            alignPlusParse(bwtExec, bwt_iter, pdDataFrame)
        
        elif bwt_iter == 1:
            SequenceToAlign = pdDataFrame[pdDataFrame['SeqLength'] > 25].index.tolist()
            createFastaInput(SequenceToAlign, bwtInput, bwt_iter)
            indexName  = str(args.organism_name) + str(indexNames[bwt_iter]) + str(ref_db)
            #indexName  = str(args.organism_name) + "_hairpin_"+ str(ref_db)
            indexFiles = Path(args.libraries_path)/args.organism_name/"index.Libs"/indexName
            bwtExec = str(bwtCommand) + " " + str(indexFiles) + str(parameters[bwt_iter]) + str(args.threads) + " " + str(bwtInput) 
            #bwtExec = str(bwtCommand) + " " + str(indexFiles) + " -n 1 -f --norc -S --threads " + str(args.threads) + " " + str(bwtInput) 
            alignPlusParse(bwtExec, bwt_iter, pdDataFrame)

        else:
            SequenceToAlign = pdDataFrame[(pdDataFrame['annotFlag'] == '0')].index.tolist()
            createFastaInput(SequenceToAlign, bwtInput, bwt_iter)
            if bwt_iter == 8: 
                indexName  = str(args.organism_name) + str(indexNames[bwt_iter]) + str(ref_db)
            else:
                indexName  = str(args.organism_name) + str(indexNames[bwt_iter])
            indexFiles = Path(args.libraries_path)/args.organism_name/"index.Libs"/indexName
            bwtExec = str(bwtCommand) + " " + str(indexFiles) + str(parameters[bwt_iter]) + str(args.threads) + " " + str(bwtInput) 
            alignPlusParse(bwtExec, bwt_iter, pdDataFrame)
            #exit()
    finish = time.perf_counter()
    print(f'Alignment completed in {round(finish-begningTime, 4)} second(s)\n')
    if not args.spikeIn:
        pdDataFrame = pdDataFrame.drop(columns=['spike-in'])
    
    os.remove(bwtInput)
    pdDataFrame = pdDataFrame.fillna('')
    return pdDataFrame
