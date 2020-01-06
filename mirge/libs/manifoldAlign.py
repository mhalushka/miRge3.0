import subprocess
from pathlib import Path
import pandas as pd
import time
import os

from libs.miRgeEssential import UID

def createFastaInput(SequenceToAlign, bwtInput):
    with open(bwtInput, 'w') as wseq:
        for sequences in SequenceToAlign:
            wseq.write(">"+str(sequences)+"\n")
            wseq.write(str(sequences)+"\n")
    return True


def alignPlusParse(bwtExec, iter_number, pdDataFrame):
    colnames = list(pdDataFrame.columns)
    colToAct = 2 + int(iter_number)
    bowtie = subprocess.run(str(bwtExec), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
    if bowtie.returncode==0:
        bwtOut = bowtie.stdout
        bwtErr = bowtie.stderr
    lines = bwtOut.split('\n')
    sam_rows = [line.split('\t') for line in lines if not line.startswith('@') ]
    for nl in sam_rows:
        if nl != ['']:
            if nl[2] != "*":
                #print(str(iter_number) + "\t" + str(nl))
                pdDataFrame.at[nl[0], colnames[colToAct]] = nl[2]
                pdDataFrame.at[nl[0], colnames[1]] = 1
    return pdDataFrame

def bwtAlign(args,pdDataFrame,workDir,ref_db):
    """
    THIS FUNCTION COLLECTS DATAFRAME AND USER ARGUMENTS TO MAP TO VARIOUS DATABASES USING BOWTIE.
    """
    begningTime = time.perf_counter()
    bwtCommand = Path(args.bowtie_path)/"bowtie " if args.bowtie_path else "bowtie "
    bwtInput = Path(workDir)/"bwtInput.fasta"
    #outSam = Path(workDir)/"SeqToAnnot.sam"
    print("Alignment in progress ...")
    indexNames = ['_mirna_', '_hairpin_', '_mrna', '_mature_trna', '_pre_trna','_snorna','_rrna','_ncrna_others','_genome','_spike-in']
    parameters = [' -n 0 -f --norc -S --threads ', ' -n 1 -f --norc -S --threads ', ' -v 1 -f -a --best --strata --norc -S --threads ', ' -n 1 -f --norc -S --threads ', ' -n 1 -f --norc -S --threads ', 
            ' -n 1 -f --norc -S --threads ', ' -n 0 -f --norc -S --threads ', ' -n 0 -f --norc -S --threads ', ' -5 1 -3 2 -v 2 -f --norc --best -S --threads ', ' -n 0 -f --norc -S --threads ']
    if args.spikeIn:
        iterations = 9
    else:
        iterations = 8
    for bwt_iter in range(iterations):
        if bwt_iter == 0:
            SequenceToAlign = pdDataFrame[pdDataFrame['SeqLength'] <= 25].index.tolist()
            createFastaInput(SequenceToAlign, bwtInput)
            indexName  = str(args.organism_name) + str(indexNames[bwt_iter]) + str(ref_db)
            #indexName  = str(args.organism_name) + "_mirna_"+ str(ref_db)
            indexFiles = Path(args.libraries_path)/args.organism_name/"index.Libs"/indexName
            bwtExec = str(bwtCommand) + " " + str(indexFiles) + str(parameters[bwt_iter]) + str(args.threads) + " " + str(bwtInput) 
            #bwtExec = str(bwtCommand) + " " + str(indexFiles) + " -n 0 -f --norc -S --threads " + str(args.threads) + " " + str(bwtInput) 
            alignPlusParse(bwtExec, bwt_iter, pdDataFrame)
        
        elif bwt_iter == 1:
            SequenceToAlign = pdDataFrame[pdDataFrame['SeqLength'] > 25].index.tolist()
            createFastaInput(SequenceToAlign, bwtInput)
            indexName  = str(args.organism_name) + str(indexNames[bwt_iter]) + str(ref_db)
            #indexName  = str(args.organism_name) + "_hairpin_"+ str(ref_db)
            indexFiles = Path(args.libraries_path)/args.organism_name/"index.Libs"/indexName
            bwtExec = str(bwtCommand) + " " + str(indexFiles) + str(parameters[bwt_iter]) + str(args.threads) + " " + str(bwtInput) 
            #bwtExec = str(bwtCommand) + " " + str(indexFiles) + " -n 1 -f --norc -S --threads " + str(args.threads) + " " + str(bwtInput) 
            alignPlusParse(bwtExec, bwt_iter, pdDataFrame)

        else:
            SequenceToAlign = pdDataFrame[(pdDataFrame['annotFlag'] == '0')].index.tolist()
            createFastaInput(SequenceToAlign, bwtInput)
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
    pdMapped = pdDataFrame[pdDataFrame.annotFlag == 1]
    pdMapped = pdMapped.drop(columns=['SeqLength'])
    pdUnmapped = pdDataFrame[pdDataFrame.annotFlag == '0']
    pdUnmapped = pdUnmapped.drop(columns=['SeqLength'])
    fileToCSV = Path(workDir)/"miRge3_collapsed.csv"
    mappedfileToCSV = Path(workDir)/"mapped.csv"
    unmappedfileToCSV = Path(workDir)/"unmapped.csv"
    pdDataFrame.to_csv(fileToCSV)
    pdMapped.to_csv(mappedfileToCSV)
    pdUnmapped.to_csv(unmappedfileToCSV)
        
