#!/usr/bin/env python
from xopen import xopen
import dnaio
import io
import os
import sys
import time
import concurrent.futures
import pandas as pd
from pathlib import Path

from cutadapt.adapters import warn_duplicate_adapters
from cutadapt.parser import AdapterParser
from cutadapt.modifiers import (LengthTagModifier, SuffixRemover, PrefixSuffixAdder,
        ZeroCapper, QualityTrimmer, UnconditionalCutter, NEndTrimmer, AdapterCutter,
        PairedAdapterCutterError, PairedAdapterCutter, NextseqQualityTrimmer, Shortener)



def parse_cutoffs(s):
    """
    FUNCTION ADOPTED FROM CUTADAPT 2.7, TO CREATE INPUT PARAMETERS ACCORDING TO CUTADAPT.
    Parse a string INT[,INT] into a two-element list of integers
    >>> parse_cutoffs("5") => [0, 5]
    >>> parse_cutoffs("6,7") => [6, 7]
    """
    try:
        cutoffs = [int(value) for value in s.split(",")]
    except ValueError as e:
        exit("Quality cutoff value not recognized: {}".format(e))
    
    if len(cutoffs) == 1:
        cutoffs = [0, cutoffs[0]]
    elif len(cutoffs) != 2:
        exit("Expected one value or two values separated by comma for the quality cutoff")
        
    return cutoffs



def add_unconditional_cutters(pipeline_add, cut1):
    """
    FUNCTION ADOPTED FROM CUTADAPT 2.7, TO CREATE INPUT PARAMETERS ACCORDING TO CUTADAPT.
    """
    for i, cut_arg in enumerate([cut1]):
        # cut_arg is a list
        if not cut_arg:
            continue
        
        if len(cut_arg) > 2:
            exit("You cannot remove bases from more than two ends.")
        if len(cut_arg) == 2 and cut_arg[0] * cut_arg[1] > 0:
            exit("You cannot remove bases from the same end twice.")
        for c in cut_arg:
            if c == 0:
               continue
            if i == 0:  # R1
                pipeline_add(UnconditionalCutter(c))


def stipulate(args):
    """
    REQUIRED TO CREATE ITERABLE FUNCTIONS TO RUN IN CUTADAPT 2.7. THIS FUNCTION IS CALLED ONLY ONE TIME. 
    """
    modifiers=[]
    pipeline_add = modifiers.append
    adapter_parser = AdapterParser(
            max_error_rate=args.error_rate,
            min_overlap=args.overlap,
            read_wildcards=args.match_read_wildcards,
            adapter_wildcards=args.match_adapter_wildcards,
            indels=args.indels,
         )
    adapters = adapter_parser.parse_multi(args.adapters)
    warn_duplicate_adapters(adapters)

    if args.nextseq_trim is not None:
        pipeline_add(NextseqQualityTrimmer(args.nextseq_trim, args.phred64))
    if args.quality_cutoff is not None:
        cutoffs = parse_cutoffs(args.quality_cutoff)
        pipeline_add(QualityTrimmer(cutoffs[0], cutoffs[1], args.phred64))

    adapter_cutter = None
    if adapters:
        adapter_cutter = AdapterCutter(adapters, args.times, args.action)
        pipeline_add(adapter_cutter)
    if args.trim_n:
        pipeline_add(NEndTrimmer())
    add_unconditional_cutters(pipeline_add, args.cut)
        
    return modifiers



def baking(args, inFileArray, inFileBaseArray, workDir):
    """
    THIS FUNCTION IS CALLED FIRST FROM THE miRge3.0. 
    THIS FUNCTION PREPARES FUNCTIONS REQUIRED TO RUN IN CUTADAPT 2.7 AND PARSE ONE FILE AT A TIME. 
    """
    global ingredients, threads, buffer_size, trimmed_reads, fasta, fileTowriteFasta, min_len
    global fileToparse,fp,fastaOut
    numlines=10000
    fasta = args.fasta
    threads = args.threads
    buffer_size = args.buffer_size
    min_len = args.minimum_length
    ingredients = stipulate(args)
    df_mirged=pd.DataFrame()
    complete_set=pd.DataFrame()
    begningTime = time.perf_counter()
    for index, FQfile in enumerate(inFileArray):
        start = time.perf_counter()
        finish2=finish3=finish4=finish5=0
        if fasta:
            fileTowriteFasta = Path(workDir)/(inFileBaseArray[index]+".fasta")
            fastaOut = open(fileTowriteFasta, "w+")
        future=[]
        fileToparse = Path(workDir)/(inFileBaseArray[index]+".txt")
        #fp = open(fileToparse, 'a+')

        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            with dnaio.open(FQfile, mode='r') as readers:
                readobj=[]
                count=trimmed=0
                for reads in readers:
                    count+=1
                    readobj.append(reads)
                    if len(readobj) == 1000000:
                        future = [executor.submit(cutadapt, readobj[i:i+numlines]) for i in range(0, len(readobj), numlines)]
                        readobj=[]
                future.extend([executor.submit(cutadapt, readobj[i:i+numlines]) for i in range(0, len(readobj), numlines)])
                readobj=[]

        finish2 = time.perf_counter()
        print(f'Cutadapt finished file {inFileBaseArray[index]} in {round(finish2-start, 4)} second(s)')
        
        if fasta:
            fastaOut.close()
        readDict={}
        with open(fileToparse, 'r') as fpin:
            for treads in fpin:
                if treads.strip() in readDict:
                    readDict[treads.strip()]+=1
                else:
                    readDict[treads.strip()]=1
        collapsed_df = pd.DataFrame(list(readDict.items()), columns=['Sequence', inFileBaseArray[index]])
        collapsed_df.set_index('Sequence',inplace = True)
        if len(inFileBaseArray) == 1:
            complete_set = collapsed_df 
            collapsed_df = pd.DataFrame() 
        elif len(inFileBaseArray) > 1:
            complete_set = complete_set.join(collapsed_df, how='outer')
            collapsed_df = pd.DataFrame() 
        complete_set = complete_set.fillna(0).astype(int)
        finish3 = time.perf_counter()
        print(f'Collapsing finished file {inFileBaseArray[index]} in {round(finish3-finish2, 4)} second(s)\n')
    #https://stackoverflow.com/questions/24039023/add-column-with-constant-value-to-pandas-dataframe
    initialFlags = ['annotFlag','exact miRNA','hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','isomiR miRNA']
    complete_set = complete_set.assign(**dict.fromkeys(initialFlags, '0'))

    finalColumns = initialFlags + inFileBaseArray 
    complete_set = complete_set.reindex(columns=finalColumns)

    #fileToCSV = Path(workDir)/"miRge3_collapsed.csv"
    #complete_set.to_csv(fileToCSV)
    finish4 = time.perf_counter()
    print(f'Matrix creationg finished in {round(finish4-finish3, 4)} second(s)\n')
    EndTime = time.perf_counter()
    print(f'Completed in {round(EndTime-begningTime, 4)} second(s)\n')
    for temp in inFileBaseArray:
        temp = temp+".txt"
        temp = Path(workDir)/temp
        os.remove(temp)
    return(complete_set)

        



# THIS IS WHERE EVERYTHIHNG HAPPENS - Modifiers, filters etc...
def cutadapt(fq):
    #trimmed_reads=[]
    fp = open(fileToparse, 'a+')
    for fqreads in fq:
        matches=[]
        for modifier in ingredients:
            fqreads = modifier(fqreads, matches)
        #print(fqreads.sequence)
        #trimmed_reads.append(fqreads.sequence)
        if int(len(fqreads.sequence)) >= int(min_len):
            fp.write(f"{fqreads.sequence}\n")

        if fasta:
            if min_len:
                if int(len(fqreads.sequence)) >= int(min_len):
                    fastaOut.write(f">{fqreads.name}\n{fqreads.sequence}\n")
            else:
                fastaOut.write(f">{fqreads.name}\n{fqreads.sequence}\n")
        #if int(len(fqreads.sequence)) >= int(min_len):
    #return trimmed_reads
    #return trimmed_seq_dict

