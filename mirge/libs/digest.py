#!/usr/bin/env python
from xopen import xopen
import dnaio
import io
import sys
import time
import concurrent.futures

from cutadapt.adapters import warn_duplicate_adapters
from cutadapt.parser import AdapterParser
from cutadapt.modifiers import (LengthTagModifier, SuffixRemover, PrefixSuffixAdder,
        ZeroCapper, QualityTrimmer, UnconditionalCutter, NEndTrimmer, AdapterCutter,
        PairedAdapterCutterError, PairedAdapterCutter, NextseqQualityTrimmer, Shortener)


start = time.perf_counter()

def parse_cutoffs(s):
    """
    FUNCTION ADOPTED FROM CUTADAPT 2.7, TO CREATE INPUT PARAMETERS ACCORDING TO CUTADAPT.
    Parse a string INT[,INT] into a two-element list of integers
    >>> parse_cutoffs("5")
    [0, 5]
    >>> parse_cutoffs("6,7")
    [6, 7]
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



def baking(args, inFileArray, inFileBaseArray):
    """
    THIS FUNCTION IS CALLED FIRST FROM THE miRge3.0. 
    THIS FUNCTION PREPARES FUNCTIONS REQUIRED TO RUN IN CUTADAPT 2.7 AND PARSE ONE FILE AT A TIME. 
    """
    global ingredients
    global threads
    global buffer_size
    numlines=10000
    threads = args.threads
    buffer_size = args.buffer_size
    ingredients = stipulate(args)
    for FQfile in inFileArray:
        readobj=[]
        with dnaio.open(FQfile, mode='r') as readers:
            for fqs in readers:
                readobj.append(fqs)

        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:

            future = [executor.submit(cutadapt, readobj[i:i+numlines]) for i in range(0, len(readobj), numlines)]
            #future = [executor.submit(distribute, readobj[i:i+numlines]) for i in range(0, len(readobj), numlines)]
#        for re in readobj:
#            cutadapt(re)

#        print(len(file_chunk_array))
#        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
#            eachFileChunkObject = [executor.submit(readfile, i) for i in file_chunk_array]
#            for outChunks in concurrent.futures.as_completed(eachFileChunkObject):
#                for fq in outChunks.result():
#                    print(fq)
        
        
        
    finish = time.perf_counter()
    print(f'Finished in {round(finish-start, 4)} second(s)')



# THIS IS WHERE EVERYTHIHNG HAPPENS - Modifiers, filters etc...
def cutadapt(fq):
    #print(fq)
    for fqreads in fq:
        matches=[]
        for modifier in ingredients:
            fqreads = modifier(fqreads, matches)
        print(fqreads.sequence)
    #return(fq.sequence)
