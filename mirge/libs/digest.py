#!/usr/bin/env python
import dnaio
import os
import time
import concurrent.futures
import pandas as pd
from pathlib import Path
import numpy as np

from cutadapt.adapters import warn_duplicate_adapters
#from cutadapt.parser import AdapterParser
from cutadapt.parser import make_adapters_from_specifications
from cutadapt.modifiers import (LengthTagModifier, SuffixRemover, PrefixSuffixAdder,
        ZeroCapper, QualityTrimmer, UnconditionalCutter, NEndTrimmer, AdapterCutter,
        PairedAdapterCutterError, PairedAdapterCutter, NextseqQualityTrimmer, Shortener)
from mirge.classes.exportHTML import FormatJS


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
    if int(args.cutadaptVersion[0]) < 3:
        pass
        adapter_parser = AdapterParser(
            max_error_rate=args.error_rate,
            min_overlap=args.overlap,
            read_wildcards=args.match_read_wildcards,
            adapter_wildcards=args.match_adapter_wildcards,
            indels=args.indels,
         )
    else:
        search_parameters = dict(
                max_errors=args.error_rate,
                min_overlap=args.overlap,
                read_wildcards=args.match_read_wildcards,
                adapter_wildcards=args.match_adapter_wildcards,
                indels=args.indels,
                )
        #adapter_parser = AdapterParser(
        #    min_overlap=args.overlap,
        #    read_wildcards=args.match_read_wildcards,
        #    adapter_wildcards=args.match_adapter_wildcards,
        #    indels=args.indels,
        # )
    adapters = make_adapters_from_specifications(args.adapters, search_parameters)    
    #adapters = adapter_parser.parse_multi(args.adapters)
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
    #print(modifiers)    
    return modifiers



def baking(args, inFileArray, inFileBaseArray, workDir):
    """
    THIS FUNCTION IS CALLED FIRST FROM THE miRge3.0. 
    THIS FUNCTION PREPARES FUNCTIONS REQUIRED TO RUN IN CUTADAPT 2.7 AND PARSE ONE FILE AT A TIME. 
    """
    global ingredients, threads, buffer_size, trimmed_reads, fasta, fileTowriteFasta, min_len, umi, qiagenumi, qiaAdapter, ModificationInfo, cu_ver
    cu_ver = int(args.cutadaptVersion[0]) # Cutadapt version (cu_ver)
    if int(args.cutadaptVersion[0]) >= 3:
        from cutadapt.modifiers import ModificationInfo
    numlines=10000
    umi = args.uniq_mol_ids
    qiagenumi = args.qiagenumi
    fasta = args.fasta
    threads = args.threads
    buffer_size = args.buffer_size
    min_len = args.minimum_length
    if qiagenumi:
        qiaAdapter = str(args.adapters[0][1])
    ingredients = stipulate(args)
    df_mirged=pd.DataFrame()
    complete_set=pd.DataFrame()
    begningTime = time.perf_counter()
    sampleReadCounts={}
    trimmedReadCounts={}
    trimmedReadCountsUnique={}
    digestReadCounts={}
    visual_treat = {'rlen':{}, 'hist':{}}
    runlogFile = Path(workDir)/"run.log"
    outlog = open(str(runlogFile),"a+")
    for index, FQfile in enumerate(inFileArray):
        start = time.perf_counter()
        finish2=finish3=finish4=finish5=0
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            with dnaio.open(FQfile, mode='r') as readers:
                readobj=[]
                count=trimmed=0
                completeDict = {}
                for reads in readers:
                    count+=1
                    readobj.append(reads)
                    if len(readobj) == 1000000:
                        future = [executor.submit(cutadapt, readobj[i:i+numlines]) for i in range(0, len(readobj), numlines)] # sending bunch of reads (#1000000) for parallel execution
                        for fqres_pairs in concurrent.futures.as_completed(future): 
                            for each_list in fqres_pairs.result(): # retreving results from parallel execution
                                varx = list(each_list)
                                try:
                                    visual_treat['rlen'][str(inFileBaseArray[index])].append(len(varx[0]))
                                except KeyError:
                                    visual_treat['rlen'][str(inFileBaseArray[index])]= [len(varx[0])]
                                if varx[0] in completeDict: # Collapsing, i.e., counting the occurance of each read for each data 
                                    completeDict[varx[0]] += int(varx[1])
                                    trimmed+=int(varx[1])
                                else:
                                    completeDict[varx[0]] = int(varx[1])
                                    trimmed+=int(varx[1])
                        readobj=[]
                future=[]
                future.extend([executor.submit(cutadapt, readobj[i:i+numlines]) for i in range(0, len(readobj), numlines)]) # sending remaining reads for parallel execution
                for fqres_pairs in concurrent.futures.as_completed(future):
                    for each_list in fqres_pairs.result(): # retreving results from parallel execution
                        varx = list(each_list)
                        try:
                            visual_treat['rlen'][str(inFileBaseArray[index])].append(len(varx[0]))
                        except KeyError:
                            visual_treat['rlen'][str(inFileBaseArray[index])]= [len(varx[0])]
                        if varx[0] in completeDict: # Collapsing, i.e., counting the occurance of each read for each data 
                            completeDict[varx[0]] += int(varx[1])
                            trimmed+=int(varx[1])
                        else:
                            completeDict[varx[0]] = int(varx[1])
                            trimmed+=int(varx[1])
                readobj=[]
        if umi:
            trimmed=0
            umicompleteDict=dict()
            umi_cut = umi.split(",")
            if not args.umiDedup:
                for s, c in completeDict.items():
                    pureSeq, cutumiSeq = UMIParser(s, int(umi_cut[0]), int(umi_cut[1]))
                    try:
                        visual_treat['hist'][str(inFileBaseArray[index])].append(c)
                    except KeyError:
                        visual_treat['hist'][str(inFileBaseArray[index])] = [c]
                    if len(pureSeq) >= int(min_len):
                        if pureSeq in umicompleteDict:
                            umicompleteDict[pureSeq] += c
                            trimmed += c
                        else:
                            umicompleteDict[pureSeq] = c
                            trimmed += c
                completeDict = umicompleteDict 
                digestReadCounts = {inFileBaseArray[index]:trimmed}
            elif args.umiDedup:
                intermediate_umi = inFileBaseArray[index]+"_umiCounts.csv"
                temp_umiFile = Path(workDir)/intermediate_umi
                iumiFile = open(temp_umiFile,"a+")
                iumiFile.write("UMISeq" + "," +"transcriptSeq," + "UMICounts" +"\n")
                for s, c in completeDict.items():
                    pureSeq, cutumiSeq = UMIParser(s, int(umi_cut[0]), int(umi_cut[1]))
                    try:
                        visual_treat['hist'][str(inFileBaseArray[index])].append(c)
                    except KeyError:
                        visual_treat['hist'][str(inFileBaseArray[index])] = [c]
                    if len(pureSeq) >= int(min_len):
                        iumiFile.write(str(cutumiSeq) + "," + str(pureSeq) + "," + str(c) +"\n")
                        if pureSeq in umicompleteDict:
                            umicompleteDict[pureSeq] += 1
                            trimmed += 1
                        else:
                            umicompleteDict[pureSeq] = 1
                            trimmed += 1
                completeDict = umicompleteDict 
                digestReadCounts = {inFileBaseArray[index]:trimmed}
                iumiFile.close()
        else:
            visual_treat['hist'][str(inFileBaseArray[index])] = []
            digestReadCounts = {inFileBaseArray[index]:trimmed}
            # umi_seq = umi_seq[:umi_cut[1]]
            # max_ad = 19 + int(umi_cut[1])
            # umi_seq = umi_seq[:max_ad][-int(umi_cut[1]):]
        uniqTrimmedReads = {inFileBaseArray[index]:len(completeDict)}
        #digestReadCounts = {inFileBaseArray[index]:sum(completeDict.values())}
        inputReadCounts = {inFileBaseArray[index]:count}
        sampleReadCounts.update(inputReadCounts)
        trimmedReadCounts.update(digestReadCounts)
        trimmedReadCountsUnique.update(uniqTrimmedReads)
        finish2 = time.perf_counter()
        if not args.quiet:
            print(f'Cutadapt finished for file {inFileBaseArray[index]} in {round(finish2-start, 4)} second(s)')
        outlog.write(f'Cutadapt finished for file {inFileBaseArray[index]} in {round(finish2-start, 4)} second(s)\n')
        """
        CREATING PANDAS MATRIX FOR ALL THE SAMPLES THAT CAME THROUGH 
        WILL BE EDITED TO A FUNCTION, ONCE UMI COMES IN PICTURE
        """
        if args.tcf_out:
            fno = str(inFileBaseArray[index]) + '.trim.collapse.fa'
            fo_tcf = Path(workDir)/fno
            header_count = 1
            with open(fo_tcf,'w') as fo:
                for seqs in sorted(completeDict, key=completeDict.get, reverse=True):
                    head_r = ">seq"+str(header_count)+"_"+str(completeDict[seqs])+"\n"
                    fo.write(head_r)
                    fo.write(str(seqs)+"\n")
                    header_count+=1

        collapsed_df = pd.DataFrame(list(completeDict.items()), columns=['Sequence', inFileBaseArray[index]])
        collapsed_df.set_index('Sequence',inplace = True)
        if len(inFileBaseArray) == 1:
            complete_set = collapsed_df 
            collapsed_df = pd.DataFrame() 
        elif len(inFileBaseArray) > 1:
            complete_set = complete_set.join(collapsed_df, how='outer')
            collapsed_df = pd.DataFrame() 
        complete_set = complete_set.fillna(0).astype(int)

        finish3 = time.perf_counter()
        if not args.quiet:
            print(f'Collapsing finished for file {inFileBaseArray[index]} in {round(finish3-finish2, 4)} second(s)\n')
        outlog.write(f'Collapsing finished for file {inFileBaseArray[index]} in {round(finish3-finish2, 4)} second(s)\n')
    
    #complete_set['SeqLength'] = complete_set.index.str.len()
    initialFlags = ['exact miRNA','hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','isomiR miRNA','spike-in'] # keeping other columns ready for next assignment
    complete_set = complete_set.assign(**dict.fromkeys(initialFlags, ''))
    annotFlags = ['annotFlag']
    complete_set = complete_set.assign(**dict.fromkeys(annotFlags, '0'))
    #lengthCol = ['SeqLength']
    finalColumns = annotFlags +initialFlags + inFileBaseArray # rearranging the columns as we want  
    #finalColumns = lengthCol + annotFlags +initialFlags + inFileBaseArray # rearranging the columns as we want  
    complete_set = complete_set.reindex(columns=finalColumns)
    complete_set = complete_set.astype({"annotFlag": int})
    finish4 = time.perf_counter()
    if not args.quiet:
        print(f'Matrix creation finished in {round(finish4-finish3, 4)} second(s)\n')
    outlog.write(f'Matrix creation finished in {round(finish4-finish3, 4)} second(s)\n')
    #for x, y in visual_treat.items():
    #    print(str(x)+"\n") 
    #    print("Arun\n")
        #print(y)
    histData = FormatJS(workDir)
    div_idnum = 1
    for sample_files in inFileBaseArray:
        rlenDistID = "readLengthID_" + str(div_idnum)
        val = visual_treat['rlen'][sample_files]
        #print(sample_files)
        maxVal = sorted(val)[-1]
        minVal = sorted(val)[0]
        hist, bins = np.histogram(val, bins=(maxVal-minVal))
        hist = hist.tolist()
        #print(hist, bins)
        histData.readLenDist(rlenDistID, sample_files, str(hist), str(list(bins)))
        if umi:
            umiDistID = "umiDivID_" + str(div_idnum)
            val = visual_treat['hist'][sample_files]
            #val = visual_treat['rlen'][sample_files]
            #print(sample_files)
            maxVal = sorted(val)[-1]
            #maxVal = val.sort()[-1]
            minVal = sorted(val)[0]
            hist, bins = np.histogram(val, bins=(maxVal-minVal))
            hist = hist.tolist()
            bins = bins.tolist()
            #print(hist, bins)
            histData.sampleUMIDist(umiDistID, sample_files, str(hist), str(list(bins)))
        div_idnum += 1
    #print(visual_treat['hist'])
    EndTime = time.perf_counter()
    if not args.quiet:
        print(f'Data pre-processing completed in {round(EndTime-begningTime, 4)} second(s)\n')
    outlog.write(f'\nData pre-processing completed in {round(EndTime-begningTime, 4)} second(s)\n\n')
    outlog.close()
    return(complete_set, sampleReadCounts, trimmedReadCounts, trimmedReadCountsUnique)


def UMIParser(s, f, b):
    front = ""
    end = ""
    front = s[:f]
    if int(b) != 0:
        center = s[f:-b]
    else:
        center = s[f:]
    end = s[-b:]
    umiseqreq = front + end
    return (center, umiseqreq)
    #return (front, center, end)


# THIS IS WHERE EVERYTHIHNG HAPPENS - Modifiers, filters etc...
def cutadapt(fq):
    readDict={}
    for fqreads in fq:
        if int(cu_ver) >= 3:
            info = ModificationInfo(None)
            info.matches=[]
        else:
            matches=[]
        if qiagenumi:
            currentSeq = fqreads.sequence
            umi_seq = ""
            for modifier in ingredients:
                if int(cu_ver) >= 3:
                    fqreads = modifier(fqreads, info)
                else:
                    fqreads = modifier(fqreads, matches)
            try:
                umi_seq = currentSeq.split(str(fqreads.sequence))[1]
                umi_cut = umi.split(",")
                max_ad = len(qiaAdapter) + int(umi_cut[1])
                umi_seq = umi_seq[:max_ad][-int(umi_cut[1]):]
            except ValueError:
                umi_seq = ""
            final_seq = fqreads.sequence + umi_seq
            if int(len(final_seq)) >= int(min_len):
                if str(final_seq) in readDict:
                    readDict[str(final_seq)]+=1
                else:
                    readDict[str(final_seq)]=1
        else:
            for modifier in ingredients:
                if int(cu_ver) >= 3:
                    fqreads = modifier(fqreads, info)
                else:
                    fqreads = modifier(fqreads, matches)
            if int(len(fqreads.sequence)) >= int(min_len):
                #print(fqreads.sequence)
                if str(fqreads.sequence) in readDict:
                    readDict[str(fqreads.sequence)]+=1
                else:
                    readDict[str(fqreads.sequence)]=1
    trimmed_pairs = list(readDict.items())
    return trimmed_pairs
