#!/usr/bin/env python
import time
import math
import re
import pandas as pd
from pathlib import Path
import numpy as np
import subprocess
from Bio import SeqIO
#from difflib import unified_diff, Differ
#from libs.miRgeEssential import UID
import _pickle as cPickle
#import cPickle
import os, sys
import matplotlib
matplotlib.use('agg')
from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold, chi2, f_classif
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.metrics import make_scorer, confusion_matrix, matthews_corrcoef, roc_auc_score, roc_curve, auc
#from sklearn.model_selection import cross_validation
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.spatial import cKDTree
# from sklearn.externals import joblib 
# /home/arun/.local/lib/python3.8/site-packages/sklearn/externals/joblib/__init__.py:15: FutureWarning: sklearn.externals.joblib is deprecated in 0.21 and will be removed in 0.23. Please import this functionality directly from joblib, which can be installed with: pip install joblib. If this warning is raised when loading pickled models, you may need to re-serialize those models with scikit-learn 0.21+.
# warnings.warn(msg, category=FutureWarning)


def convert2Fasta(pdUnmapped, infile, minLength, maxLength, countCutoff, outputdir2, species, speciesNameDic, base_names):
    pdUnmapped = pdUnmapped.reset_index(level=['Sequence'])
    outf1 = open((Path(outputdir2)/"unmapped_mirna_raw.fa"), "w+")
    outf2 = open((Path(outputdir2)/"unmapped_mirna.fa"), "w+")
    cols1 = ["Sequence"] + base_names
    unmap = pd.DataFrame(pdUnmapped, columns= cols1)
    unmap['Sum'] =  unmap[base_names].sum(axis=1)
    finalColumns = ['Sequence','Sum']+base_names
    unmap = unmap.reindex(columns=finalColumns) # RAW dataframe
    unmap.index = unmap.index + 1 # Setting to start index from 1 instead of 0
    filtered_unmap = unmap[(unmap['Sequence'].str.len() >= minLength) & (unmap['Sequence'].str.len() <= maxLength) & (unmap['Sum'] >= countCutoff)] #Filtered Dataframe
    unmapList = unmap.reset_index().values.tolist()
    for index, row in enumerate(unmapList): 
        outf1.write(">mir"+str(row[0])+"_"+str(row[2])+"\n"+row[1]+"\n")
    filtered_unmapList = filtered_unmap.reset_index().values.tolist()
    for index, row in enumerate(filtered_unmapList):
        outf2.write(">mir"+str(row[0])+"_"+str(row[2])+"\n"+row[1]+"\n")
    for each_sample in base_names:
        current_fname_raw = "unmapped_mirna_"+each_sample+"_raw.fa"
        current_fname_filtered = "unmapped_mirna_"+each_sample+".fa"
        outfn_raw = open((Path(outputdir2)/current_fname_raw), "w+") 
        outfn_filtered = open((Path(outputdir2)/current_fname_filtered), "w+")
        Sample_DF_raw = pd.DataFrame(unmap[unmap[each_sample] >= 1], columns=['Sequence', each_sample]) 
        Sample_DF_filtered = pd.DataFrame(filtered_unmap[filtered_unmap[each_sample] >= countCutoff], columns=['Sequence', each_sample])
        Sample_DF_rawList = Sample_DF_raw.reset_index().values.tolist() 
        Sample_DF_filteredList = Sample_DF_filtered.reset_index().values.tolist()
        for index, row in enumerate(Sample_DF_rawList):
            outfn_raw.write(">mir"+str(row[0])+"_"+str(row[2])+"\n"+row[1]+"\n")
        outfn_raw.close()
        for index, row in enumerate(Sample_DF_filteredList):
            outfn_filtered.write(">mir"+str(row[0])+"_"+str(row[2])+"\n"+row[1]+"\n")
        outfn_filtered.close()
    step5e = time.perf_counter()
    outf1.close()
    outf2.close()
    


def predict_nmir(args, workDir, ref_db, base_names, pdUnmapped):
    predict_start_time = time.perf_counter()
    print('\nPerforming prediction of novel miRNAs...')
    print('Start to predict')
    samtoolsBinary = Path(args.samtools_path)/"samtools " if args.samtools_path else "samtools " 
    rnafoldBinary = Path(args.RNAfold_path)/"RNAfold " if args.RNAfold_path else "RNAfold "
    species = args.organism_name
    genomRepeats = species+"_genome_repeats.pckl"
    genomFasta = species+"_genome.pckl"
    genom_miRDB = species+"_mature_"+ref_db+".fa"
    genome_repeats = Path(args.libraries_path)/args.organism_name/"annotation.Libs"/genomRepeats
    genome_fa = Path(args.libraries_path)/args.organism_name/"fasta.Libs"/genomFasta
    mature_miRNA_fa = Path(args.libraries_path)/args.organism_name/"fasta.Libs"/genom_miRDB
    nameAbbrNameDic ={'human':'hsa', 'zebrafish':'dre', 'mouse':'mmu', 'rat':'rno', 'fruitfly':'dme', 'nematode':'cel'}
    abbrName = nameAbbrNameDic[species]
    speciesNameDic = {species:abbrName}
    if args.organism_name  not in ['human', 'mouse']:
        speciesType = 'others'
    else:
        speciesType = args.organism_name
    #wideSampleListFile = args.wideSampleListFile
    minLength = int(args.minLength)
    maxLength = int(args.maxLength)
    countCutoff = int(args.minReadCounts)
    mapping_loc = int(args.maxMappingLoci)
    seedLength = int(args.seedLength)
    # Here, the optimal overlapLenCutoff should be set to be 14, based on the trial of 8, 10 , 12
    overlapLenCutoff = int(args.overlapLenCutoff)
    clusterSeqLenCutoff = int(args.clusterLength)
    genome_idx_file = species+'_genome.1.ebwt'
    genome_idx = species+'_genome'
    indexPath = Path(args.libraries_path)/args.organism_name/'index.Libs'/genome_idx_file
    genome_index = Path(args.libraries_path)/args.organism_name/'index.Libs'/genome_idx
    if not Path(indexPath).exists():
        print(f'ERROR: The bowtie index file of {species}_genome.*.ebwt is not located at {indexPath}, please check it.')
        exit()
    try:
        with open(genome_repeats, 'rb') as f:
            repEleChrCoordinateDic = cPickle.load(f)
    except IOError:
        repEleChrCoordinateDic = {}
    with open(genome_fa, 'rb') as f:
        chrSeqDic = cPickle.load(f)
    chrSeqLenDic = {}
    for key in chrSeqDic.keys():
        chrSeqLenDic.update({key:len(chrSeqDic[key])})
    #Load genome miRNA fasta file
    exactmiRNASeqDic = {}
    for record in SeqIO.parse(mature_miRNA_fa, "fasta"):
        exactmiRNASeqDic.update({record.id:str(record.seq)})
    #Load the genome miRNA coordinate file
    mc_name = args.organism_name+'_'+ref_db+'.gff3'
    miRNA_coordinate = Path(args.libraries_path)/args.organism_name/"annotation.Libs"/mc_name
    miRNAchrCoordivateDic = {}
    strandLabel=""
    with open(miRNA_coordinate,"r") as inf1:
        for line1 in inf1:
            if line1[0] != "#":
                content = line1.strip().split("\t")
                if content[2] == "miRNA":
                    chr = content[0]
                    if ref_db == 'miRBase':
                        miRNAName = content[-1].split(";")[2].split("=")[-1]
                    else:
                        miRNAName = content[-1].split(";")[0].split("=")[-1]
                        startPos = int(content[3])
                        endPos = int(content[4])
                        strandLabel = content[6]
                    if chr not in miRNAchrCoordivateDic.keys(): 
                        #in miRNAchrCoordivateDic[chr], the fist list store the information about "+" forward strand match.
                        #and the second one store the information about "-" reverse strand match
                        miRNAchrCoordivateDic.update({chr:[[],[]]})
                    else:
                        pass
                    if strandLabel == "+":
                        miRNAchrCoordivateDic[chr][0].append((startPos, endPos, miRNAName))
                    elif strandLabel == "-":
                        miRNAchrCoordivateDic[chr][1].append((startPos, endPos, miRNAName))
                    else:
                        pass
                else:
                    pass
            else:
                pass
    ### Arguments preprocession is done!
    # Perform novel miRNA detection if selecting predict mode.

    # To modify the variable infile
    #infile = args.infile
    #Load the output file 'maaped.csv' and 'unmapped.csv' for downstream analysis representatively.
    outputdir2 = Path(workDir)/"unmapped_tmp"
    os.mkdir(outputdir2)

    bwtCmdTmp = Path(args.bowtie_path)/"bowtie " if args.bowtie_path else "bowtie "
    bwtBuildCmdTmp = Path(args.bowtie_path)/"bowtie-build " if args.bowtie_path else "bowtie-build " 
    samtoolsCmdTmp = Path(args.samtools_path)/"samtools " if args.samtools_path else "samtools "
    rnafoldCmdTmp = Path(args.RNAfold_path)/"RNAfold " if args.RNAfold_path else "RNAfold "
    infile = "unmapped.log"
    logFile = Path(outputdir2)/infile
    with open(logFile, 'w') as outfLog:
        time1 = time.perf_counter()
        outfLog.write(f'Reading and trimming unmapped reads from Pandas Library\n')
        convert2Fasta(pdUnmapped, infile, minLength, maxLength, countCutoff, outputdir2, species, speciesNameDic, base_names)
        #time.sleep(.555)
        time2 = time.perf_counter()
        outfLog.write('Execution time for reading and trimming unmapped Pandas Library took: %.4fs\n'%(time2-time1))
        print('Execution time for reading and trimming unmapped Pandas Library took: %.4fs\n'%(time2-time1))
        outfLog.write('********************\n')
        outfLog.flush()
