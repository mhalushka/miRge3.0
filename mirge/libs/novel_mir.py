#!/usr/bin/env python
import time
import math
import re
import pandas as pd
from pathlib import Path
import numpy as np
import subprocess
from Bio import SeqIO
import _pickle as cPickle
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
import re
from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna
import numpy as np
from mirge.libs.processSam import split_fasta_from_sam, combineSam, decorateSam, parse_refine_sam
from mirge.libs.generate_featureFiles import generate_featureFiles, get_precursors, renameStrFile
from mirge.libs.screen_precusor_candidates import screen_precusor_candidates
from mirge.libs.preprocess_featureFiles import preprocess_featureFiles, model_predict
from mirge.libs.write_novel_report import write_novel_report
from mirge.classes.exportHTML import FormatJS
# from sklearn.externals import joblib 
# /home/arun/.local/lib/python3.8/site-packages/sklearn/externals/joblib/__init__.py:15: FutureWarning: sklearn.externals.joblib is deprecated in 0.21 and will be removed in 0.23. Please import this functionality directly from joblib, which can be installed with: pip install joblib. If this warning is raised when loading pickled models, you may need to re-serialize those models with scikit-learn 0.21+.
# warnings.warn(msg, category=FutureWarning)


def convert2Fasta(pdUnmapped, infile, minLength, maxLength, countCutoff, outputdir2, species, speciesNameDic, base_names, rawReadCounts, filteredReadCounts):
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
        #outfn_raw = open((Path(outputdir2)/current_fname_raw), "w+") 
        outfn_filtered = open((Path(outputdir2)/current_fname_filtered), "w+")
        Sample_DF_raw = pd.DataFrame(unmap[unmap[each_sample] >= 1], columns=['Sequence', each_sample]) 
        rawReadCounts[each_sample] = Sample_DF_raw.shape[0]
        Sample_DF_filtered = pd.DataFrame(filtered_unmap[filtered_unmap[each_sample] >= countCutoff], columns=['Sequence', each_sample])
        filteredReadCounts[each_sample] = Sample_DF_filtered.shape[0]
        Sample_DF_rawList = Sample_DF_raw.reset_index().values.tolist() 
        Sample_DF_filteredList = Sample_DF_filtered.reset_index().values.tolist()
        #for index, row in enumerate(Sample_DF_rawList):
        #    outfn_raw.write(">mir"+str(row[0])+"_"+str(row[2])+"\n"+row[1]+"\n")
        #outfn_raw.close()
        for index, row in enumerate(Sample_DF_filteredList):
            outfn_filtered.write(">mir"+str(row[0])+"_"+str(row[2])+"\n"+row[1]+"\n")
        outfn_filtered.close()
    outf1.close()
    outf2.close()
    
def cluster_basedon_location(outfile2, overlapLenThresholdTmp,files, outputdir2, cluster_file):
    # The arguments of this function are: unmapped_mirna_vs_genome_sorted.sam(Or mapped_mirna_vs_genome_sorted.sam) overlapLenThreshold(12, this value must be less than  16)
    overlapLenThreshold = int(overlapLenThresholdTmp)
    inf = open(outfile2, 'r')
    chrList = []
    chrSeqDic = {}
    for line in inf:
        if line[0] != '@':
            contentTemp = line.strip().split('\t')
            chr = contentTemp[2]
            if 'chr' in chr:
                if chr not in chrList:
                    chrList.append(chr)
                    chrSeqDic.update({chr:[[],[]]})
                    #in chrSeqDic[chr], the fist list store the information about "+" forward strand match. #and the second one store the information about "-" reverse strand match
                flag = contentTemp[1]
                startPos = int(contentTemp[3])
                seq = contentTemp[9]
                endPos = startPos+len(seq)-1
                mirID = contentTemp[0]
                if flag == "0" and len(chrSeqDic[chr][0]) == 0:
                    chrSeqDic[chr][0].append([startPos, endPos, seq, [mirID]])
                elif flag == "16" and len(chrSeqDic[chr][1]) == 0:
                    chrSeqDic[chr][1].append([startPos, endPos, seq, [mirID]])
                else:
                    if flag == "0":
                        temp = chrSeqDic[chr][0][-1]
                        if startPos >= temp[0] and startPos <= temp[1] and temp[1]-startPos+1 >= overlapLenThreshold:
                            mirIDList = temp[3]
                            mirIDList.append(mirID)
                            startPosNew = temp[0]
                            if endPos <= temp[1]:
                                endPosNew = temp[1]
                                seqNew = temp[2]
                            else:
                                endPosNew = endPos
                                seqNew = temp[2] + seq[temp[1]-startPos+1:]
                            chrSeqDic[chr][0][-1] = [startPosNew, endPosNew, seqNew, mirIDList]
                        else:
                            chrSeqDic[chr][0].append([startPos, endPos, seq,[mirID]])
                    elif flag == "16":
                        temp = chrSeqDic[chr][1][-1]
                        if startPos >= temp[0] and startPos <= temp[1]and temp[1]-startPos+1 >= overlapLenThreshold:
                            mirIDList = temp[3]
                            mirIDList.append(mirID)
                            startPosNew = temp[0]
                            if endPos <= temp[1]:
                                endPosNew = temp[1]
                                seqNew = temp[2]
                            else:
                                endPosNew = endPos
                                seqNew = temp[2] + seq[temp[1]-startPos+1:]
                            chrSeqDic[chr][1][-1] = [startPosNew, endPosNew, seqNew, mirIDList]
                    else:
                        chrSeqDic[chr][1].append([startPos, endPos, seq,[mirID]])
    inf.close()
    with open(cluster_file,'w') as outf:
        i = 1
        outf.write('miRClusterID\tChr\tStrand\tStart\tEnd\tSequence\tSequenceLenght\tCoutOfReads\tCountOfMembers\tMembers\n')
        for chr in chrList:
            for content in chrSeqDic[chr][0]:
                readCount = 0
                for item in content[3]:
                    readCount = readCount + int(item.split('_')[1])
                outf.write('\t'.join([files+':'+'miRCluster_'+str(i)+'_'+str(len(content[2])),chr, "+", str(content[0]), str(content[1]), content[2],str(len(content[2])), str(readCount), str(len(content[3])),','.join(content[3])])+ '\n')
                i = i + 1
            for content in chrSeqDic[chr][1]:
                readCount = 0
                for item in content[3]:
                    readCount = readCount + int(item.split('_')[1])
                outf.write('\t'.join([files+':'+'miRCluster_'+str(i)+'_'+str(len(content[2])),chr, "-", str(content[0]), str(content[1]), content[2],str(len(content[2])), str(readCount), str(len(content[3])),','.join(content[3])])+ '\n')
                i = i + 1
    return i

def preTrimClusteredSeq(CoordinateDic, cluster_file, files, clusterSeqLenCutoff, clusterTrimedFile, clusterTrimedFile_orig_FASTA):
    inf1 = open(cluster_file, 'r')
    outf = open(clusterTrimedFile, 'w')
    outFA = open(clusterTrimedFile_orig_FASTA, 'w')
    line = inf1.readline()
    outf.write("\t".join(['miRClusterID\tOriginalSeq\tFlag\trepetitiveElementName']+line.split('\t')[1:]))
    line = inf1.readline()    
    trimmed_count=1
    while line != "":
        # flag = "0" means the seq will be abandoned. flag = "1" means the seq will be retained.
        flag = "0"
        repetitiveElementName = '*'
        content = line.strip().split('\t')
        # remove representative seqs with ength larger than 30
        if content[2] == '-':
            seq = Seq(content[5])
            orginalSeq = str(seq.reverse_complement())
        else:
            orginalSeq = content[5]
        if int(content[6]) <= int(clusterSeqLenCutoff):
            chr = content[1]
            start = int(content[3])
            end = int(content[4])
            # remove representative seqs with "AAAAAA" in the 3' and "TTTTTT" in the 5'
            if (re.search('A{6,}$', orginalSeq) is None) and (re.search('^T{6,}', orginalSeq) is None):
                if chr in CoordinateDic.keys():
                    dist, idx =CoordinateDic[chr][0][0].query([(start, 0)], 2)
                    if len(idx) == 1 and np.isfinite(dist[0][0]) and np.isfinite(dist[0][1]):
                        candidate1 = CoordinateDic[chr][1][idx[0][0]]
                        candidate2 = CoordinateDic[chr][1][idx[0][1]]
                        list1 = range(candidate1[0],candidate1[1]+1)
                        list2 = range(candidate2[0],candidate2[1]+1)
                        listSeq = range(start, end+1)
                        if len(set(list1) & set(listSeq)) == 0 and len(set(list2) & set(listSeq)) == 0:
                            flag = "1"
                        else:
                            if len(set(list1) & set(listSeq)) > 0:
                                repetitiveElementName = candidate1[2]
                            else:
                                repetitiveElementName = candidate2[2]
                    else:
                        candidate1 = CoordinateDic[chr][1][idx[0][0]]
                        list1 = range(candidate1[0],candidate1[1]+1)
                        listSeq = range(start, end+1)
                        if len(set(list1) & set(listSeq)) == 0:
                            flag = "1"
                        else:
                            repetitiveElementName = candidate1[2]
                else:
                    flag = "1"
        outf.write('\t'.join([content[0], orginalSeq, flag, repetitiveElementName]+content[1:])+'\n')
        if flag == "1":
            trimmed_count+=1
            outFA.write('>'+content[0] + ":" + chr + ":" + str(start)+"_"+str(end) + content[2]+"\n"+orginalSeq+"\n")
        line = inf1.readline()
    inf1.close()
    outf.close()
    outFA.close()
    return trimmed_count


def predict_nmir(args, workDir, ref_db, base_names, pdUnmapped):
    htmlJS = FormatJS(workDir)
    htmlJS.openNovelmiRJSData()
    predict_start_time = time.perf_counter()
    runlogFile = Path(workDir)/"run.log"
    outlog = open(str(runlogFile),"a+")
    if not args.quiet:
        print('\nPerforming prediction of novel miRNAs...')
        print('Start to predict')
    outlog.write('\nPerforming prediction of novel miRNAs...')
    outlog.write('Start to predict')
    samtoolsBinary = Path(args.samtools_path)/"samtools " if args.samtools_path else "samtools " 
    rnafoldBinary = Path(args.RNAfold_path)/"RNAfold " if args.RNAfold_path else "RNAfold "
    species = args.organism_name
    genomRepeats = species+"_genome_repeats.pckl"
    genomFasta = species+"_genome.pckl"
    genom_miRDB = species+"_mature_"+ref_db+".fa"
    genome_repeats = Path(args.libraries_path)/args.organism_name/"annotation.Libs"/genomRepeats
    genome_fa = Path(args.libraries_path)/args.organism_name/"fasta.Libs"/genomFasta
    mature_miRNA_fa = Path(args.libraries_path)/args.organism_name/"fasta.Libs"/genom_miRDB
    nameAbbrNameDic ={'human':'hsa', 'zebrafish':'dre', 'mouse':'mmu', 'rat':'rno', 'fruitfly':'dme', 'nematode':'cel', 'hamster':'mau'}
    if species in nameAbbrNameDic.keys():
        abbrName = nameAbbrNameDic[species]
        speciesNameDic = {species:abbrName}
    else:
        speciesNameDic = {species:"CUS"}

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
        outlog.write(f'ERROR: The bowtie index file of {species}_genome.*.ebwt is not located at {indexPath}, please check it.')
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

    bwtCmdTmp = str(Path(args.bowtie_path)/"bowtie ") if args.bowtie_path else "bowtie "
    if args.bowtieVersion == "False": # That is if version is v1.3.0
        bwtCmdTmp += " -x "
    bwtBuildCmdTmp = Path(args.bowtie_path)/"bowtie-build " if args.bowtie_path else "bowtie-build " 
    samtoolsCmdTmp = Path(args.samtools_path)/"samtools " if args.samtools_path else "samtools "
    rnafoldCmdTmp = Path(args.RNAfold_path)/"RNAfold " if args.RNAfold_path else "RNAfold "
    infile = "unmapped.log"
    logFile = Path((Path(outputdir2).resolve().parents[0]))/infile
    #logFile = Path(outputdir2)/infile
    rawReadCounts={}
    filteredReadCounts={}
    with open(logFile, 'w') as outfLog:
        time1 = time.perf_counter()
        outfLog.write(f'Reading and trimming unmapped reads from Pandas Library\n')
        convert2Fasta(pdUnmapped, infile, minLength, maxLength, countCutoff, outputdir2, species, speciesNameDic, base_names,rawReadCounts,filteredReadCounts)
        #print(rawReadCounts, filteredReadCounts)
        #time.sleep(.555)
        time2 = time.perf_counter()
        outfLog.write('Execution time for reading and trimming unmapped Pandas Library took: %.4fs\n'%(time2-time1))
        outfLog.write('********************\n')
        outfLog.flush()
        for files in base_names:
            errorTrue = 0
            outfLog.write(f'Processing {files}\n')
            outfLog.write(f'**There are {str(rawReadCounts[files])} collapsed reads in the raw fasta file\n')
            outfLog.write(f'**After filtering, there are {str(filteredReadCounts[files])} reads left\n')
            outfLog.write('Mapping reads to %s genome\n'%(species))
            time3 = time.perf_counter()
            current_fname_filtered = "unmapped_mirna_"+ files +".fa"
            fileNameTemp = Path(outputdir2)/current_fname_filtered
            outfile1 = Path(outputdir2)/("unmapped_mirna_"+ files +"_vs_genome.sam")
            bwtExec = str(bwtCmdTmp) + str(genome_index) + " " + str(fileNameTemp) + " -f -n 0 --best -a --threads " + str(args.threads) + " -m " + str(mapping_loc) + " -l "+ str(seedLength) + " -S " + str(outfile1)
            bowtie = subprocess.run(str(bwtExec), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
            # SORT SAM FILE
            outfile2 = Path(outputdir2)/("unmapped_mirna_"+ files +"_vs_genome_sorted.sam") 
            samyExec = str(samtoolsCmdTmp) + " sort --threads "+ str(args.threads) + " -O sam -T sample.sort -o " + str(outfile2) + " " + str(outfile1)
            samysort = subprocess.run(str(samyExec), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
            time4 = time.perf_counter()
            outfLog.write('Mapping reads to humna genome time: %.4fs\n'%(time4-time3))
            outfLog.write('Clustering the reads based on the coordinate in the genome\n')
            cluster_file = Path(outputdir2)/(files+"_clusters.tsv")
            clusters_log = cluster_basedon_location(outfile2, overlapLenCutoff, files, outputdir2, cluster_file)
            time5 = time.perf_counter()
            outfLog.write('Clustering the reads based on the coordinate in the genome time: %.4fs\n'%(time5-time4))
            outfLog.write('Trimming the clustered sequences (preliminary filtering)\n')
            outfLog.write(f'**Before preliminary trimming, there are {clusters_log} clustered sequences\n')
            clusterTrimedFile = Path(outputdir2)/(files+"clusters_trimmed.tsv")
            clusterTrimedFile_orig_FASTA = Path(outputdir2)/(files+"_clusters_trimmed_orig.fa")
            trimedClusters_log = preTrimClusteredSeq(repEleChrCoordinateDic, cluster_file, files, clusterSeqLenCutoff, clusterTrimedFile, clusterTrimedFile_orig_FASTA)
            time6 = time.perf_counter()
            outfLog.write('Trimming the clustered sequences time: %.fs\n'%(time6-time5))
            outfLog.write(f'**After preliminary trimming, there are {trimedClusters_log} the clustered sequences\n')
            outfLog.write('Generating the orignial seqeneces of the trimmed clustered sequences\n')
            outfLog.write('Building bowtie index file and mapping reads to the clustered seqneces\n')
            time7 = time.perf_counter()
            #os.system('%s -f %s_clusters_trimmed_orig.fa %s_representative_seq >> %s 2>&1'%(bwtBuildCmdTmp, outfile2, outfile3, outLogTemp))
            outfile3 = Path(outputdir2)/(files+"_representative_seq")
            bwtBuildExec = str(bwtBuildCmdTmp) +" -f "+ str(clusterTrimedFile_orig_FASTA) + " " + str(outfile3) + " --threads " + str(args.threads) 
            try:
                bowtie = subprocess.run(str(bwtBuildExec), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
            except subprocess.CalledProcessError:
                errorTrue = 1
                pass
            # Align the total reads to the cluster sequences with exact match. '--norc' only considers matching to the forwad reference strand.
            outfile4 = Path(outputdir2)/(files+"_tmp1.sam")
            bwt2Exec = str(bwtCmdTmp) + str(outfile3) + " " + str(fileNameTemp) + " -f -n 0 --best -a --norc --threads " + str(args.threads) + " -m " + str(mapping_loc) + " -l "+ str(seedLength) + " -S " + str(outfile4)
            try:
                bowtie = subprocess.run(str(bwt2Exec), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
            except subprocess.CalledProcessError:
                errorTrue = 1
                pass
            #Generate the reads that can't align to the cluster reads.
            try:
                imperfect_FASTA=Path(outputdir2)/(files+"_imperfectMath2Cluster.fa")
                split_fasta_from_sam(outfile4, fileNameTemp, str(imperfect_FASTA))
            except IOError:
                errorTrue = 1
                outlog.write(f'No cluster sequences are generated and prediction is aborted for {files}.\n')
                #os.system('rm -r %s'%(outputdir2))
                #sys.exit(1)
            # Align the rest part of reads to the cluster sequences with up to 1 mismatches with special 5 vs 3 prime considerations.
            # The alignments in the best stratum are those having the least number of mismatches.
            # '--norc' only considers matching to the forwad reference strand.
            if errorTrue != 1:
                outfile4_tmp2 = Path(outputdir2)/(files+"_tmp2.sam")
                bwt3Exec = str(bwtCmdTmp) + str(outfile3) + " " + str(imperfect_FASTA) + " -f -n 1 -l 15 -5 1 -3 3 --best --strata -a --norc --threads " + str(args.threads) + " -S " + str(outfile4_tmp2)
                bowtie = subprocess.run(str(bwt3Exec), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
                # Combine the aligned result of the two type of reads: perfect matched reads and imperfect matched reads.
                combined_Sam = Path(outputdir2)/(files+".sam")  
                combineSam(str(outfile4), str(outfile4_tmp2), str(combined_Sam))
                time8 = time.perf_counter()
                outfLog.write('Building bowtie index file and mapping reads to the clustered seqneces time: %.1fs\n'%(time8-time7))
                outfLog.write('Decorating bowtie output sam files\n')
                outfile_modifiedSam = Path(outputdir2)/(files+"_modified.sam")
                outfile_modSam_RepSeq = Path(outputdir2)/(files+"_RepSeq_modified.sam")
                decorateSam(str(combined_Sam), str(fileNameTemp), str(outfile_modifiedSam), str(clusterTrimedFile_orig_FASTA)) #Decorate SAM handles two different inputs and produce two outputs!!!
                decorateSam(str(outfile2), str(fileNameTemp), str(outfile_modSam_RepSeq))
                outfLog.write('Decorating is done\n')
                outfileSelectTSV = Path(outputdir2)/(files+"_selected.tsv")
                outfileRevKeptTSV = Path(outputdir2)/(files+"_selected_reverseKept.tsv")
                parse_refine_sam(str(outfile_modifiedSam), str(outfileSelectTSV), str(outfileRevKeptTSV))
                #os.system('rm %s.bam %s.sam %s.bai %s.bam'%(outfile1, outfile1, outfile2, outfile2))
                os.system('sort -k6,6 -k1,1 %s > %s'%(str(outfileSelectTSV), str(Path(outputdir2)/(files+"_modified_selected_sorted.tsv"))))
                os.system('sort -k6,6 -k1,1 %s > %s'%(str(outfileRevKeptTSV), str(Path(outputdir2)/(files+"_modified_selected_reverseKept_sorted.tsv"))))
                # Trimming the clustered seuences based on the alligned results of all the reads (secondary filtering)
                #generate_featureFiles(outfile4+'_modified_selected_sorted.tsv', chrSeqDic, chrSeqLenDic, miRNAchrCoordivateDic, exactmiRNASeqDic)
                generate_featureFiles(str(Path(outputdir2)), files, chrSeqDic, chrSeqLenDic, miRNAchrCoordivateDic, exactmiRNASeqDic)
                outfLog.write('feature file is generated\n')
                outfLog.write('Generating the precusors based on the clustered sequences\n')
                time9 = time.perf_counter()
                precursor_count = get_precursors(str(Path(outputdir2)), files, chrSeqDic)
                time10 = time.perf_counter()
                outfLog.write('Generating the precusors based on the clustered sequences time: %.1fs\n'%(time10-time9))
                outfLog.write(f'**There are {str(precursor_count)} precusors generated\n')
                outfLog.write('Folding precursors with RNAfold\n')
                time11 = time.perf_counter()
                infile_pre = str(Path(outputdir2)/(files+"_precursor.fa"))
                outfile_str = str(Path(outputdir2)/(files+"_precursor_tmp.str"))
                rnafld_exec = str(rnafoldCmdTmp) + " " + str(infile_pre) + " --noPS --noLP > " + str(outfile_str)
                rnafldRun = subprocess.run(str(rnafld_exec), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
                strFileOut= str(Path(outputdir2)/(files+"_precursor.str"))
                renameStrFile(infile_pre, outfile_str, strFileOut)
                time12 = time.perf_counter()
                outfLog.write('Folding precursors with RNAfold time: %.1fs\n'%(time12-time11))
                screen_precusor_candidates(str(Path(outputdir2)), files, str(Path(outputdir2)/(files+"_features.tsv")), strFileOut, str(rnafoldCmdTmp))
                outfLog.write('********************\n')
                outfLog.flush()
                modelDirTmp = Path(__file__).resolve().parents[1]
                modelDir = str(Path(modelDirTmp)/'models')
                fileToPredict = Path(outputdir2)/(files+'_updated_stableClusterSeq_15.tsv')
                preprocess_featureFiles(str(Path(outputdir2)), files, fileToPredict, str(Path(modelDir)/'total_features_namelist.txt'))
                speciesType = args.organism_name
                if args.organism_name == "human" or args.organism_name == "mouse":
                    mf = str(Path(modelDir)/(speciesType+'_svc_model.pkl'))
                else:
                    mf = str(Path(modelDir)/"others_svc_model.pkl")
                model_predict(str(outputdir2), files, mf)
                #model_predict(str(outputdir2), files, str(Path(modelDir)/(speciesType+'_svc_model.pkl')))
                novelmiRNALListFile = str(Path(outputdir2)/(files+'_novel_miRNAs_miRge2.0.csv'))
                featureFile = fileToPredict
                clusterFile = str(Path((outputdir2)/(files+'_cluster.txt')))
                write_novel_report(novelmiRNALListFile, featureFile, clusterFile, str(rnafoldCmdTmp), str(Path(outputdir2)), files)
    if errorTrue ==1: 
        print(f'No cluster sequences are generated and prediction is aborted.')
    predict_end_time = time.perf_counter()
    os.system('rm -r %s'%(outputdir2))
    htmlJS.closeNovelmiRJSData()
    if not args.quiet:
        print('Prediction of novel miRNAs Completed (%.2f sec)'%(predict_end_time - predict_start_time))
    outlog.write('Prediction of novel miRNAs Completed (%.2f sec)'%(predict_end_time - predict_start_time))
    outlog.close()
