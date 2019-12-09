#!/usr/bin/env python

# CHANGE LOG:
# 2-17-2017 YL conduct calling Bowtie to mapping.
# 4-14-2017 YL reorgonize the code and transform the scattered script files into unify classes and functions.
# 6-23-2017 YL modify the input interface and integrate the predictive model.
# 9-22-2017 YL add a '-ai' option.
# 1-01-2018 YL add a '-gff' option.
# 2-22-2018 YL add a '-ex' option.
# 5-06-2018 YL rebuild the miRNA libraries from the newly released miRBase v22 and MirGeneDB v2.0 and tRNA libraries (including mature tRNA and primary tRNA).
# 6-26-2018 YL added the function of detecting tRFs

__author__ = 'Marc Halushka, Yin Lu'
__copyright__ = 'Copyright 2018, Johns Hopkins University'
__credits__ = ['Marc Halushka', 'Yin Lu']
__license__ = 'GPL'
__version__ = '2.0'
__maintainer__ = 'Yin Lu'
__email__ = 'ylu61@jhmi.edu'

# System imports
import os
import sys
import argparse
import subprocess
import time
import re
import cPickle
import shutil
import commands
import warnings

from distutils.spawn import find_executable
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from scipy.spatial import cKDTree

# Class imports
from mirge.classes.readCluster import ReadCluster
from mirge.classes.readPrecusor import ReadPrecusor

# Fuction modules imports
from mirge.utils.parseArgument import *
from mirge.utils.trim_file import *
from mirge.utils.quantReads import *
from mirge.utils.runAnnotationPipeline import *
from mirge.utils.summarize import *
from mirge.utils.miRNAmerge import *
from mirge.utils.filter import *
from mirge.utils.generateReport import *
from mirge.utils.writeDataToCSV import *
from mirge.utils.writeDataToCSV import addDashNew
from mirge.utils.convert2Fasta import *
from mirge.utils.cluster_basedon_location import *
from mirge.utils.preTrimClusteredSeq import *
from mirge.utils.generate_Seq import *
from mirge.utils.processSam import *
from mirge.utils.get_precursors import *
from mirge.utils.generate_featureFiles import *
from mirge.utils.fileExist import *
from mirge.utils.renameStrFile import *
from mirge.utils.screen_precusor_candidates import *
from mirge.utils.preprocess_featureFiles import *
from mirge.utils.model_predict import *
from mirge.utils.write_novel_report import *
from mirge.utils.extractPreMiRName import *

# If infFile is mapped.csv, the testMode will be turned on automatically or manually.
# test mode will call cluster_analysis.py.
# testMode = True

def main():
	warnings.filterwarnings("ignore")
	# Parse the input arguments
	args = parseArgument()
	# Preprocess the inout arguments
	sampleListTmp = args.sampleList
	dir = os.path.abspath(args.output_dir)
	if args.miRNA_database.lower() in ['mirbase', 'mirgenedb']:
		if args.miRNA_database.lower() == 'mirbase':
			miRNA_database = 'miRBase'
		else:
			miRNA_database = 'MirGeneDB'
	else:
		print >> sys.stderr, "The value of parameter '-d' is invalid. Please check it"
		sys.exit(1)
	bowtieBinary = args.bowtieBinary
	# Check the binarys directory to make sure it can work
	executable1 = find_executable(os.path.join(bowtieBinary, 'bowtie'))
	if not executable1:
		print >> sys.stderr, "bowtie can't be found in %s. Please check it."%(bowtieBinary)
		sys.exit(1)
	if sys.argv[1] == 'predict':
		samtoolsBinaryTmp = args.samtoolsBinary
		rnafoldBinaryTmp = args.rnafoldBinary
		executable2 = find_executable(os.path.join(samtoolsBinaryTmp, 'samtools'))
		executable3 = find_executable(os.path.join(rnafoldBinaryTmp, 'RNAfold'))
		if (not executable2) or (not executable3):
			if not executable2:
				print >> sys.stderr, "samtools can't be found in %s. Please check it."%(samtoolsBinaryTmp)
			else:
				pass
			if not executable3:
				print >> sys.stderr, "RNAfold can't be found in %s. Please check it."%(rnafoldBinaryTmp)
			sys.exit(1)

	libraryPath = args.libraryPath
	species = args.species
	indexPath = os.path.join(libraryPath, species, 'index.Libs')
	miRNA_fa = os.path.join(libraryPath, species, 'fasta.Libs', species+'_mirna_SNP_pseudo_'+miRNA_database+'.fa')
	mergeLibFile = os.path.join(libraryPath, species, 'annotation.Libs', species+'_merges_'+miRNA_database+'.csv')

	mirMergedNameDic = {}
	# key is the name of miRNAs (including SNPs), and value is the merged name.
	with open(mergeLibFile, 'r') as inf:
		for line in inf:
			line_content = line.strip().split(',')
			for item in line_content[1:]:
				mirMergedNameDic.update({item:line_content[0]})

	canoRatio =args.canoRatio
	adapter = args.adapter
	if adapter == 'illumina':
		adapter = 'TGGAATTCTCGGGTGCCAAGGAACTCCAG'
	elif adapter == 'ion':
		adapter = '11'
	else:
		pass

	phred = args.phred64
	spikeIn = args.spikeIn
	trimmed_collapsed_fa = args.trimmed_collapsed_fa
	isomirDiff = args.diff_isomirs
	numCPU = args.cpu
	a_to_i = args.a_to_i
	miRNAs_in_repetitive_element = os.path.join(libraryPath, species, 'annotation.Libs', species+'_miRNAs_in_repetitive_element_'+miRNA_database+'.csv')
	removedMiRNA_ai_List = []
	if a_to_i:
		if os.path.isfile(miRNAs_in_repetitive_element):
			with open(miRNAs_in_repetitive_element, 'r') as inf:
				for line in inf:
					if line.strip().split(',')[0] not in removedMiRNA_ai_List:
						removedMiRNA_ai_List.append(line.strip().split(',')[0])
		#else:
			#print >> sys.stderr,'please input the file that contains the miRNA name located at repetetive element regions  when -ai is turned on.'
			#sys.exit(1)

	gff_output = args.gff_output
	miRNA_coordinate = os.path.join(libraryPath, species, 'annotation.Libs', species+'_'+miRNA_database+'.gff3')
	if gff_output:
		miRNamePreNameDic = extractPreMiRName(miRNA_coordinate, miRNA_database)
		isomiRContentDic = {}
		# the keys of isomiRContentDic are miRNA and isomiR sequences. The value is a dictionary and the keys of the dictionary are 'miRName', 'preMiRName', 'type', 'pre_start',
		# 'pre_end', 'strand', 'variant', 'cigar', 'expression'
	else:
		miRNamePreNameDic = None
		isomiRContentDic = None

	trf_output = args.trf_output
	if trf_output and species != "human":
		print >> sys.stderr, "tRF detection is only supported for the species of human. Please check it."
		sys.exit(1)
	if trf_output:
		# trna_stru_file is the structure of mature tRNA
		trna_stru_file = os.path.join(libraryPath, species, 'annotation.Libs', species+'_trna.str')
		trnaStruDic = {}
		try:
			with open(trna_stru_file, 'r') as inf:
				line = inf.readline()
				while line != '':
					trnaName = line.strip()[1:]
					line = inf.readline()
					trnaSeq = line.strip()
					line = inf.readline()
					trnaStru = line.strip()
					# 1-based position
					anticodonStart = trnaStru.index('XXX')+1
					anticodonEnd = anticodonStart+2
					trnaStruDic.update({trnaName:{'seq':trnaSeq, 'stru':trnaStru, 'anticodonStart':anticodonStart, 'anticodonEnd':anticodonEnd}})
					line = inf.readline()
		except IOError:
			print >> sys.stderr, "%s does not exsit. Please check it."%(trna_stru_file)
			sys.exit(1)
		
		trna_aa_anticodon_file = os.path.join(libraryPath, species, 'annotation.Libs', species+'_trna_aminoacid_anticodon.csv')
		trnaAAanticodonDic = {}
		try:
			with open(trna_aa_anticodon_file, 'r') as inf:
				for line in inf:
					contentTmp = line.strip().split(',')
					trnaAAanticodonDic.update({contentTmp[0]:{'aaType':contentTmp[1], 'anticodon':contentTmp[2]}})
		except IOError:
			print >> sys.stderr, "%s does not exsit. Please check it."%(trna_aa_anticodon_file)
			sys.exit(1)
		
		trfContentDic = {}
		# the keys of trfContentDic are tRF sequences. The value is a dictionary and the keys of the dictionary are 'uid', 'count', 'RPM', 'location'
		trna_duplicated_list_file = os.path.join(libraryPath, species, 'annotation.Libs', species+'_trna_deduplicated_list.csv')
		duptRNA2UniqueDic = {}
		try:
			with open(trna_duplicated_list_file, 'r') as inf:
				line = inf.readline()
				line = inf.readline()
				while line != '':
					contentTmp = line.strip().split(',')
					for item in contentTmp[1].split('/'):
						duptRNA2UniqueDic.update({item.strip():contentTmp[0].strip()})
					line = inf.readline()
		except IOError:
			print >> sys.stderr, "%s does not exsit. Please check it."%(trna_duplicated_list_file)
			sys.exit(1)
		
		# Load predefined tRF loci information file
		tRNAtrfDic = {}
		tRF_infor_file = os.path.join(libraryPath, species, 'annotation.Libs', species+'_tRF_infor.csv')
		try:
			with open(tRF_infor_file, 'r') as  inf:
				line = inf.readline()
				line = inf.readline()
				while line != '':
					content = line.strip().split(',')
					tRNAName = content[0].split('_Cluster')[0]
					tRNAClusterName = content[0]
					seq = content[4]
					tRNAlength = len(content[5])
					start = int(content[3].split('-')[0])
					end = int(content[3].split('-')[1])
					if tRNAName not in tRNAtrfDic.keys():
						tRNAtrfDic.update({tRNAName:{}})
					tRNAtrfDic[tRNAName].update({addDashNew(seq, tRNAlength, start, end):tRNAClusterName})
					line = inf.readline()
		except IOError:
			print >> sys.stderr, "%s does not exsit. Please check it."%(tRF_infor_file)
			sys.exit(1)
		
		# Load predifined tRF merged file
		trfMergedNameDic = {}
		trfMergedList = []
		tRF_merge_file = os.path.join(libraryPath, species, 'annotation.Libs', species+'_tRF_merges.csv')
		try:
			with open(tRF_merge_file, 'r') as inf:
				for line in inf:
					tmp = line.strip().split(',')
					mergedName = tmp[0]
					trfMergedList.append(mergedName)
					for item in tmp[1].split('/'):
						trfMergedNameDic.update({item:mergedName})
		except IOError:
			print >> sys.stderr, "%s does not exsit. Please check it."%(tRF_merge_file)
			sys.exit(1)
	else:
		trnaStruDic = None
		trnaAAanticodonDic = None
		trfContentDic = None
		duptRNA2UniqueDic = None
		tRNAtrfDic = None
		trfMergedNameDic = None
		trfMergedList = None
	
	# Check index files 
	for type in ['mirna_'+miRNA_database, 'hairpin_'+miRNA_database, 'mrna', 'mature_trna', 'pre_trna', 'snorna', 'rrna', 'ncrna_others']:
		if os.path.isfile(os.path.join(indexPath, species+'_'+type+'.1.ebwt')):
			pass
		else:
			print 'The bowtie index file of %s_%s.*.ebwt is not located at %s, please check it.'%(species, type, indexPath)
			sys.exit(1)

	mirna_index = os.path.join(indexPath, species+'_mirna_'+miRNA_database)
	hairpin_index = os.path.join(indexPath, species+'_hairpin_'+miRNA_database)
	mrna_index = os.path.join(indexPath, species+'_mrna')
	mature_tRNA_index = os.path.join(indexPath, species+'_mature_trna')
	pre_tRNA_index = os.path.join(indexPath, species+'_pre_trna')
	snoRNA_index = os.path.join(indexPath, species+'_snorna')
	rRNA_index = os.path.join(indexPath, species+'_rrna')
	ncrna_others_index = os.path.join(indexPath, species+'_ncrna_others')
	genome_index = os.path.join(indexPath, species+'_genome')
	if spikeIn:
		spikeIn_index = os.path.join(indexPath, species+'_spike-in')
	else:
		spikeIn_index = None

	# set the outout directory
	tStamp = time.strftime('%Y-%m-%d_%H-%M-%S',time.localtime(time.time()))
	outputdir = os.path.join(dir,'miRge.'+tStamp)
	os.mkdir(outputdir)
	
	# Read the sampleListFile.
	dir_sample = os.path.split(os.path.abspath(sampleListTmp[0]))[0]
	sampleListRaw = []
	if all(os.path.basename(sampleTmp).split('.')[-1] == 'fastq' or '.'.join(os.path.basename(sampleTmp).split('.')[-2:]) == 'fastq.gz' for sampleTmp in sampleListTmp):
		for sampleTmp in sampleListTmp:
			if os.path.isfile(os.path.abspath(sampleTmp)):
				pass
			else:
				print "%s can't be found in current directory. Please check it."%(sampleTmp)
				sys.exit(1)
		sampleListRaw = sampleListTmp
	elif (len(sampleListTmp) == 1 and (os.path.basename(sampleListTmp[0]).split('.')[-1] != 'fastq' or '.'.join(os.path.basename(sampleListTmp[0]).split('.')[-2:]) != 'fastq.gz')):
		try:
			with open(sampleListTmp[0], 'r') as inf:
				for line in inf:
					if line.strip() not in sampleListRaw:
						if os.path.isfile(os.path.abspath(line.strip())):
							sampleListRaw.append(line.strip())
						else:
							print '%s cannot be found, please check the path of the sample file.'%(line.strip())
							sys.exit(1)
		except IOError,e:
			print '%s is not a file, please check it.'%(sampleListTmp[0])
			sys.exit(1)
	else:
		print "The format of input argument '-s' is wrong, please check it."
		sys.exit(1)
	
	### Arguments preprocession is done!
	# Perform annotation analysis
	seqDic = {}
	#  key is nucleotide sequence, each record is a dictionary with
	# 'annot' a boolean array where element [0] is overall matched status and each element [i] the matched status from each respective round of alignment
	# 'quant' a integer array of the counts from each sample
	mirDic = {} 
	# key is annot string (the name of the miRNA) from miRNA library records 
	mirNameSeqDic = {}
	# key is annot string (the name of the miRNA, including merged miRNAs), value is the sequence.
	readLengthDic = {}
	# key is length
	logDic = {}
	logDic.update({'quantStats':[]})
	logDic.update({'annotStats':[]})
	# two primary keys 'quantStats' and 'annotStats'
	if spikeIn:
		annotNameList = ['exact miRNA', 'hairpin miRNA', 'mature tRNA', 'primary tRNA', 'snoRNA', 'rRNA', 'ncrna others', 'mRNA', 'isomiR miRNA', 'spike-in']
	else:
		annotNameList = ['exact miRNA', 'hairpin miRNA', 'mature tRNA', 'primary tRNA', 'snoRNA', 'rRNA', 'ncrna others', 'mRNA', 'isomiR miRNA']
	
	# perform quantitation analysis.
	time_0 = time.time()
	sampleList_tmp = [os.path.basename(item) for item in sampleListRaw]
	sampleList = []
	for item in sampleList_tmp:
		if item[-3:] == '.gz':
			sampleList.append(item[:-3])
		else:
			sampleList.append(item)
	for i, sampleFile in enumerate(sampleListRaw):
		sampleFileName = os.path.basename(sampleFile)
		if sampleFileName[-3:] == '.gz':
			sampleFileName = sampleFileName[:-3]
		print 'Performing quantitation analysis of %s...'%(sampleFileName)
		dicRmp = {}
		dicRmp.update({'filename' : sampleFileName})
		time_1 = time.time()
		#samplePrefix = os.path.join(dir, sampleFile)
		samplePrefix = os.path.abspath(sampleFile)
		cleanedReads = os.path.join(outputdir, os.path.splitext(sampleFileName)[0]+'.trim.fastq')
		# Trim sequence file
		phred, processed_count, kept_count = trim_file(samplePrefix, adapter, cleanedReads, int(numCPU))
		phred64 = phred==64
		dicRmp.update({'totalReads':processed_count}) 
		dicRmp.update({'trimmedReads':kept_count}) 
		time_2 = time.time()
		dicRmp.update({'cpuTime-trim' : time_2-time_1})
		# Unique sequence file
		quantReads(cleanedReads, seqDic, readLengthDic, len(sampleList), i, sampleList, trimmed_collapsed_fa, spikeIn)
		time_3 = time.time()
		dicRmp.update({'cpuTime-uniq' : time_3-time_2})
		logDic['quantStats'].append(dicRmp)
		print 'It takes: %.2fs'%(time_3-time_1)
	
	# perform annotation
	print '\nPerforming annotation for all of the collasped sequences...'
	time_4 = time.time()
	#numCPU_new = '1'
	runAnnotationPipeline(bowtieBinary, seqDic, numCPU, phred64, annotNameList, outputdir, logDic, mirna_index, hairpin_index, mature_tRNA_index, pre_tRNA_index, snoRNA_index, rRNA_index, ncrna_others_index, mrna_index, spikeIn, spikeIn_index, gff_output, miRNamePreNameDic, isomiRContentDic, miRNA_database, trf_output, trnaStruDic, trfContentDic, sampleList)
	time_5 = time.time()
	print 'All annotation cycles completed (%.2f sec).\n'%(time_5-time_4)

	# Summarize and tabulate the results
	print "Summarizing and tabulating results..."
	summarize(seqDic, sampleList, logDic, mirDic, mirna_index, outputdir, spikeIn, bowtieBinary)

	miRNAmerge(mergeLibFile, sampleList, mirDic, miRNA_fa, mirNameSeqDic)

	filter(mirDic, sampleList, logDic, canoRatio)

	generateReport(outputdir, sampleList, readLengthDic, logDic, annotNameList, seqDic, spikeIn)
	
	writeDataToCSV(outputdir, annotNameList, sampleList, isomirDiff, a_to_i, logDic, seqDic, mirDic, mirNameSeqDic, mirMergedNameDic, bowtieBinary, genome_index, numCPU, phred64, removedMiRNA_ai_List, spikeIn, gff_output, isomiRContentDic, miRNA_database, trf_output, trfContentDic, trnaStruDic, pre_tRNA_index, duptRNA2UniqueDic, trnaAAanticodonDic, tRNAtrfDic, trfMergedNameDic, trfMergedList)
	time_6 = time.time()
	print 'Summary Complete (%.2f sec)'%(time_6-time_5)
	print 'Annotation of miRge2.0 Completed (%.2f sec)'%(time_6-time_0)

	if sys.argv[1] == 'predict':
		time_7 = time.time()
		print '\nPerforming prediction of novel miRNAs...'
		print 'Start to predict'
		samtoolsBinary = args.samtoolsBinary
		rnafoldBinary = args.rnafoldBinary
		genome_repeats = os.path.join(libraryPath, species, 'annotation.Libs', species+'_genome_repeats.pckl')
		genome_fa = os.path.join(libraryPath, species, 'fasta.Libs', species+'_genome.pckl')
		mature_miRNA_fa = os.path.join(libraryPath, species, 'fasta.Libs', species+'_mature_'+miRNA_database+'.fa')
		nameAbbrNameDic ={'human':'hsa', 'zebrafish':'dre', 'mouse':'mmu', 'rat':'rno', 'fruitfly':'dme', 'nematode':'cel'}
		abbrName = nameAbbrNameDic[species]
		speciesNameDic = {species:abbrName}
		wideSampleListFile = args.wideSampleListFile
		minLength = args.minLength
		maxLength = args.maxLength
		countCutoff = args.countCutoff
		mapping_loc = args.mapping_loc
		seedLength = args.seedLength
		# Here, the optimal overlapLenCutoff should be set to be 14, based on the trial of 8, 10 , 12
		overlapLenCutoff = args.overlapLenCutoff
		clusterSeqLenCutoff = args.clusterSeqLenCutoff
		if os.path.isfile(os.path.join(indexPath, species+'_genome.1.ebwt')):
			pass
		else:
			print 'The bowtie index file of %s_genome.*.ebwt is not located at %s, please check it.'%(species, indexPath)
			sys.exit(1)
		genome_index = os.path.join(indexPath, species+'_genome')
		# In default, wideSampleListFile is 'none'(Null) so that all of the samples in the temparary output file 'mapped.csv' or 'unmapped.csv' will be all analyzed.
		# If wideSampleListFile is designated, only the samples in 'mapped.csv' or 'unmapped.csv' that are covered by wideSampleListFile will be analyzed for novel miRNA predition.
		if wideSampleListFile == 'none':
			wideSampleList = None
		else:
			try:
				with open(wideSampleListFile, 'r') as inf:
					wideSampleList = [item.strip() for item in inf.readlines() if item.strip() != '']
			except IOError,e:
				print '%s is not a file, please check it.'%(wideSampleListFile)
				sys.exit(1)
		#Load the coordinate file of repeat region in Genome which is GRCH38_genome_repeats_sorted.pckl'
		try:
			with open(genome_repeats, 'rb') as f:
				repEleChrCoordinateDic = cPickle.load(f)
		except IOError:
			repEleChrCoordinateDic = {}
		#Load the genome fasta file: human_genome.pckl
		with open(genome_fa, 'r') as f:
			chrSeqDic = cPickle.load(f)
		chrSeqLenDic = {}
		for key in chrSeqDic.keys():
			chrSeqLenDic.update({key:len(chrSeqDic[key])})
		#Load genome miRNA fasta file
		exactmiRNASeqDic = {}
		for record in SeqIO.parse(mature_miRNA_fa, "fasta"):
			exactmiRNASeqDic.update({record.id:str(record.seq)})
		#Load the genome miRNA coordinate file
		miRNAchrCoordivateDic = {}
		with open(miRNA_coordinate,"r") as inf1:
			for line1 in inf1:
				if line1[0] != "#":
					content = line1.strip().split("\t")
					if content[2] == "miRNA":
						chr = content[0]
						if miRNA_database == 'miRBase':
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
	if sys.argv[1] == 'predict':
		#for fileType in ['mapped.csv', 'unmapped.csv']:
		for fileType in ['unmapped.csv']:
			infile = os.path.join(outputdir, fileType)
			# Set the output directory of 'mapped.csv' or 'unmapped.csv'
			outputdir2 = os.path.join(outputdir,fileType.split('.')[0]+'_tmp')
			os.mkdir(outputdir2)

			bwtCmdTmp = os.path.join(bowtieBinary, 'bowtie')
			bwtBuildCmdTmp = os.path.join(bowtieBinary, 'bowtie-build')
			samtoolsCmdTmp = os.path.join(samtoolsBinary, 'samtools')
			rnafoldCmdTmp = os.path.join(rnafoldBinary, 'RNAfold')
			logFile = os.path.join(outputdir2, os.path.splitext(infile)[0]+'.log')
			with open(os.path.join(outputdir2, os.path.splitext(infile)[0]+'.log'), 'w') as outfLog:
				time1 = time.time()
				outfLog.write('Reading and trimming %s\n'%(infile))
				convert2Fasta(infile, minLength, maxLength, countCutoff, outputdir2, species, speciesNameDic, spikeIn)
				time2 = time.time()
				outfLog.write('Reading and trimming %s time: %.1fs\n'%(infile, time2-time1))
				outfLog.write('********************\n')
				outfLog.flush()
				# remove the samples if they are not in wideSampleList
				if wideSampleList == None:
					pass
				else:
					for file in os.listdir(outputdir2):
						if file not in ['mapped_mirna.fa', 'mapped_mirna_raw.fa', 'mapped_nonMirna.fa', 'mapped_nonMirna_raw.fa', 'mapped.log',
							'unmapped_mirna_raw.fa', 'unmapped_mirna.fa', 'unmapped.log']:
							if not fileExist(file, wideSampleList):
								os.remove(os.path.join(outputdir2, file))

				for file in os.listdir(outputdir2):
					fileName = os.path.join(outputdir2,file)
					if os.path.isfile(fileName) and fileExist(file, sampleList) and ('_raw' not in file):
						#print file
						outputdirTemp = os.path.join(outputdir2,os.path.splitext(file)[0]+'_tmp')
						os.mkdir(outputdirTemp)
						shutil.move(fileName, outputdirTemp)
						fileNameTemp = os.path.join(outputdirTemp,file)
						outfLog.write('Processing %s\n'%(file))
						outfLog.write('**There are %d collapsed reads in the raw fasta file\n'%((int(commands.getstatusoutput('cat %s_raw.fa|wc -l'%(os.path.splitext(fileName)[0]))[1]))/2))
						shutil.move(os.path.splitext(fileName)[0]+'_raw.fa', outputdirTemp)
						outfLog.write('**After filtering, there are %d reads left\n'%((int(commands.getstatusoutput('cat %s|wc -l'%(fileNameTemp))[1]))/2))
						outLogTemp = os.path.join(outputdirTemp, os.path.splitext(file)[0]+'.log')
						outfile1 = os.path.join(outputdirTemp, os.path.splitext(file)[0]+'_vs_genome')
						outfile2 = os.path.join(outputdirTemp, os.path.splitext(file)[0]+'_vs_genome_sorted')
						outfile3 = os.path.join(outputdirTemp, os.path.splitext(file)[0])
						outfile4 = os.path.join(outputdirTemp, os.path.splitext(file)[0]+'_vs_representative_seq')

						outfLog.write('Mapping reads to %s genome\n'%(species))
						time3 = time.time()
						#print '%s %s %s -f -n 0 -m %s -l %s -S -a --best %s.sam >> %s 2>&1'%(bwtCmdTmp, genome_index, fileNameTemp, mapping_loc, seedLength, outfile1, outLogTemp)
						os.system('%s %s %s -f -n 0 -m %s -l %s -S -a --best %s.sam >> %s 2>&1'
								%(bwtCmdTmp, genome_index, fileNameTemp, mapping_loc, seedLength, outfile1, outLogTemp))
						os.system('%s view -Sb %s.sam > %s.bam 2>%s'%(samtoolsCmdTmp, outfile1, outfile1, outLogTemp))
						os.system('%s sort %s.bam -o %s.bam'%(samtoolsCmdTmp, outfile1, outfile2))
						os.system('%s index %s.bam %s.bai'%(samtoolsCmdTmp, outfile2, outfile2))
						os.system('%s view -h %s.bam > %s.sam'%(samtoolsCmdTmp, outfile2, outfile2))
						time4 = time.time()
						outfLog.write('Mapping reads to humna genome time: %.1fs\n'%(time4-time3))

						outfLog.write('Clustering the reads based on the coordinate in the genome\n')
						cluster_basedon_location(outfile2+'.sam', overlapLenCutoff)
						time5 = time.time()
						outfLog.write('Clustering the reads based on the coordinate in the genome time: %.1fs\n'%(time5-time4))
						os.system('rm %s.bam %s.sam %s.bai %s.bam'%(outfile1, outfile1, outfile2, outfile2))

						#Trim the clustered sequences (preliminary filtering)
						
						outfLog.write('Trimming the clustered sequences (preliminary filtering)\n')
						outfLog.write('**Before preliminary trimming, there are %d the clustered sequences\n'%(int(commands.getstatusoutput('cat %s_clusters.tsv|wc -l'%(outfile2))[1])-1))
						preTrimClusteredSeq(repEleChrCoordinateDic, outfile2+'_clusters.tsv', clusterSeqLenCutoff)
						time6 = time.time()
						outfLog.write('Trimming the clustered sequences time: %.1fs\n'%(time6-time5))

						#Generate the orignial seqeneces of the trimmed clustered sequences.
						generate_Seq(outfile2+'_clusters_trimmed.tsv')
						outfLog.write('**After preliminary trimming, there are %d the clustered sequences\n'%(int(commands.getstatusoutput('cat %s_clusters_trimmed_orig.fa|wc -l'%(outfile2))[1])/2))
						outfLog.write('Generating the orignial seqeneces of the trimmed clustered sequences\n')

						outfLog.write('Building bowtie index file and mapping reads to the clustered seqneces\n')
						time7 = time.time()
						os.system('%s -f %s_clusters_trimmed_orig.fa %s_representative_seq >> %s 2>&1'
								%(bwtBuildCmdTmp, outfile2, outfile3, outLogTemp))

						# Align the total reads to the cluster sequences with exact match. '--norc' only considers matching to the forwad reference strand.
						os.system('%s %s %s -f -n 0 -l %s -S -a --best --norc %s_tmp1.sam >> %s 2>&1'
								%(bwtCmdTmp, os.path.join(outputdirTemp, os.path.splitext(file)[0]+'_representative_seq'), fileNameTemp, seedLength, outfile4, outLogTemp))

						# Generate the reads that can't align to the cluster reads.
						try:
							split_fasta_from_sam(outfile4+'_tmp1.sam', fileNameTemp)
						except IOError:
							print 'No cluster sequences are generated and prediction is aborted.'
							os.system('rm -r %s'%(outputdir2))
							sys.exit(1)

						# Align the rest part of reads to the cluster sequences with up to 1 mismatches with special 5 vs 3 prime considerations.
						# The alignments in the best stratum are those having the least number of mismatches.
						# '--norc' only considers matching to the forwad reference strand.
						os.system('%s %s %s_imperfectMath2Cluster.fa -f -n 1 -l 15 -5 1 -3 3 -S -a --best --strata --norc %s_tmp2.sam >> %s 2>&1'%(bwtCmdTmp, os.path.join(outputdirTemp, os.path.splitext(file)[0]+'_representative_seq'), fileNameTemp[:-3], outfile4, outLogTemp))

						# Combine the aligned result of the two type of reads: perfect matched reads and imperfect matched reads.
						combineSam(outfile4+'_tmp1.sam', outfile4+'_tmp2.sam')

						time8 = time.time()
						outfLog.write('Building bowtie index file and mapping reads to the clustered seqneces time: %.1fs\n'%(time8-time7))

						outfLog.write('Decorating bowtie output sam files\n')
						decorateSam(outfile4+'.sam', fileNameTemp, outfile2+'_clusters_trimmed_orig.fa')
						decorateSam(outfile2+'.sam', fileNameTemp)
						outfLog.write('Decorating is done\n')
						parse_refine_sam(outfile4+'_modified.sam')

						os.system('sort -k6,6 -k1,1 %s_modified_selected.tsv > %s_modified_selected_sorted.tsv'%(outfile4, outfile4))
						os.system('sort -k6,6 -k1,1 %s_modified_selected_reverseKept.tsv > %s_modified_selected_reverseKept_sorted.tsv'%(outfile4, outfile4))

						# Trimming the clustered seuences based on the alligned results of all the reads (secondary filtering)
						generate_featureFiles(outfile4+'_modified_selected_sorted.tsv', chrSeqDic, chrSeqLenDic, miRNAchrCoordivateDic, exactmiRNASeqDic)
						outfLog.write('feature file is generated\n')
						outfLog.write('Generating the precusors based on the clustered sequences\n')
						time9 = time.time()
						get_precursors(outfile4+'_modified_selected_sorted_features.tsv', chrSeqDic)
						time10 = time.time()
						outfLog.write('Generating the precusors based on the clustered sequences time: %.1fs\n'%(time10-time9))
						outfLog.write('**There are %d precusors generated\n'%(int(commands.getstatusoutput('cat %s_modified_selected_sorted_features_precusor.fa|wc -l'%(outfile4))[1])/2))
						
						# Folding precursors with RNAfold
						outfLog.write('Folding precursors with RNAfold\n')
						time11 = time.time()
						os.system('%s < %s_modified_selected_sorted_features_precusor.fa --noPS --noLP > %s_modified_selected_sorted_features_precusor_tmp.str'%(rnafoldCmdTmp, outfile4, outfile4))
						renameStrFile(outfile4+'_modified_selected_sorted_features_precusor.fa', outfile4+'_modified_selected_sorted_features_precusor_tmp.str', outfile4+'_modified_selected_sorted_features_precusor.str')
						time12 = time.time()
						outfLog.write('Folding precursors with RNAfold time: %.1fs\n'%(time12-time11))
						# Computing randfold p-values
						#os.system('randfold -s %s_modified_selected_sorted_features_precusor.fa 50 > %s_modified_selected_sorted_features_precusor.rand'%(outfile4, outfile4))
						# Compute the structual features of the precusor for each cluter sequence.
						# Here, the reatained one putative precusor should be the optimal one filtered by some criteria. 
						#for i in range(0, 21, 5):
						for i in [15]:
							screen_precusor_candidates(outfile4+'_modified_selected_sorted_features.tsv', outfile4+'_modified_selected_sorted_features_precusor.str', i, rnafoldCmdTmp)
						outfLog.write('********************\n')
						outfLog.flush()
						# Invoke the predictive model to predict
						
						if fileType == 'unmapped.csv':
							#PopenTmp = subprocess.Popen(['which', sys.argv[0]], stdout=subprocess.PIPE)
							#result = PopenTmp.communicate()[0]
							#print result
							modelDir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'models')
							fileToPredict = outfile4+'_modified_selected_sorted_features_updated_stableClusterSeq_15.tsv'
							sampleNameTmp = '_'.join(os.path.basename(fileToPredict).split('_')[2:-10])
							if species not in ['human', 'mouse']:
								speciesType = 'others'
							else:
								speciesType = species
							if speciesType == 'others':
								print 'Notice: For %s, the predictive model is trained on human and mouse data.'%(species)
							novelOutputDir = os.path.join(outputdir, '_'.join(os.path.splitext(file)[0].split('_')[2:])+'_novel_miRNAs')
							os.mkdir(novelOutputDir)
							# Refine the feature files 
							preprocess_featureFiles(fileToPredict, os.path.join(modelDir, 'total_features_namelist.txt'))
							model_predict(fileToPredict, os.path.join(modelDir, speciesType+'_svc_model.pkl'))
							novelmiRNALListFile = os.path.join(os.path.dirname(fileToPredict), sampleNameTmp+'_novel_miRNAs_miRge2.0.csv')
							featureFile = fileToPredict
							clusterFile = outfile4+'_modified_selected_sorted_cluster.txt'
							write_novel_report(novelmiRNALListFile, featureFile, clusterFile, rnafoldCmdTmp)
							os.system('cp %s/*.pdf %s 2>/dev/null'%(outputdirTemp, novelOutputDir))
							os.system('cp %s/*_novel_miRNAs_report.csv %s 2>/dev/null'%(outputdirTemp, novelOutputDir))
							time_8 = time.time()
							print 'Prediction of novel miRNAs Completed (%.2f sec)'%(time_8-time_7)
			if fileType == 'unmapped.csv':
				os.system('rm -r %s'%(outputdir2))

if __name__ == '__main__':
	main()
