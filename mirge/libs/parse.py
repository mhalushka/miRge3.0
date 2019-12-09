import os
import sys
import argparse
import subprocess

def parseArg():
	parser = argparse.ArgumentParser()
	parser.parse_args()
	print("I am called")

"""
def parseArg():
	#dir_tmp = os.path.expanduser('~')
	dir_tmp = os.getcwd()
	version = '2.0'
	usageTmp = '\r{}\n\
##                                                                              ##\n\
##      miRge2.0 (Comprehensive analysis of miRNA sequencing Data)              ##\n\
##                                                                              ##\n\
##      last change: 06/26/2018                                                 ##\n\
##                                                                              ##\n\
##                                                                              ##\n\
##################################################################################\n\
\n'.format('##################################################################################'.ljust(len('usage:')))
	#create top-level parser
	usage = usageTmp+'Usage: %(prog)s <command> [<args>]\n\
The two functions of miRge2.0are:\n\
   annotate      Annotate the reads from miRNA sequencing data\n\
                 Type "miRge2.0 annotate -h" to show the help message of this funtion\n\
   predict       Detect novel miRNAs from miRNA sequencing data\n\
                 Type "miRge2.0 predict -h" to show the help message of this funtion\n'
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=usage, usage=argparse.SUPPRESS)
	subparsers = parser.add_subparsers(help='sub-command help')
	#create the parser for the 'annotate' command
	usage1 = usageTmp+'Usage: miRge2.0 annotate [-h] [<args>]\n\nExample:\nmiRge2.0 annotate -s seq_file.fastq -d miRBase -pb /home/yin/tools/bowtie-1.1.1 -lib /home/yin/miRge.Libs -sp human -ad illumina -ai -gff -trf -cpu 4\n'
	#print usage1
	parser_1 = subparsers.add_parser('annotate', usage = usage1)
	parser_1.add_argument('-s', nargs='*', required=True, dest='sampleList', metavar='sample <required>', help='two options: 1. A file where each row represents one sample name;  2. *.fastq *.fastq ... Or *.fastq.gz *.fastq.gz ...')
	parser_1.add_argument('-o', default=dir_tmp, dest='output_dir', metavar='<dir>', help='the directory of the outputs (default: current directory)')
	parser_1.add_argument('-d', default='miRBase',dest='miRNA_database', metavar='<string required>', help="the miRNA database (default: miRBase. miRGeneDB is optional)")
	parser_1.add_argument('-pb', required=True, dest='bowtieBinary', metavar='<dir required>', help="the path to the system's bowtie binary")
	parser_1.add_argument('-lib', required=True, dest= 'libraryPath', metavar='<dir required>', help="the path to the miRge libraries" )
	parser_1.add_argument('-sp', required=True, dest='species', metavar='<string required>', help="the species can be human, mouse, fruitfly, nematode, rat and zebrafish (novel miRNA detection is confined in human and mouse)")
	#parser_1.add_argument('-pi', required=True, dest='indexPath', metavar='<dir required>', help="The path to bowtie's index files, (index files of miRNA, hairpin, mrna, snoRNA and trna are compulsory, spike-ins is necessory when turning it on and index files of genome is necessory when launching novel miRNA detection)")
	#parser_1.add_argument('-mf', required=True, dest='mergeLibFile', metavar='<file required>', help='The file that contains the merged miRNAs for a specific sepeices')
	#parser_1.add_argument('-mif', required=True, dest='miRNA_fa', metavar='<file required>', help='the fasta file of miRNA including miRNA SNPs')
	#parser_1.add_argument('-sp', required=True, dest='species', metavar='<string required>', help="The species where the sequencing data are from. (novel miRNA detection is confined in human)")
	parser_1.add_argument('-ex', default='0.1', dest ='canoRatio', metavar='<float>', help='the threshold of the proportion of canonical reads for the miRNAs to determine whether keeping them or not when counting. Users can set it between 0 and 0.5. (default: 0.1)')
	parser_1.add_argument('-ad', default='none', dest='adapter', metavar='<string>', help='the adapter need to be removed which could be illumina, ion or a defined sequence (default: none)')
	parser_1.add_argument('-phred64', action = 'store_true', help='phred64 format (default: 64)')
	parser_1.add_argument('-spikeIn', dest='spikeIn', action = 'store_true', help="switch to annotate spike-ins if the bowtie index files are loacted at the path of bowtie's index files (default: off)")
	parser_1.add_argument('-tcf', dest='trimmed_collapsed_fa', action = 'store_true', help='switch to write trimmed and collapsed fasta file (default: off)')
	parser_1.add_argument('-di', dest='diff_isomirs', action = 'store_true', help='switch to calculate of isomirs entropy (default: off)')
	parser_1.add_argument('-cpu', dest='cpu', metavar='<int>', default='1', help='the number of processors to use for trimming, qc, and alignment (default: 1)')
	parser_1.add_argument('-ai', dest='a_to_i', action = 'store_true', help='switch to calculate of A to I editing (default: off)')
	#parser_1.add_argument('-mre', default='none', dest='miRNAs_in_repetitive_element', metavar='<file>', help='The file that contains the miRNA name located at repetetive element regions (compulsory when -ai is turned on). There will be no repetive element infromation if omitted')
	parser_1.add_argument('-gff', dest='gff_output', action = 'store_true', help='switch to output results in gff format (default: off)')
	parser_1.add_argument('-trf', dest='trf_output', action = 'store_true', help='switch to analyze tRNA fragment (default: off)')
	parser_1.add_argument('--version', action='version', version='%s'%(version))
	#create the parser for the 'predict' command
	usage2 = usageTmp+'Usage: miRge2.0 predict [-h] [<args>]\n\nExample:\nmiRge2.0 predict -s seq_file.fastq -d miRBase -pb /home/yin/tools/bowtie-1.1.1 -lib /home/yin/miRge.Libs -ps /usr/local/bin -pr /usr/local/bin -sp human -ad illumina -ai -gff -trf -cpu 4\n'
	parser_2 = subparsers.add_parser('predict', usage =usage2)
	parser_2.add_argument('-s', nargs='*', required=True, dest='sampleList', metavar='sample <required>', help='two options: 1. A file where each row represents one sample name;  2. *.fastq *.fastq ... Or *.fastq.gz *.fastq.gz ...')
	parser_2.add_argument('-o', default=dir_tmp, dest='output_dir', metavar='<dir>', help='the directory of the outputs (default: current directory)')
	parser_2.add_argument('-d', default='miRBase',dest='miRNA_database', metavar='<string required>', help="the miRNA database (default: miRBase. miRGeneDB is optional)")
	parser_2.add_argument('-pb', required=True, dest='bowtieBinary', metavar='<dir required>', help="the path to the system's bowtie binary")
	parser_2.add_argument('-lib', required=True, dest= 'libraryPath', metavar='<dir required>', help="the path to the miRge libraries" )
	parser_2.add_argument('-sp', required=True, dest='species', metavar='<string required>', help="the species can be human, mouse, fruitfly, nematode, rat and zebrafish (novel miRNA detection is confined in human and mouse)")
	#parser_2.add_argument('-pi', required=True, dest='indexPath', metavar='<dir required>', help="The path to bowtie's index files, (index files of miRNA, hairpin, mrna, snoRNA and trna are compulsory, and index files of genome is necessory when launching novel miRNA detection)")
	#parser_2.add_argument('-mf', required=True, dest='mergeLibFile', metavar='<file required>', help='The file that contains the merged miRNAs for a specific sepeices')
	#parser_2.add_argument('-mif', required=True, dest='miRNA_fa', metavar='<file required>', help='the fasta file of miRNA including miRNA SNPs')
	#parser_2.add_argument('-sp', required=True, dest='species', metavar='<string required>', help="The species where the sequencing data are from. (novel miRNA detection is confined in human)")
	#parser_2.add_argument('-gr', required=True, dest='genome_repeats', metavar='<file required>', help="The path of human genome repeats file (GRCH38_genome_repeats_sorted.pckl")
	#parser_2.add_argument('-gf', required=True, dest='genome_fa', metavar='<file required>', help="The path of human genome fasta file (human_genome.pckl)")
	#parser_2.add_argument('-maf', required=True, dest='mature_miRNA_fa', metavar='<file required>', help='the fasta file of human mature miRNA fasta file (hsa_mature.fa)')
	#parser_2.add_argument('-mic', required=True, dest='miRNA_coordinate', metavar='<file required>', help='the directory of human genome miRNA coordinate file (hsa_miRNA.gff3)')
	#parser_2.add_argument('-abbr', required=True, dest='abbrName', metavar='<string required>', help="The abbreviation name of species in miRBase")
	parser_2.add_argument('-ps', required=True, dest='samtoolsBinary', metavar='<dir required>', help="the path to the system's samtools binary")
	parser_2.add_argument('-pr', required=True, dest='rnafoldBinary', metavar='<dir required>', help="the path to the system's rnafold binary")
	parser_2.add_argument('-ex', default='0.1', dest ='canoRatio', metavar='<float>', help='the threshold of the proportion of canonical reads for the miRNAs to determine whether keeping them or not when counting. Users can set it between 0 and 0.5. (default: 0.1)')
	parser_2.add_argument('-ad', default='none', dest='adapter', metavar='<string>', help='the adapter need to be removed which could be illumina, ion or a defined sequence (default: none)')
	parser_2.add_argument('-phred64', action = 'store_true', help='phred64 format(default: 64)')
	parser_2.add_argument('-spikeIn', dest='spikeIn', action = 'store_true', help="switch to annotate spike-ins if the bowtie index files are loacted at the path of bowtie's index files (default: off)")
	parser_2.add_argument('-tcf', dest='trimmed_collapsed_fa', action = 'store_true', help='switch to write trimmed and collapsed fasta file (default: off)')
	parser_2.add_argument('-di', dest='diff_isomirs', action = 'store_true', help='switch to calculate of isomirs entropy (default: off)')
	parser_2.add_argument('-cpu', dest='cpu', metavar='<int>', default='1', help='the number of processors to use for trimming, qc, and alignment (default: 1)')
	parser_2.add_argument('-ai', dest='a_to_i', action = 'store_true', help='switch to calculate of A to I editing (default: off)')
	#parser_2.add_argument('-mre', default='none', dest='miRNAs_in_repetitive_element', metavar='<file>', help='The file that contains the miRNA name located at repetetive element regions (compulsory when -ai is turned on)')
	parser_2.add_argument('-gff', dest='gff_output', action = 'store_true', help='switch to output results in gff format (default: off)')
	parser_2.add_argument('-trf', dest='trf_output', action = 'store_true', help='switch to analyze tRNA fragment (default: off)')
	parser_2.add_argument('-ws',  default='none', dest ='wideSampleListFile', metavar='<file>', help='the file containing the overall samples to analysis for novel miRNA prediction. No header, just a list of *.fastq file names in a column. Names of files can be to your choosing (e.g. filestochecknovel.txt)')
	parser_2.add_argument('-minl', default='16', dest ='minLength', metavar='<int>', help='the minimum length of the reatined reads for novel miRNA detection (default: 16)')
	parser_2.add_argument('-maxl', default='25', dest='maxLength', metavar='<int>', help='the maximum length of the reatined reads for novel miRNA detection (default: 25)')
	parser_2.add_argument('-cc', default='2', dest='countCutoff', metavar='<int>', help='the maximum read count of the reatined reads for novel miRNA detection (default: 2)')
	parser_2.add_argument('-ml', default='3', dest='mapping_loc', metavar='<int>', help='the maximum number of mapping loci for the retained reads for novel miRNA detection (default: 3)')
	parser_2.add_argument('-sl', default='25', dest='seedLength', metavar='<int>', help='the seed length when invoking Bowtie for novel miRNA detection (default: 25)')
	parser_2.add_argument('-olc', default='14', dest='overlapLenCutoff', metavar='<int>', help='the length of overlapped seqence when joining reads into longer sequences based on the coordinate on the genome for novel miRNA detection (default: 14)')
	parser_2.add_argument('-clc', default='30', dest='clusterSeqLenCutoff', metavar='<int>', help='the maximum length of the clustered sequences for novel miRNA detection (default: 30)')
	parser_2.add_argument('--version', action='version', version='%s'%(version))
	return parser.parse_args()

"""
