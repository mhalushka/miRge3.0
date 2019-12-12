import os
import sys
import argparse
import subprocess

def parseArg():
    version = '3.0'
    parser = argparse.ArgumentParser(description='miRge3.0 (Comprehensive analysis of miRNA sequencing Data)',usage='miRge3.0 [options]',formatter_class=argparse.RawTextHelpFormatter,)
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    group = parser.add_argument_group("Options")
    group.add_argument('-s', nargs='*', required=True, dest='sampleList', metavar='samples', help='two options: 1. A file where each row represents one sample name;  2. *.fastq *.fastq ... Or *.fastq.gz *.fastq.gz ...')
    group.add_argument('-d', default='miRBase',dest='miRNA_database', metavar='<string required>', help="the miRNA database (default: miRBase. miRGeneDB is optional)")
    group.add_argument('-lib', required=True, dest= 'libraryPath', metavar='<dir required>', help="the path to the miRge libraries" )
    group.add_argument('-sp', required=True, dest='species', metavar='<string required>', help="the species can be human, mouse, fruitfly, nematode, rat and zebrafish (novel miRNA detection is confined in human and mouse)")
    group.add_argument('-ex', default='0.1', dest ='canoRatio', metavar='<float>', help='the threshold of the proportion of canonical reads for the miRNAs to determine whether keeping them or not when counting. Users can set it between 0 and 0.5. (default: 0.1)')
    #group.add_argument('-ad', default='none', dest='adapter', metavar='<string>', help='the adapter need to be removed which could be illumina, ion or a defined sequence (default: none)')
    group.add_argument('-phred64', action = 'store_true', help='phred64 format (default: 64)')
    group.add_argument('-spikeIn', dest='spikeIn', action = 'store_true', help="switch to annotate spike-ins if the bowtie index files are loacted at the path of bowtie's index files (default: off)")
    group.add_argument('-tcf', dest='trimmed_collapsed_fa', action = 'store_true', help='switch to write trimmed and collapsed fasta file (default: off)')
    group.add_argument('-di', dest='diff_isomirs', action = 'store_true', help='switch to calculate of isomirs entropy (default: off)')
    group.add_argument('-cpu', dest='cpu', metavar='<int>', default='1', help='the number of processors to use for trimming, qc, and alignment (default: 1)')
    group.add_argument('-ai', dest='a_to_i', action = 'store_true', help='switch to calculate of A to I editing (default: off)')
    group.add_argument('-gff', dest='gff_output', action = 'store_true', help='switch to output results in gff format (default: off)')
    group.add_argument('-trf', dest='trf_output', action = 'store_true', help='switch to analyze tRNA fragment (default: off)')
    group.add_argument('-o', default="dir_tmp", dest='output_dir', metavar='<dir>', help='the directory of the outputs (default: current directory)')
    group.add_argument('--version', action='version', version='%s'%(version))

    group1 = parser.add_argument_group("Data pre-processing")
    group1.add_argument("-a", "--adapter", type=lambda x: ("back", x), action="append",
        default=[], metavar="ADAPTER", dest="adapters", help="""Sequence of a 3' adapter. The adapter and subsequent bases are trimmed.""")
    group1.add_argument("-g", "--front", type=lambda x: ("front", x), action="append",
        default=[], metavar="ADAPTER", dest="adapters", help="""Sequence of a 5' adapter. The adapter and any preceding bases are trimmed.""")
    group1.add_argument("-u", "--cut", action='append', default=[], type=int, metavar="LENGTH",
        help="""Remove bases from each read. If LENGTH is positive, remove bases from the beginning. If LENGTH is negative, remove bases from the end.""")
    group1.add_argument("--nextseq-trim", type=int, default=None, metavar="3'CUTOFF",
        help="""NextSeq-specific quality trimming (each read). Trims also dark cycles appearing as high-quality G bases.""")
    group1.add_argument("-q", "--quality-cutoff", default=None, metavar="[5'CUTOFF,]3'CUTOFF",
        help="""Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. If one value is given, only the 3' end is trimmed.
	 If two comma-separated cutoffs are given, the 5' end is trimmed with the first cutoff, the 3' end with the second.""")
    group1.add_argument("--length", "-l", type=int, default=None, metavar="LENGTH",
            help="""Shorten reads to LENGTH. Positive values remove bases at the end while negative ones remove bases at the beginning. This and the following 
	modifications are applied after adapter trimming.""")
    group1.add_argument("--trim-n", action='store_true', default=False, help="""Trim N's on ends of reads.""")
    group1.add_argument("-m", "--minimum-length", default=None, metavar="LEN[:LEN2]", help="""Discard reads shorter than LEN. Default: 0""")
    group1.add_argument('-minl', default='16', dest ='minLength', metavar='<int>', help='the minimum length of the reatined reads for novel miRNA detection (default: 16)')
    ## group - 3 ##

    group2 = parser.add_argument_group('Predicting novel miRNAs',description='''-nmir, --novel_miRNA        include prediction of novel miRNAs
-c,    --minReadCounts      the minimum read counts supporting novel miRNA detection (default: 2)
-mloc, --maxMappingLoci     the maximum number of mapping loci for the retained reads for novel miRNA detection (default: 3)
-sl,   --seedLength         the seed length when invoking Bowtie for novel miRNA detection (default: 25)
-olc,  --overlapLength      the length of overlapped seqence when joining reads into longer sequences based on the coordinate 
                            on the genome for novel miRNA detection (default: 14)
-clc,  --clusterLength      the maximum length of the clustered sequences for novel miRNA detection (default: 30)
''')
    group2.add_argument('-nmir', '--novel_miRNA', help=argparse.SUPPRESS, action='store_true')
    group2.add_argument('-c', '--minReadCounts', default='2', metavar='', help=argparse.SUPPRESS)
    group2.add_argument('-mloc', '--maxMappingLoci', default='3', metavar='', help=argparse.SUPPRESS)
    group2.add_argument('-sl','--seedLength', default='25', metavar='', help=argparse.SUPPRESS)
    group2.add_argument('-olc', default='14', dest='overlapLenCutoff', metavar='', help=argparse.SUPPRESS)
    group2.add_argument('-clc', default='30', dest='clusterSeqLenCutoff', metavar='', help=argparse.SUPPRESS)


    ## group - 3 ##

    group3 = parser.add_argument_group('Optional PATH arguments',
	description='''-pbwt, --bowtie-path        the path to system's directory containing bowtie binary
-psam, --samtools-path      the path to system's directory containing samtools binary
-prf,  --RNAfold-path       the path to system's directory containing RNAfold binary
''')
    group3.add_argument('-pbwt', '--bowtie-path', metavar="", help=argparse.SUPPRESS)
    group3.add_argument('-psam', '--samtools-path', metavar="", help=argparse.SUPPRESS)
    group3.add_argument('-prf', '--RNAfold-path', metavar="", help=argparse.SUPPRESS)
    
    return parser.parse_args()

