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
    group.add_argument('-s', nargs='*', required=True, dest='sampleList', metavar='sample <required>', help='two options: 1. A file where each row represents one sample name;  2. *.fastq *.fastq ... Or *.fastq.gz *.fastq.gz ...')
    group.add_argument('-d', default='miRBase',dest='miRNA_database', metavar='<string required>', help="the miRNA database (default: miRBase. miRGeneDB is optional)")
    group.add_argument('-pb', required=True, dest='bowtieBinary', metavar='<dir required>', help="the path to the system's bowtie binary")
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
    group.add_argument('--version', action='version', version='%s'%(version))
    group.add_argument('-o', default="dir_tmp", dest='output_dir', metavar='<dir>', help='the directory of the outputs (default: current directory)')

    group1 = parser.add_argument_group("Finding adapters",
        description="""Parameters -a, -g, -b specify adapters to be removed from 
            each read (or from the first read in a pair if data is paired). 
            If specified multiple times, only the best matching adapter is 
            trimmed (but see the --times option). When the special notation 
            'file:FILE' is used, adapter sequences are read from the given 
            FASTA file.""")
    group1.add_argument("-a", "--adapter", type=lambda x: ("back", x), action="append",
        default=[], metavar="ADAPTER", dest="adapters",
        help="""Sequence of an adapter ligated to the 3' end (paired data: of the 
            first read). The adapter and subsequent bases are trimmed. If a 
            '$' character is appended ('anchoring'), the adapter is only 
            found if it is a suffix of the read.""")
    group1.add_argument("-g", "--front", type=lambda x: ("front", x), action="append",
        default=[], metavar="ADAPTER", dest="adapters",
        help="""Sequence of an adapter ligated to the 5' end (paired data: of the 
            first read). The adapter and any preceding bases are trimmed. 
            Partial matches at the 5' end are allowed. If a '^' character is 
            prepended ('anchoring'), the adapter is only found if it is a 
            prefix of the read.""")
    group1.add_argument("-b", "--anywhere", type=lambda x: ("anywhere", x), action="append",
        default=[], metavar="ADAPTER", dest="adapters",
        help="""Sequence of an adapter that may be ligated to the 5' or 3' end 
            (paired data: of the first read). Both types of matches as 
            described under -a und -g are allowed. If the first base of the 
            read is part of the match, the behavior is as with -g, otherwise 
            as with -a. This option is mostly for rescuing failed library 
            preparations - do not use if you know which end your adapter was 
            ligated to!""")
    group1.add_argument("-e", "--error-rate", type=float, default=0.1, metavar="RATE",
        help="""Maximum allowed error rate as value between 0 and 1 (no. of 
            errors divided by length of matching region). Default: %(default)s (=10%%)""")
    group1.add_argument("--no-indels", action='store_false', dest='indels', default=True,
        help="""Allow only mismatches in alignments. 
            Default: allow both mismatches and indels""")
    group1.add_argument("-n", "--times", type=int, metavar="COUNT", default=1,
        help="""Remove up to COUNT adapters from each read. Default: %(default)s""")
    group1.add_argument("-O", "--overlap", type=int, metavar="MINLENGTH", default=3,
        help="""Require MINLENGTH overlap between read and adapter for an adapter 
            to be found. Default: %(default)s""")
    group1.add_argument("--match-read-wildcards", action="store_true", default=False,
        help="""Interpret IUPAC wildcards in reads. Default: %(default)s""")
    group1.add_argument("-N", "--no-match-adapter-wildcards", action="store_false",
        default=True, dest='match_adapter_wildcards',
        help="""Do not interpret IUPAC wildcards in adapters.""")
    group1.add_argument("--action", choices=('trim', 'mask', 'lowercase', 'none'), default='trim',
        help="""What to do with found adapters. 
            mask: replace with 'N' characters; 
            lowercase: convert to lowercase; 
            none: leave unchanged (useful with 
            --discard-untrimmed). Default: %(default)s""")

    group1 = parser.add_argument_group("Additional read modifications")
    group1.add_argument("-u", "--cut", action='append', default=[], type=int, metavar="LENGTH",
        help="""Remove bases from each read. 
            If LENGTH is positive, remove bases from the beginning. 
            If LENGTH is negative, remove bases from the end. 
            Can be used twice if LENGTHs have different signs. 
            This is applied *before* adapter trimming.""")
    group1.add_argument("--nextseq-trim", type=int, default=None, metavar="3'CUTOFF",
        help="""NextSeq-specific quality trimming (each read). Trims also dark 
            cycles appearing as high-quality G bases.""")
    group1.add_argument("-q", "--quality-cutoff", default=None, metavar="[5'CUTOFF,]3'CUTOFF",
        help="""Trim low-quality bases from 5' and/or 3' ends of each read before 
            adapter removal. Applied to both reads if data is paired. If one 
            value is given, only the 3' end is trimmed. If two 
            comma-separated cutoffs are given, the 5' end is trimmed with 
            the first cutoff, the 3' end with the second.""")
    group1.add_argument("--quality-base", type=int, default=33, metavar='N',
        help="""Assume that quality values in FASTQ are encoded as ascii(quality 
            + N). This needs to be set to 64 for some old Illumina 
            FASTQ files. Default: %(default)s""")
    group1.add_argument("--length", "-l", type=int, default=None, metavar="LENGTH",
            help="""Shorten reads to LENGTH. Positive values remove bases at the end 
            while negative ones remove bases at the beginning. This and the 
            following modifications are applied after adapter trimming.""")
    group1.add_argument("--trim-n", action='store_true', default=False,
        help="""Trim N's on ends of reads.""")
    group1.add_argument("--length-tag", metavar="TAG",
        help="""Search for TAG followed by a decimal number in the description 
            field of the read. Replace the decimal number with the correct 
            length of the trimmed read. For example, use --length-tag 'length=' 
            to correct fields like 'length=123'.""")
    group1.add_argument("--strip-suffix", action='append', default=[],
        help="""Remove this suffix from read names if present. Can be given multiple times.""")
    group1.add_argument("-x", "--prefix", default='',
        help="""Add this prefix to read names. Use {name} to insert the name of the matching 
            adapter.""")
    group1.add_argument("-y", "--suffix", default='',
        help="""Add this suffix to read names; can also include {name}""")
    group1.add_argument("--zero-cap", "-z", action='store_true', default=False,
        help="""Change negative quality values to zero.""")

    group1 = parser.add_argument_group("Filtering of processed reads",
        description="""Filters are applied after above read modifications. 
            "Paired-end reads are always discarded pairwise (see also 
            --pair-filter).""")
    group1.add_argument("-m", "--minimum-length", default=None, metavar="LEN[:LEN2]",
        help="""Discard reads shorter than LEN. Default: 0""")
    group1.add_argument("-M", "--maximum-length", default=None, metavar="LEN[:LEN2]",
        help="""Discard reads longer than LEN. Default: no limit""")
    group1.add_argument("--max-n", type=float, default=None, metavar="COUNT",
        help="""Discard reads with more than COUNT 'N' bases. If COUNT is a number 
             between 0 and 1, it is interpreted as a fraction of the read length.""")
    group1.add_argument("--discard-trimmed", "--discard", action='store_true', default=False,
        help="""Discard reads that contain an adapter. Use also -O to avoid 
            discarding too many randomly matching reads.""")
    group1.add_argument("--discard-untrimmed", "--trimmed-only", action='store_true', default=False,
        help="""Discard reads that do not contain an adapter.""")
    group1.add_argument("--discard-casava", action='store_true', default=False,
        help="""Discard reads that did not pass CASAVA filtering (header has :Y:).""")

    group1 = parser.add_argument_group("Output")
    group1.add_argument("--quiet", default=False, action='store_true',
        help="""Print only error messages.""")
    group1.add_argument("--report", choices=('full', 'minimal'), default=None,
        help="""Which type of report to print: 'full' or 'minimal'. Default: full""")
    group1.add_argument("--fasta", default=False, action='store_true',
        help="""Output FASTA to standard output even on FASTQ input.""")
    group1.add_argument("-Z", action="store_const", const=1, dest="compression_level",
        help="""Use compression level 1 for gzipped output files (faster, but uses more space)""")
    group1.add_argument("--info-file", metavar="FILE",
        help="""Write information about each read and its adapter matches into FILE. 
            See the documentation for the file format.""")
    group1.add_argument("-r", "--rest-file", metavar="FILE",
        help="""When the adapter matches in the middle of a read, write the 
            rest (after the adapter) to FILE.""")
    group1.add_argument("--wildcard-file", metavar="FILE",
        help="""When the adapter has N wildcard bases, write adapter bases 
            matching wildcard positions to FILE. (Inaccurate with indels.)""")
    group1.add_argument("--too-short-output", metavar="FILE",
        help="""Write reads that are too short (according to length specified by 
        -m) to FILE. Default: discard reads""")
    group1.add_argument("--too-long-output", metavar="FILE",
        help="""Write reads that are too long (according to length specified by
        -M) to FILE. Default: discard reads""")
    group1.add_argument("--untrimmed-output", default=None, metavar="FILE",
        help="""Write reads that do not contain any adapter to FILE. Default: 
            output to same file as trimmed reads""")
    

    group2 = parser.add_argument_group('Predicting novel miRNAs')
    group2.add_argument("-t", "--table", metavar='', help="Specify table name to use")
    group2.add_argument("-txto", "--txtout", metavar='', help="Writes the output of the query to a file speficied. Format (-fmt) is a tab-delimited text file by default")
    group2.add_argument('-ps', required=True, dest='samtoolsBinary', metavar='<dir required>', help="the path to the system's samtools binary")
    group2.add_argument('-pr', required=True, dest='rnafoldBinary', metavar='<dir required>', help="the path to the system's rnafold binary")
    group2.add_argument('-ws',  default='none', dest ='wideSampleListFile', metavar='<file>', help=
    """the file containing the overall samples to analysis for novel miRNA prediction. No header, just a list of *.fastq file names in a column. 
Names of files can be to your choosing (e.g. filestochecknovel.txt)""")
    group2.add_argument('-minl', default='16', dest ='minLength', metavar='<int>', help='the minimum length of the reatined reads for novel miRNA detection (default: 16)')
    group2.add_argument('-maxl', default='25', dest='maxLength', metavar='<int>', help='the maximum length of the reatined reads for novel miRNA detection (default: 25)')
    group2.add_argument('-cc', default='2', dest='countCutoff', metavar='<int>', help='the maximum read count of the reatined reads for novel miRNA detection (default: 2)')
    group2.add_argument('-ml', default='3', dest='mapping_loc', metavar='<int>', help='the maximum number of mapping loci for the retained reads for novel miRNA detection (default: 3)')
    group2.add_argument('-sl', default='25', dest='seedLength', metavar='<int>', help='the seed length when invoking Bowtie for novel miRNA detection (default: 25)')
    group2.add_argument('-olc', default='14', dest='overlapLenCutoff', metavar='<int>', help='the length of overlapped seqence when joining reads into longer sequences based on the coordinate on the genome for novel miRNA detection (default: 14)')
    group2.add_argument('-clc', default='30', dest='clusterSeqLenCutoff', metavar='<int>', help='the maximum length of the clustered sequences for novel miRNA detection (default: 30)')

    
    return parser.parse_args()

