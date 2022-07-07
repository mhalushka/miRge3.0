import os
import sys
import argparse
import subprocess

def parseArg():
    version = '3.0'
    parser = argparse.ArgumentParser(description='miRge3.0 (Comprehensive analysis of small RNA sequencing Data)',usage='miRge3.0 [options]',formatter_class=argparse.RawTextHelpFormatter,)
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    parser.add_argument('--version', action='version', version='%s'%(version))
    group = parser.add_argument_group("Options",description='''-s,    --samples            list of one or more samples separated by comma or a file with list of samples separated by new line (accepts *.fastq, *.fastq.gz) 
-db,   --mir-DB             the reference database of miRNA. Options: miRBase and miRGeneDB (Default: miRBase) 
-lib,  --libraries-path     the path to miRge libraries 
-on,   --organism-name      the organism name can be human, mouse, fruitfly, nematode, rat or zebrafish
-ex,   --crThreshold        the threshold of the proportion of canonical reads for the miRNAs to retain. Range for ex (0 - 0.5), (Default: 0.1)
-phr,  --phred64            phred64 format (Default: 33)
-spk,  --spikeIn            switch to annotate spike-ins if spike-in bowtie index files are located at the path of bowtie's index files (Default: off)
-ie,   --isoform-entropy    switch to calculate isomir entropy (default: off)
-cpu,  --threads            the number of processors to use for trimming, qc, and alignment (Default: 1)
-ai,   --AtoI               switch to calculate A to I editing (Default: off)
-tcf   --tcf-out            switch to write trimmed and collapsed fasta file (Default: off)
-gff   --gff-out            switch to output isomiR results in gff format (Default: off) 
-bam   --bam-out            switch to output results in bam format (Default: off) 
-trf   --tRNA-frag          switch to analyze tRNA fragment and halves (Default: off)
-o     --outDir             the directory of the outputs (Default: current directory) 
-dex   --diffex             perform differential expression with DESeq2 (Default: off)
-mdt   --metadata           the path to metadata file (Default: off, require '.csv' file format if -dex is opted)
-cms   --chunkmbs           chunk memory in megabytes per thread to use during bowtie alignment (Default: 256)
-spl   --save-pkl           save collapsed reads in binary format for later runs (Default: off)
-rr    --resume             resume from collapsed reads (Default: off)
-shh   --quiet              enable quiet/silent mode, only show warnings and errors (Default: off)
''')
    group.add_argument('-s','--samples', nargs='*', required=True, help=argparse.SUPPRESS)
    group.add_argument('-db', '--mir-DB', default='miRBase', required=True, help=argparse.SUPPRESS) 
    group.add_argument('-lib', '--libraries-path', required=True, help=argparse.SUPPRESS)
    group.add_argument('-on','--organism-name', required=True, help=argparse.SUPPRESS)
    group.add_argument('-ex', '--crThreshold', default='0.1', help=argparse.SUPPRESS)
    group.add_argument('-phr', '--phred64', type=int, default=33, help=argparse.SUPPRESS)
    group.add_argument('-spk', '--spikeIn', action='store_true', default=False, help=argparse.SUPPRESS)
    group.add_argument('-ie', '--isoform-entropy', action='store_true', default=False, help=argparse.SUPPRESS)
    group.add_argument('-cpu', '--threads', type=int, default='0', help=argparse.SUPPRESS)
    group.add_argument('-ai', '--AtoI', action='store_true', default=False, help=argparse.SUPPRESS)
    group.add_argument('-tcf', '--tcf-out', action='store_true', default=False, help=argparse.SUPPRESS)
    group.add_argument('-bam', '--bam-out', action='store_true', default=False, help=argparse.SUPPRESS)
    group.add_argument('-gff', '--gff-out', action='store_true', default=False, help=argparse.SUPPRESS)
    group.add_argument('-trf', '--tRNA-frag', action='store_true', default=False, help=argparse.SUPPRESS)
    group.add_argument('-o', '--outDir', help=argparse.SUPPRESS)
    group.add_argument('-dex', '--diffex', action='store_true', default=False, help=argparse.SUPPRESS)
    group.add_argument('-mdt', '--metadata', help=argparse.SUPPRESS)
    group.add_argument('-onam', '--outDirName', help=argparse.SUPPRESS)
    group.add_argument('-cms', '--chunkmbs', type=int, default=256, help=argparse.SUPPRESS)
    group.add_argument('-spl', '--save-pkl', action='store_true', default=False, help=argparse.SUPPRESS)
    group.add_argument('-rr', '--resume', action='store_true', default=False, help=argparse.SUPPRESS)
    group.add_argument('-shh',"--quiet", default=False, action='store_true', help=argparse.SUPPRESS)

    group1 = parser.add_argument_group("Data pre-processing", description='''-a,    --adapter            Sequence of a 3' adapter. The adapter and subsequent bases are trimmed
-g,    --front              Sequence of a 5' adapter. The adapter and any preceding bases are trimmed
-u,    --cut                Remove bases from each read. If LENGTH is positive, remove bases from the beginning. If LENGTH is negative, remove bases from the end
-nxt,  --nextseq-trim       NextSeq-specific quality trimming (each read). Trims also dark cycles appearing as high-quality G bases
-q,    --quality-cutoff     Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. If one value is given, only the 3' end is trimmed
                            If two comma-separated cutoffs are given, the 5' end is trimmed with the first cutoff, the 3' end with the second
-l,    --length             Shorten reads to LENGTH. Positive values remove bases at the end while negative ones remove bases at the beginning. This and the following
                            modifications are applied after adapter trimming
-NX,   --trim-n             Trim N's on ends of reads
-m,    --minimum-length     Discard reads shorter than LEN. (Default: 16)
-umi,  --uniq-mol-ids       Trim nucleotides of specific length at 5’ and 3’ ends of the read, after adapter trimming. eg: 4,4 or 0,4. (Use -udd to remove PCR duplicates)  
-udd,  --umiDedup           Specifies argument to removes PCR duplicates (Default: False); if TRUE it will remove UMI and remove PCR duplicates otherwise it only remove UMI and keep the raw counts (Require -umi option)
-qumi, --qiagenumi          Removes PCR duplicates of reads obtained from Qiagen platform (Default: Illumina; "-umi x,y " Required)


''')
    group1.add_argument("-a", "--adapter", type=lambda x: ("back", x), action="append",
        default=[], dest="adapters", help=argparse.SUPPRESS)
    group1.add_argument("-g", "--front", type=lambda x: ("front", x), action="append",
        default=[], dest="adapters", help=argparse.SUPPRESS)
    group1.add_argument("-u", "--cut", action='append', default=[], type=int, metavar="LENGTH", help=argparse.SUPPRESS)
    group1.add_argument("-nxt","--nextseq-trim", type=int, default=None, metavar="3'CUTOFF", help=argparse.SUPPRESS)
    group1.add_argument("-q", "--quality-cutoff", default="10", metavar="[5'CUTOFF,]3'CUTOFF", help=argparse.SUPPRESS)
    group1.add_argument("--length", "-l", type=int, default=3, metavar="LENGTH", help=argparse.SUPPRESS)
    group1.add_argument("-NX", "--trim-n", action='store_true', default=False, help=argparse.SUPPRESS)
    group1.add_argument("-m", "--minimum-length", default=16, help=argparse.SUPPRESS)
    group1.add_argument("-umi", "--uniq-mol-ids", default=None, help=argparse.SUPPRESS)
    group1.add_argument("-qumi", "--qiagenumi", action='store_true', default=False, help=argparse.SUPPRESS)
    group1.add_argument("-udd", "--umiDedup", action='store_true', default=False, help=argparse.SUPPRESS)
    #### we use none of the following cutadapt options but are required to pass default values for miRNA and cutadapt pipeline
    group1.add_argument("-op", "--output", metavar="FILE", help=argparse.SUPPRESS) #"Default: write to standard output"
    group1.add_argument("--compression-level", type=int, default=6, help=argparse.SUPPRESS)
    group1.add_argument("--overlap", type=int, metavar="MINLENGTH", default=3, help=argparse.SUPPRESS)
    group1.add_argument("--error-rate", type=float, default=0.12, metavar="RATE",help=argparse.SUPPRESS)
    group1.add_argument("--gc-content", type=float, default=50,  help=argparse.SUPPRESS)
    group1.add_argument("-cuv","--cutadaptVersion", help=argparse.SUPPRESS)
    group1.add_argument("-buv","--bowtieVersion", help=argparse.SUPPRESS)
    group1.add_argument("--action", choices=('trim', 'mask', 'lowercase', 'none'), default='trim', help=argparse.SUPPRESS)
    group1.add_argument("-n", "--times", type=int, metavar="COUNT", default=1, help=argparse.SUPPRESS)
    group1.add_argument("--match-read-wildcards", action="store_true", default=False,help=argparse.SUPPRESS)
    group1.add_argument("-N", "--no-match-adapter-wildcards", action="store_false", default=True, dest='match_adapter_wildcards',help=argparse.SUPPRESS)
    group1.add_argument("--fasta", default=False, action='store_true',help=argparse.SUPPRESS)
    group1.add_argument("--buffer-size", type=int, default=4000000, help=argparse.SUPPRESS)
    group1.add_argument("--no-indels", action='store_false', dest='indels', default=True, help=argparse.SUPPRESS)
    group1.add_argument("-M", "--maximum-length", default=None, type=int, metavar="LEN[:LEN2]", help=argparse.SUPPRESS)
    group1.add_argument("--numba-pll", default=None, action='store_false', help=argparse.SUPPRESS)
    group1.add_argument("--numba-cuda", default=None, action='store_false', help=argparse.SUPPRESS)
    
    group4 = parser.add_argument_group('miRNA Error Correction', description='''microRNA correction method for single base substitutions due to sequencing errors (Note: Refines reads at the expense of time)
-mEC,  --miREC              Enable miRNA error correction (miREC)
-kh,   --threshold          the value for frequency threshold τ (Default kh = 5)
-ks,   --kmer-start         kmer range start value (k_1, default 15) 
-ke,   --kmer-end           kmer range end value (k_end, default 20)
''')
    group4.add_argument('-mEC','--miREC', help=argparse.SUPPRESS, action='store_true')
    group4.add_argument('-kh','--threshold', help=argparse.SUPPRESS, default='5', metavar='', type=int)
    group4.add_argument('-ks','--kmer-start', help=argparse.SUPPRESS, default='15', metavar='', type=int)
    group4.add_argument('-ke','--kmer-end', help=argparse.SUPPRESS, default='20', metavar='', type=int)
    ## group - 2 ##

    group2 = parser.add_argument_group('Predicting novel miRNAs',description='''The predictive model for novel miRNA detection is trained on human and mouse!
-nmir, --novel-miRNA        include prediction of novel miRNAs
-minl, --minLength          the minimum length of the retained reads for novel miRNA detection (default: 16)
-maxl, --maxLength          the maximum length of the retained reads for novel miRNA detection (default: 25)
-c,    --minReadCounts      the minimum read counts supporting novel miRNA detection (default: 2)
-mloc, --maxMappingLoci     the maximum number of mapping loci for the retained reads for novel miRNA detection (default: 3)
-sl,   --seedLength         the seed length when invoking Bowtie for novel miRNA detection (default: 25)
-olc,  --overlapLenCutoff   the length of overlapped seqence when joining reads into longer sequences based on the coordinate 
                            on the genome for novel miRNA detection (default: 14)
-clc,  --clusterLength      the maximum length of the clustered sequences for novel miRNA detection (default: 30)
''')
    group2.add_argument('-nmir', '--novel-miRNA', help=argparse.SUPPRESS, action='store_true')
    group2.add_argument('-minl', '--minLength', type=int, default='16', metavar='', help=argparse.SUPPRESS)
    group2.add_argument('-maxl', '--maxLength', type=int, default='25', metavar='', help=argparse.SUPPRESS)
    group2.add_argument('-c', '--minReadCounts', default='2', metavar='', help=argparse.SUPPRESS)
    group2.add_argument('-mloc', '--maxMappingLoci', default='3', metavar='', help=argparse.SUPPRESS)
    group2.add_argument('-sl','--seedLength', default='25', metavar='', help=argparse.SUPPRESS)
    group2.add_argument('-olc', default='14', dest='overlapLenCutoff', metavar='', help=argparse.SUPPRESS)
    group2.add_argument('-clc', default='30', dest='clusterLength', metavar='', help=argparse.SUPPRESS)


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

    

#group.add_argument('-ad', default='none', dest='adapter', metavar='<string>', help='the adapter need to be removed which could be illumina, ion or a defined sequence (default: none)')
