```
usage: miRge3.0 [options]

miRge3.0 (Comprehensive analysis of small RNA sequencing Data)

optional arguments:
  -h, --help  show this help message and exit
  --version   show program's version number and exit

Options:
  -s,    --samples            list of one or more samples separated by comma or a file with list of samples separated by new line (accepts *.fastq, *.fastq.gz) 
  -db,   --mir-DB             the reference database of miRNA. Options: miRBase and miRGeneDB (Default: miRBase) 
  -lib,  --libraries-path     the path to miRge libraries 
  -on,   --organism-name      the organism name can be human, mouse, fruitfly, nematode, rat or zebrafish
  -ex,   --crThreshold        the threshold of the proportion of canonical reads for the miRNAs to retain. Range for ex (0 - 0.5), (Default: 0.1)
  -phr,  --phred64           phred64 format (Default: 33)
  -spk,  --spikeIn            switch to annotate spike-ins if spike-in bowtie index files are located at the path of bowtie's index files (Default: off)
  -ie,   --isoform-entropy    switch to calculate isomir entropy (default: off)
  -cpu,  --threads            the number of processors to use for trimming, qc, and alignment (Default: 1)
  -ai,   --AtoI               switch to calculate A to I editing (Default: off)
  -tcf   --tcf-out            switch to write trimmed and collapsed fasta file (Default: off)
  -gff   --gff-out            switch to output isomiR results in gff format (Default: off) 
  -trf   --tRNA-frag          switch to analyze tRNA fragments and halves (Default: off)
  -o     --outDir             the directory of the outputs (Default: current directory) 

Data pre-processing:
  -a,    --adapter            Sequence of a 3' adapter. The adapter and subsequent bases are trimmed
  -g,    --front              Sequence of a 5' adapter. The adapter and any preceding bases are trimmed
  -u,    --cut                Remove bases from each read. If LENGTH is positive, remove bases from the beginning. If LENGTH is negative, remove bases from the end
  -nxt,  --nextseq-trim       NextSeq-specific quality trimming (each read). Trims also dark cycles appearing as high-quality G bases
  -q,    --quality-cutoff     Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. If one value is given, only the 3' end is trimmed
                              If two comma-separated cutoffs are given, the 5' end is trimmed with the first cutoff, the 3' end with the second
  -l,    --length             Shorten reads to LENGTH. Positive values remove bases at the end while negative ones remove bases at the beginning. This and the following
                              modifications are applied after adapter trimming
  -NX,   --trim-n             Trim N's on ends of reads
  -m,    --minimum-length     Discard reads shorter than LEN. (Default: 16)
  

Predicting novel miRNAs:
  novel miRNA detection is confined to human and mouse only (with "-on" argument):
  -nmir, --novel_miRNA        include prediction of novel miRNAs
  -c,    --minReadCounts      the minimum read counts supporting novel miRNA detection (default: 2)
  -mloc, --maxMappingLoci     the maximum number of mapping loci for the retained reads for novel miRNA detection (default: 3)
  -sl,   --seedLength         the seed length when invoking Bowtie for novel miRNA detection (default: 25)
  -olc,  --overlapLength      the length of overlapped seqence when joining reads into longer sequences based on the coordinate 
                              on the genome for novel miRNA detection (default: 14)
  -clc,  --clusterLength      the maximum length of the clustered sequences for novel miRNA detection (default: 30)

Optional PATH arguments:
  -pbwt, --bowtie-path        the path to system's directory containing bowtie binary
  -psam, --samtools-path      the path to system's directory containing samtools binary
  -prf,  --RNAfold-path       the path to system's directory containing RNAfold binary
```
