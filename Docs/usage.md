```
[ahanuma2@jhu.edu@bc-login01 mirge]$ python mirge.py -h
usage: miRge3.0 [options]

miRge3.0 (Comprehensive analysis of miRNA sequencing Data)

optional arguments:
  -h, --help            show this help message and exit

Options:
  -s [sample <required> [sample <required> ...]]
                        two options: 1. A file where each row represents one sample name;  2. *.fastq *.fastq ... Or *.fastq.gz *.fastq.gz ...
  -d <string required>  the miRNA database (default: miRBase. miRGeneDB is optional)
  -pb <dir required>    the path to the system's bowtie binary
  -lib <dir required>   the path to the miRge libraries
  -sp <string required>
                        the species can be human, mouse, fruitfly, nematode, rat and zebrafish (novel miRNA detection is confined in human and mouse)
  -ex <float>           the threshold of the proportion of canonical reads for the miRNAs to determine whether keeping them or not when counting. Users can set it between 0 and 0.5. (default: 0.1)
  -phred64              phred64 format (default: 64)
  -spikeIn              switch to annotate spike-ins if the bowtie index files are loacted at the path of bowtie's index files (default: off)
  -tcf                  switch to write trimmed and collapsed fasta file (default: off)
  -di                   switch to calculate of isomirs entropy (default: off)
  -cpu <int>            the number of processors to use for trimming, qc, and alignment (default: 1)
  -ai                   switch to calculate of A to I editing (default: off)
  -gff                  switch to output results in gff format (default: off)
  -trf                  switch to analyze tRNA fragment (default: off)
  --version             show program's version number and exit
  -o <dir>              the directory of the outputs (default: current directory)

CutAdapt commands:
  -a ADAPTER, --adapter ADAPTER
                        Sequence of an adapter ligated to the 3' end. The adapter and subsequent bases are trimmed. If a '$' character is appended ('anchoring'), the adapter is only
                                    found if it is a suffix of the read.
  -g ADAPTER, --front ADAPTER
                        Sequence of an adapter ligated to the 5' end. The adapter and any preceding bases are trimmed. Partial matches at the 5' end are allowed. If a '^' character is
                                    prepended ('anchoring'), the adapter is only found if it is a prefix of the read.
  -b ADAPTER, --anywhere ADAPTER
                        Sequence of an adapter that may be ligated to the 5' or 3' end. Both types of matches as described under -a und -g are allowed. If the first base of the read is
                                    part of the match, the behavior is as with -g, otherwise as with -a. This option is mostly for rescuing failed library preparations - do not use if you know
                                    which end your adapter was ligated to!
  -u LENGTH, --cut LENGTH
                        Remove bases from each read. If LENGTH is positive, remove bases from the beginning. If LENGTH is negative, remove bases from the end.
                                    Can be used twice if LENGTHs have different signs. This is applied *before* adapter trimming.
  --nextseq-trim 3'CUTOFF
                        NextSeq-specific quality trimming (each read). Trims also dark cycles appearing as high-quality G bases.
  -q [5'CUTOFF,]3'CUTOFF, --quality-cutoff [5'CUTOFF,]3'CUTOFF
                        Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. Applied to both reads if data is paired. If one
                                value is given, only the 3' end is trimmed. If two comma-separated cutoffs are given, the 5' end is trimmed with the first cutoff,
                                the 3' end with the second.
  --length LENGTH, -l LENGTH
                        Shorten reads to LENGTH. Positive values remove bases at the end while negative ones remove bases at the beginning. This and the following
                                modifications are applied after adapter trimming.
  --trim-n              Trim N's on ends of reads.
  -m LEN[:LEN2], --minimum-length LEN[:LEN2]
                        Discard reads shorter than LEN. Default: 0

Predicting novel miRNAs:
  -ps <dir required>    the path to the system's samtools binary
  -pr <dir required>    the path to the system's rnafold binary
  -ws <file>            the file containing the overall samples to analysis for novel miRNA prediction. No header, just a list of *.fastq file names in a column.
                        Names of files can be to your choosing (e.g. filestochecknovel.txt)
  -minl <int>           the minimum length of the reatined reads for novel miRNA detection (default: 16)
  -maxl <int>           the maximum length of the reatined reads for novel miRNA detection (default: 25)
  -cc <int>             the maximum read count of the reatined reads for novel miRNA detection (default: 2)
  -ml <int>             the maximum number of mapping loci for the retained reads for novel miRNA detection (default: 3)
  -sl <int>             the seed length when invoking Bowtie for novel miRNA detection (default: 25)
  -olc <int>            the length of overlapped seqence when joining reads into longer sequences based on the coordinate on the genome for novel miRNA detection (default: 14)
  -clc <int>            the maximum length of the clustered sequences for novel miRNA detection (default: 30)
```
