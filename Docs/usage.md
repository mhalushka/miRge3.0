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

Finding adapters:
 [REMOVE] Parameters -a, -g, -b specify adapters to be removed from
              each read (or from the first read in a pair if data is paired).
              If specified multiple times, only the best matching adapter is
              trimmed (but see the --times option). When the special notation
              'file:FILE' is used, adapter sequences are read from the given
              FASTA file.

  -a ADAPTER, --adapter ADAPTER
                        Sequence of an adapter ligated to the 3' end (paired data: of the
                                    first read). The adapter and subsequent bases are trimmed. If a
                                    '$' character is appended ('anchoring'), the adapter is only
                                    found if it is a suffix of the read.
  -g ADAPTER, --front ADAPTER
                        Sequence of an adapter ligated to the 5' end (paired data: of the
                                    first read). The adapter and any preceding bases are trimmed.
                                    Partial matches at the 5' end are allowed. If a '^' character is
                                    prepended ('anchoring'), the adapter is only found if it is a
                                    prefix of the read.
  -b ADAPTER, --anywhere ADAPTER
                        Sequence of an adapter that may be ligated to the 5' or 3' end
                                    (paired data: of the first read). Both types of matches as
                                    described under -a und -g are allowed. If the first base of the
                                    read is part of the match, the behavior is as with -g, otherwise
                                    as with -a. This option is mostly for rescuing failed library
                                    preparations - do not use if you know which end your adapter was
                                    ligated to!
  [REMOVE] -e RATE, --error-rate RATE
                        Maximum allowed error rate as value between 0 and 1 (no. of
                                    errors divided by length of matching region). Default: 0.1 (=10%)
  [REMOVE] --no-indels           Allow only mismatches in alignments.
                                    Default: allow both mismatches and indels
  [REMOVE]-n COUNT, --times COUNT
                        Remove up to COUNT adapters from each read. Default: 1
  [REMOVE]-O MINLENGTH, --overlap MINLENGTH
                        Require MINLENGTH overlap between read and adapter for an adapter
                                    to be found. Default: 3
  [REMOVE] --match-read-wildcards
                        Interpret IUPAC wildcards in reads. Default: False
  [REMOVE] -N, --no-match-adapter-wildcards
                        Do not interpret IUPAC wildcards in adapters.
  [REMOVE]--action {trim,mask,lowercase,none}
                        What to do with found adapters.
                                    mask: replace with 'N' characters;
                                    lowercase: convert to lowercase;
                                    none: leave unchanged (useful with
                                    --discard-untrimmed). Default: trim

Additional read modifications:
  -u LENGTH, --cut LENGTH
                        Remove bases from each read.
                                    If LENGTH is positive, remove bases from the beginning.
                                    If LENGTH is negative, remove bases from the end.
                                    Can be used twice if LENGTHs have different signs.
                                    This is applied *before* adapter trimming.
  --nextseq-trim 3'CUTOFF
                        NextSeq-specific quality trimming (each read). Trims also dark
                                    cycles appearing as high-quality G bases.
  -q [5'CUTOFF,]3'CUTOFF, --quality-cutoff [5'CUTOFF,]3'CUTOFF
                        Trim low-quality bases from 5' and/or 3' ends of each read before
                                    adapter removal. Applied to both reads if data is paired. If one
                                    value is given, only the 3' end is trimmed. If two
                                    comma-separated cutoffs are given, the 5' end is trimmed with
                                    the first cutoff, the 3' end with the second.
  [REMOVE]--quality-base N      Assume that quality values in FASTQ are encoded as ascii(quality
                                    + N). This needs to be set to 64 for some old Illumina
                                    FASTQ files. Default: 33
  --length LENGTH, -l LENGTH
                        Shorten reads to LENGTH. Positive values remove bases at the end
                                    while negative ones remove bases at the beginning. This and the
                                    following modifications are applied after adapter trimming.
  --trim-n              Trim N's on ends of reads.
  [REMOVE]--length-tag TAG      Search for TAG followed by a decimal number in the description
                                    field of the read. Replace the decimal number with the correct
                                    length of the trimmed read. For example, use --length-tag 'length='
                                    to correct fields like 'length=123'.
  [REMOVE]--strip-suffix STRIP_SUFFIX
                        Remove this suffix from read names if present. Can be given multiple times.
  [REMOVE]-x PREFIX, --prefix PREFIX
                        Add this prefix to read names. Use {name} to insert the name of the matching
                                    adapter.
  [REMOVE]-y SUFFIX, --suffix SUFFIX
                        Add this suffix to read names; can also include {name}
  [REMOVE]--zero-cap, -z        Change negative quality values to zero.

Filtering of processed reads:
  Filters are applied after above read modifications.
              "Paired-end reads are always discarded pairwise (see also
              --pair-filter).

  -m LEN[:LEN2], --minimum-length LEN[:LEN2]
                        Discard reads shorter than LEN. Default: 0
  [REMOVE]-M LEN[:LEN2], --maximum-length LEN[:LEN2]
                        Discard reads longer than LEN. Default: no limit
  [REMOVE]--max-n COUNT         Discard reads with more than COUNT 'N' bases. If COUNT is a number
                                     between 0 and 1, it is interpreted as a fraction of the read length.
  [REMOVE]--discard-trimmed, --discard
                        Discard reads that contain an adapter. Use also -O to avoid
                                    discarding too many randomly matching reads.
  [REMOVE]--discard-untrimmed, --trimmed-only
                        Discard reads that do not contain an adapter.
  [REMOVE]--discard-casava      Discard reads that did not pass CASAVA filtering (header has :Y:).

Output:
  [REMOVE]--quiet               Print only error messages.
  [REMOVE]--report {full,minimal}
                        Which type of report to print: 'full' or 'minimal'. Default: full
  [REMOVE]--fasta               Output FASTA to standard output even on FASTQ input.
  [REMOVE]-Z                    Use compression level 1 for gzipped output files (faster, but uses more space)
  [REMOVE]--info-file FILE      Write information about each read and its adapter matches into FILE.
                                    See the documentation for the file format.
  [REMOVE]-r FILE, --rest-file FILE
                        When the adapter matches in the middle of a read, write the
                                    rest (after the adapter) to FILE.
  [REMOVE]--wildcard-file FILE  When the adapter has N wildcard bases, write adapter bases
                                    matching wildcard positions to FILE. (Inaccurate with indels.)
  [REMOVE]--too-short-output FILE
                        Write reads that are too short (according to length specified by
                                -m) to FILE. Default: discard reads
  [REMOVE]--too-long-output FILE
                        Write reads that are too long (according to length specified by
                                -M) to FILE. Default: discard reads
  [REMOVE]--untrimmed-output FILE
                        Write reads that do not contain any adapter to FILE. Default:
                                    output to same file as trimmed reads

Predicting novel miRNAs:
  -t , --table          Specify table name to use
  -txto , --txtout      Writes the output of the query to a file speficied. Format (-fmt) is a tab-delimited text file by default
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
