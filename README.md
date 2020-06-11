(`Project_120919`)

miRge3.0 
======== 

An update to Python program to perform comprehensive analysis of miRNA sequencing Data, including miRNA annotation, A-to-I analysis, novel miRNA detection, isomiR analysis and tRF detection etc.


Documentation
-------------

* [Installation](#installation)
  * [Download libraries](#download-libraries)
  * [Install miRge3.0](#download-libraries)
  * [Troubleshooting installation](#troubleshooting-installation)
* [How to use it](#how-to-use-it)
* [Changelog](#changelog)
* [Citation](#citation)

Installation
------------

### Download libraries

miRge3.0 relies on a huge number of libraries like: <br />
1) Bowtie indexes of genome, hairping, mature miRNAs in miRBase, mature miRNAs in miRGeneDB, mRNA, rRNA, snoRNA, mature tRNA, primary tRNA, other ncRNA and spike-in sequences (optional) <br />
2) Sequences of genome, mature miRNAs (including SNP information) in miRBase and miRGeneDB <br />
3) Corrdinates of repetitive elements and mature miRNAs in the genome and miRNA merge information in miRBase and MirGeneDB <br />

Libraries of six species including ___human___, ___mouse___, ___rat___, ___zebrafish___, ___nematode___ and ___fruitfly___ can be downloaded separately at [`One Drive`](https://livejohnshopkins-my.sharepoint.com/personal/mhalush1_jh_edu/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fmhalush1%5Fjh%5Fedu%2FDocuments%2FmiRge3%2E0)
After unpacking the downloaded ***.tar.gz files to the new directory of miRge.Libs, the folder named by species contains three folders: index.Libs (libraries Part 1), fasta.Libs (libraries Part 2) and annotation.Libs (libraries Part 3). The absolute directory of miRge.Libs is used as the value of parameter ‘-lib’ in miRge3.0. <br /><br />
If the users want to build the libraries for other species, they can use scripts miRge_bowtie_build.py and miRge_pckls_build.py which can be downloaded from https://github.com/mhalushka/miRge_build, wherein miRge_bowtie_build.py is used to build bowtie index files and miRge_pckls_build.py is used to transform ***_genome_repeats.GTF and ***_genome.fa in oder to accelerate the speed of reading larg files into memory. 

### Install miRge3.0

NEEDS UPDATING

miRge3.0 is implemented as a Python program running on a Linux/Unix platform that requires pre-installation of Bowtie (v1.1.1 or v1.1.2; http://bowtie-bio.sourceforge.net/index.shtml), SAMtools (v1.5; http://samtools.sourceforge.net/) and RNAfold (v2.3.5; http://www.tbi.univie.ac.at/RNA). <br />
It was built with Python (v2.7.*) programming language and Python-related libraries, including cutadapt(v1.11 to v1.16), biopython(>= v1.68), numpy(>= v1.11.3), scipy(>= v0.17.0), matplotlib(>= v2.1.1), pandas(>= v0.21.0), sklearn(>= v0.18.1), reportlab(>= v3.3.0) and forgi(v0.20). <br />

PLEASE NOTE: miRge3.0 is currently incompatible with cutadapt v1.18.  Using v1.18 will give a "TypeError: __call__() takes exactly 3 arguments (2 given)" error.<br />

The source code is hosted at: https://github.com/mhalushka/miRge.<br />

miRge3.0 also can be installed from the source code by pip (THIS VERSION IS CURRENTLY UNDER REVISION. USE BIOCONDA FOR AN ACCURATE INSTALL UNTIL THIS MESSAGE DISAPPEARS):<br />
1) Download miRge3.0 source code from https://github.com/mhalushka/miRge and unzip the zipped file folder.<br />
2) If the package of wheel is not installed, run `pip install wheel` to install it.<br />
3) Change the directory to miRge3.0's directory and run `python setup.py bdist_wheel` to build a wheel file for the subsequent installation via pip.<br />
4) Run `pip install ./dist/mirge-0.1.40-py3-none-any.whl` to install miRge3.0.<br />


### Troubleshooting installation<br />
miRge3.0 was tested on the specific version of required softwares and python packages. Please make sure the version is correct.<br />
1) If Bowtie, SAMtools or RNAfold have been already installed in the system, please run `which bowtie`, `which samtools` or `which RNAfold` to find their installation paths. If the versions are incorrect, please install them with right version. 
2) Running `pip freeze` to check th version of current python packages. If some python packages can't work, please mannually install them by running `pip install package==*.**`.<br />
3) If the required python pacakages of the specific version can't be installed by pip or imported by python, make sure the installed python is complied by 4-byte Unicode so that pip can install UCS4 wheels (supporting cp27mu not cp27m). Type python and enter following commands `import sys` `print  sys.maxunicode`. If output is 1114111 then it is UCS4 otherwise if output is 65535 then it is UCS2.<br />
   If it is UCS2, please re-compile already installed python with 4-bype Unicode from the source code by running: a) `./configure --enable-unicode=ucs4 --prefix=***` b) `make` c) `make install` <br />
   

How to use it
-------------

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


Changelog
---------

To be updated


Links
-----

* [Documentation](<https://blank.readthedocs.io/>)
* [Source code](<https://github.com/mhalushka/mirge/>)
* [Report an issue](<https://github.com/mhalushka/mirge/issues>)
* [Project page on PyPI](<https://pypi.python.org/pypi/mirge/>)

Citation
--------

To be updated
