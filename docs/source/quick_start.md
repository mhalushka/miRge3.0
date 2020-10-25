# User guide

## Parameters

To view command-line parameters type `miRge3.0 -h`:
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
  -phr,  --phred64            phred64 format (Default: 33)
  -spk,  --spikeIn            switch to annotate spike-ins if spike-in bowtie index files are located at the path of bowtie's index files (Default: off)
  -ie,   --isoform-entropy    switch to calculate isomir entropy (default: off)
  -cpu,  --threads            the number of processors to use for trimming, qc, and alignment (Default: 1)
  -ai,   --AtoI               switch to calculate A to I editing (Default: off)
  -tcf   --tcf-out            switch to write trimmed and collapsed fasta file (Default: off)
  -gff   --gff-out            switch to output isomiR results in gff format (Default: off)
  -bam   --bam-out            switch to output isomiR results in gff format (Default: off)
  -trf   --tRNA-frag          switch to analyze tRNA fragment and halves (Default: off)
  -o     --outDir             the directory of the outputs (Default: current directory)
  -shh   --quiet              enable quiet/silent mode, only show warnings and errors (Default: off)

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
  -umi,  --uniq-mol-ids       Removes PCR duplicates and trim UMI of length by specifying two comma-separated cutoffs as 5’ cutoff,3’ bp from both ends of the read. eg: 4,4 or 0,4
  -udd,  --umiDedup           Specifies argument to removes PCR duplicates (Default: False); if TRUE it will remove UMI and remove PCR duplicates otherwise it only remove UMI and keep the raw counts
  -umiq, --umiqiagen          Removes PCR duplicates of reads obtained from Qiagen platform (Default: Illumina; "-umi x,y " Required)

Predicting novel miRNAs:
  The predictive model for novel miRNA detection is trained on human and mouse!
  -nmir, --novel-miRNA        include prediction of novel miRNAs
  -minl, --minLength          the minimum length of the reatined reads for novel miRNA detection (default: 16)
  -maxl, --maxLength          the maximum length of the reatined reads for novel miRNA detection (default: 25)
  -c,    --minReadCounts      the minimum read counts supporting novel miRNA detection (default: 2)
  -mloc, --maxMappingLoci     the maximum number of mapping loci for the retained reads for novel miRNA detection (default: 3)
  -sl,   --seedLength         the seed length when invoking Bowtie for novel miRNA detection (default: 25)
  -olc,  --overlapLenCutoff   the length of overlapped seqence when joining reads into longer sequences based on the coordinate
                              on the genome for novel miRNA detection (default: 14)
  -clc,  --clusterLength      the maximum length of the clustered sequences for novel miRNA detection (default: 30)

Optional PATH arguments:
  -pbwt, --bowtie-path        the path to system's directory containing bowtie binary
  -psam, --samtools-path      the path to system's directory containing samtools binary
  -prf,  --RNAfold-path       the path to system's directory containing RNAfold binary
    
```

## File format options
Having the right file format is important before making miRge libraries. When dealing with new species which are not available in the set of miRge3.0 libraries, it is important to prioritize what is essential.  Novel miRNAs runs scipy cKDTree during library preparation and it consumes a lot of computational resources and time depending on the genome size (up to 10 hours). Making a general build without novel miRNA detection is more straight forward and faster to build libraries.

#### General format options ####

##### Example usage #####
Example command usage:
```
miRge-build -g genome.fasta -mmf nematode_mature_miRBase.fa -hmf hairpin_miR.fa -mtf mature_trna.fasta -ptf pre_trna.fasta -snorf snorna.fasta -rrf rrna.fasta -ncof ncrna_other.fasta -mrf mrna.fasta -agff nematode_miRBase.gff3 -db miRBase -on roundworm -cpu 10  -ngrs WBcel235_genome_repeats.GTF
```
Output command line:
```
bowtie version: 1.2.3

Library indexing in progress...

Building the kdTree of roundworm_genome_repeats.GTF....

Building the kdTree of roundworm_genome_repeats.GTFtakes: 1.4s
Transforming roundworm_genome.fa takes: 0.9s

miRge-build is complete in 108.2122 second(s)
```
Output directory structure: 
```
DB = '--mir-DB'; name of the database used (miRBase or miRGeneDB)

Organism
├── annotation.Libs
│   ├── organism_DB.gff3
│   ├── organism_genome_repeats.pckl (if `-ngrs` is opted)
│   ├── organism_miRNAs_in_repetitive_element_DB.csv (if `-ngrs` is opted)
│   └── organism_merges_DB.csv
├── fasta.Libs
│   ├── organism_genome.pckl (if `-ngrs` is opted) 
│   └── organism_merges_DB.fa
└── index.Libs
    ├── organism_genome*.ebwt
    ├── organism_hairpin_DB*.ebwt
    ├── organism_mirna_DB*.ebwt
    ├── organism_mature_trna*.ebwt
    ├── organism_pre_trna*.ebwt
    ├── organism_rrna*.ebwt
    ├── organism_snorna*.ebwt
    ├── organism_mrna*.ebwt
    ├── organism_ncrna_others*.ebwt
    ├── organism_mature_trna*.ebwt
    └── organism_spike-in*.ebwt (Optional)
```

##### Name of the database #####
miRge uses miRBase or miRGeneDB as the reference database.
 So, it is mandatory to use `-db` option to either `-db miRBase` or `-db miRGeneDB`. Reference miRNA database `-db` and annotation GFF `-agff` files can be found at [miRGeneDB](https://mirgenedb.org/) and [miRBase](http://www.mirbase.org/). 

##### Name of the organism #####
miRge-build creates and stores all the libraries under the folder which is named after the organism. It is recommended to use a simple name and avoid any special character (use "_" if the name needs to be seperated by a space). Example: ` -on human `; ` -on horse `; `-on golden_lemur`; ` -on my_database ` etc.
    

##### Fasta format #####
Parameters with `-g`, `-mmf`, `-hmf`, `-mtf`, `-ptf`, `-snorf`, `-mrf`, `-spnf` should be in FASTA format as shown below. `-spnf ` or --spike-in is optional if the user is interested in adding an additional database with spike-in reads. 

*FASTA Format:*

```
>Header or Identifier
NUCLEOTIDE SEQUENCE 
Ex:
  >hsa-let-7a-5p
  TGAGGTAGTAGGTTGTATAGTT
```

**NOTE:**
```
The Header ID of hairpin miRNA FASTA should match the mature miRNA FASTA file. This is required for accurate isomiR annotation. 
miRge-build fetches 2bp upstream to 5p and 6bp downstream to 3p mature miRNA from the hairpin miRNA based on the matching ID. 
Exception: If the mature miRNA name contains XXXX-5p, XXXX-3p, XXXX-[5|3]p*,  XXXX_5p or XXXX_3p where XXXX matches the hairpin miRNA ID. 
Also, if this is not possible, miRge will not throw any errors, however, and it will proceed with the user provided files.  
```

#### Novel miRNA options ####
Novel miRNA prediction requires the genome file (which is provided in the general format) and genome repeats file in GTF format, `-ngrs`. As mentioned previously, novel miRNA analysis consumes a lot of computational resources and time.


#### Custom annotation options ####
This is **optional**, that two files under the `annotation.Libs` subdirectory requires users input manually. 

##### \_merges\_ #####
This file structured as `organism_merges_database.csv` allows users to define a miRNA family for miRNAs with similar sequences. This method is described in detail in the original miRge manuscript (Baras et al Plos One, 2015).
Below is the guide to format the file, where `hsa-miR-376b-5p/376c-5p` is the name of the miRNA family seperated by `/` followed by the family members such as `hsa-miR-376b-5p` and `hsa-miR-376c-5p` all separated by `,`. The next such miRNA family should begin in a new line. Here, four such examples are shown below. 

```
hsa-miR-376b-5p/376c-5p,hsa-miR-376b-5p,hsa-miR-376c-5p
hsa-miR-518c-3p/518f-3p,hsa-miR-518c-3p,hsa-miR-518f-3p
hsa-miR-642a-3p/642b-3p,hsa-miR-642a-3p,hsa-miR-642b-3p
hsa-miR-3155a-3p/3155b,hsa-miR-3155a-3p,hsa-miR-3155b
hsa-miR-3689b-3p/3689c,hsa-miR-3689b-3p,hsa-miR-3689c
```

##### \_miRNAs\_in\_repetitive\_element\_ #####

This file structured as `organism_miRNAs_in_repetitive_element_database.csv` allows users to define miRNAs that overlap with repeat elements in the genome. This eliminates miRNA reads to be identified as novel miRNAs or identifying one as A-to-I editing, both of which might be misleading. 

Below is the guide to format the file, where miRNA names which overlaps with repeat elements are separated by `,`. The `gene_id` and `transcript_id` of a repeat element should follow the miRNA name. See the example below: 

```
hsa-miR-28-5p,gene_id "L2c"; transcript_id "L2c_dup8856";
hsa-miR-28-3p,gene_id "L2c"; transcript_id "L2c_dup8856";
hsa-miR-95-5p,gene_id "L2c"; transcript_id "L2c_dup382";
hsa-miR-95-3p,gene_id "L2b"; transcript_id "L2b_dup437";
hsa-miR-181c-5p,gene_id "MamRTE1"; transcript_id "MamRTE1_dup11";
```

#### Resources ####
* The genome repeats can be obtained from [UCSC](https://genome-euro.ucsc.edu/cgi-bin/hgTables)
* The database sequences for other small RNA can be obtained from [UCSC](https://genome-euro.ucsc.edu/cgi-bin/hgTables) or [Ensembl](http://uswest.ensembl.org/Homo_sapiens/Info/Index)
* [Bowtie-v1.2.3](https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.3) - please pick one based on your OS.

