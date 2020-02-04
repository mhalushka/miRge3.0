# miRge-Py3

## Main goals
- [ ] Easy to implement and use
- [ ] Robust (accurate) data output
- [ ] Fast and able to run multiple datasets at once
- [ ] Attractive and meaningful output

## GitHub Issues 
- [x] To integrate miRge --build-idx (OK - Should be implemented)
- [x] Discuss spike-in and how it matters for users! (OK - Should be implemented)
[Ref1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5100345/)
[Ref2](https://www.ncbi.nlm.nih.gov/pubmed/25870415?dopt=Abstract)
[Ref3](https://www.nature.com/articles/s41598-017-06174-3)
- [ ] Find data with spike-ins and visualize/interpret the output (Pending)
- [ ] miRge --build-lib (Issue: #1) (Pending)
- [ ] miRge --allele-freq (Issue: #3) (Pending)
- [x] A zero count and an absent record in miR.Counts.csv (Issue: #20) (Complete)
- [x] Option: -ad ion (Issue: #25) - I think its mostly due to the use of different cutadapt version. (Complete).  The 11 is defined in the code below. (Ans: Ion torrent adapters should be provided by user unless modified in miRge3 in near future)
- [x] Any other type of adapters possible? (Ans: Any adapter seq can be provided on 5' `-g` and 3' `-a` end, or mention Illumina for adapters at 5' or 3')
```
adapter = args.adapter
	if adapter == 'illumina':
		adapter = 'TGGAATTCTCGGGTGCCAAGGAACTCCAG'
	elif adapter == 'ion':
		adapter = '11'
``` 
**Adapter sequences**

[Accurate Adapter Information Is Crucial for Reproducibility and Reusability in Small RNA Seq Studies](https://www.mdpi.com/2311-553X/5/4/49).
[Supplementary](https://www.mdpi.com/2311-553X/5/4/49/s1)
```
1. RealSeq®-AC =  TGGAATTCTCGGGTGCCAAGG (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6120088/)
2. NEBNext® (NEB) =  AGATCGGAAGAGCACACGTCT (New England Biolabs) {SRR6464616,SRR6464623}
3. TruSeq® (Illumina) = TGGAATTCTCGGGTGCCAAGG {SRR6464662, SRR6464661}
4. NEXTFlex™ (Bioo Scientific/PerkinElmer) = TGGAATTCTCGGGTGCCAAGG {SRR6464672,SRR6464671}
5. QIAseq (Qiagen) = AACTGTAGGCACCATCAAT {}
6. SMARTer (Takara Bio) = AAAAAAAAAA {SRR6464721, SRR6464720}
(a ligation-free technique that uses poly(A)-tailing and reverse transcriptase (RT) template switching to detect small RNAs)
7. TailorMix = TGGAATTCTCGGGTGCCAAGG
8. CleanTag = TGGAATTCTCGGGTGCCAAGG
9. Lexogen - Small RNA-Seq Library Prep Kit = TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC
10. CATS small RNA-seq Kit = GATCGGAAGAGCACACGTCTG {SRR6464674, SRR6464673}
11. Incorporation of UMIs into TruSeq adapters (TrUMIseq adapters) 
(https://www.future-science.com/doi/10.2144/000114608?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub%3Dwww.ncbi.nlm.nih.gov&)
Ref: https://github.com/LieberInstitute/miRNA_Kit_Comparison
Datasets: SRR* for different lib kit and adapters etc: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?study=SRP199350
```

## Basic functions - command line options
- [ ] User friendly parameters 
- [ ] Implement [`Read the docs`](https://readthedocs.org/)
  - [ ] Mention file type to use for loading sample information 
  - [ ] `Additional details to add to docs here`
- [ ] `Additional features here`

## Primary functions - Annotation:
- [x] Incorporate latest version of cutadapt and its features. 
- [ ] Keep track to update features of miRge automatically with newer versions of cutadapt.
- [ ] No installation required for [Bowtie 1.2.3](https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.3/). Use: `\bin\bowtie`
- [x] ~~**miRge --paste**~~. Option to comibne colummns of several \*.csv files into one matrix using panda dataframe (The output miR_Counts.csv and miR_RPM.csv can handle that automatically)
- [ ] `Additional features here`

## Secondary functions - Predict:
- [ ] RNA-fold for novel miRNAs
- [ ] Correct miRge GFF output
- [ ] tRNA halve and fragment prediction
- [ ] SVM re-implementation 
- [ ] PDF and HTML output
- [ ] Incorporate SQLite db for quick query and visualization
- [ ] A function to handle Batch effect? 
- [ ] `Add features here`

## Additional features - Customization of software:
- [ ] Bowtie to avoid mismatches in the seed region - Interact with **Mihaela Pertea**. 
- [ ] Software to use .SRA, .fastq.gz and .fastq
- [ ] `Additional features here`

## Advanced features - Packaging and testing
- [ ] Integrate Cython (better performance) 
- [ ] CPU and GPU switch usage options 
- [ ] Interactive or visually appealing outputs 
- [ ] Docker implementation 
- [ ] Integrates with other high througput Genomics tools, especially having output work in a Bioconductor environment.
- [ ] GUI - [Beeware](https://beeware.org/project/using/desktop-app/) or [Electron](https://electronjs.org/) implementaion which ever works best for cross platform
- [ ] GUI requires an installation of `Windows Subsystem for Linux`
  - [ ] [Cygwin](https://www.cygwin.com/) for Windows < 10
  - [ ] [Ubuntu](https://docs.microsoft.com/en-us/windows/wsl/install-win10) for Windows10
- [ ] `Additional features here`

#### References: 
- [miRge2.0](https://github.com/mhalushka/miRge)
- [miRge1.0](https://github.com/mhalushka/miRge-1) 
  - [sample_output](https://baraslab.github.io/miRge/miRge/miRge.exampleOutput/report.html)
- [sRNAtoolbox](https://bioinfo5.ugr.es/srnatoolbox/srnabench/)
- [Tools4miRs](https://tools4mirs.org/software/isomirs_identification/)
- [Chimira](http://wwwdev.ebi.ac.uk/enright-dev/chimira/index.php)
- [Beeware by Russell Keith-Magee](https://www.youtube.com/watch?v=qaPzlIJ57dk) 
- [ViennaRNA 2.0](https://github.com/ViennaRNA/ViennaRNA)
- [ ] `Additional features here`

## Upcoming events:
- [COVALENCE CONFERENCE 2020](https://www.covalenceconf.com/)
  - Date: January 24, 2020. (09:00AM to 07:00PM)
  - Venue:
    > Slack HQ, 
    > 500 Howard St,
    > San Francisco CA 94105.
