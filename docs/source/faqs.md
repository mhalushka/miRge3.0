# Frequently asked questions (FAQ) 


_We are very greatful and also thankful to all the users of miRge3.0 who rasied GitHub issues in the past that helped us solve few technical problems and improve miRge3.0 functionality further. We expect continued support towards this project. Here we have gathered a few frequently asked questions over the period regarding technical as well as biological/scientific questions. I hope this documentation will be useful as a ready response/solution for your queries._

*Before getting started please note; if you don't find a solution to your query in this page then create a new issue and we will get back to you at the earliest. Describe the `Title` to include the error you are facing e.g., `numbpy type error` and in the `Comment` section, it would be best if you could put the command line used, followed by the whole error. (You can delete your file names if you prefer).*

##### ***How to create an issue?***
> Click [create new issue](https://github.com/mhalushka/miRge3.0/issues/new) and in Title: "Please describe the error you think is obvious and will be general for the scientific community to recognize", and Comment: "Give us the maximum information possible regarding the error that you can see on the standard output/terminal"


## Frequent questions raised on GitHub:

1. [How to use Unique Molecular Identifiers (UMIs)?](https://mirge3.readthedocs.io/en/latest/faqs.html#how-to-use-unique-molecular-identifiers-umis)
2. [TypeError: Cannot interpret <attribute 'dtype' of 'numpy.generic' objects> as a data type](https://mirge3.readthedocs.io/en/latest/faqs.html#typeerror-cannot-interpret-attribute-dtype-of-numpy-generic-objects-as-a-data-type)
3. [UnsatisfiableError: bowtie=1.3.0 -> libgcc-ng[version='>=9.3.0'] -> __glibc[version='>=2.17']](https://mirge3.readthedocs.io/en/latest/faqs.html#unsatisfiableerror-bowtie-1-3-0-libgcc-ng-version-9-3-0-glibc-version-2-17)
4. [Is there any way to skip the adaptor trimming process? and how to determine adapter sequence of a Run?](https://mirge3.readthedocs.io/en/latest/faqs.html#is-there-any-way-to-skip-the-adaptor-trimming-process-and-how-to-determine-adapter-sequence-of-a-run)
5. [How to use and tweak data with Spike-in expirements?](https://mirge3.readthedocs.io/en/latest/faqs.html#how-to-use-and-tweak-data-with-spike-in-expirements)
6. [How to use -dex DESeq2 analysis?](https://mirge3.readthedocs.io/en/latest/faqs.html#how-to-use-dex-deseq2-analysis)
7. [What is the threshold of the proportion of canonical reads (-ex, --crThreshold)?](https://mirge3.readthedocs.io/en/latest/faqs.html#what-is-the-threshold-of-the-proportion-of-canonical-reads-ex-crthreshold)
8. [How to input paired-end sequencing data?](https://mirge3.readthedocs.io/en/latest/faqs.html#how-to-input-paired-end-sequencing-data)


##### ***How to use Unique Molecular Identifiers (UMIs)?*** 

>A detailed documentation for UMI test run is available [here](https://mirge3.readthedocs.io/en/latest/quick_start.html#running-samples-with-umi). miRge3.0 is designed to process UMIs for Illumina and Qiagen. The parameters to trim UMIs and removing PCR duplicates are different, and also, selecting Qiagen UMI needs an additional parameter. 
>
>These following issues were raised:
>
> - [#32 (comment)](https://github.com/mhalushka/miRge3.0/issues/32#issue-1149944971)<br/>
> - [#46 (comment)](https://github.com/mhalushka/miRge3.0/issues/46#issue-1273723168)<br/>
> - [#28 (comment)](https://github.com/mhalushka/miRge3.0/issues/28#issue-1077400071)<br/>


##### ***TypeError: Cannot interpret <attribute 'dtype' of 'numpy.generic' objects> as a data type***


> I suspect there is a conflict with pandas and numpy in your local machine, I want you to upgrade pandas and try the command again. You can upgrade it as shown (python3.7 if you are using py37) in the following issues:
>
>[#20 (comment)](https://github.com/mhalushka/miRge3.0/issues/20#issuecomment-942408755)
>[#47 (comment)](https://github.com/mhalushka/miRge3.0/issues/47#issue-1285035891)

##### ***UnsatisfiableError: bowtie=1.3.0 -> libgcc-ng[version='>=9.3.0'] -> __glibc[version='>=2.17']***

> The discussion on this issue is available in the following GitHub issue. Thank you [@asucrer](https://github.com/asucrer), for providing solution.
> 
> [#31 (comment)](https://github.com/mhalushka/miRge3.0/issues/31#issuecomment-1027858446) 
```
Solution suggested by the user @asucrer, please follow the steps:

conda create -n mirge     # IMPORTANT to not specify the python version in this step 
source activate mirge
conda install -c bioconda mirge3   # Every dependency (including python) is installed
conda install -c bioconda tbb=2020.2    # Solves issue associated to Bowtie installation
conda install -c bioconda openssl=1.0   # Solves issue associated to Samtools installation
```

##### ***Is there any way to skip the adaptor trimming process? and how to determine adapter sequence of a Run?*** 

>miRge3.0 allows users to skip the adapter trimming step, and there are several options on how to provide adapter sequences and the following issue provide a list of adapter sequences for various platforms. [Curation date: January 2020]. 
>
> - [#20 (comment)](https://github.com/mhalushka/miRge3.0/issues/20#issue-1023411842)
> - [#20 (looped - comment)](https://github.com/mhalushka/miRge3.0/issues/20#issuecomment-943970154)
>
>Please NOTE: To trim adapter sequences at both ends please follow the documentation [Linked-adapters](https://mirge3.readthedocs.io/en/master/quick_start.html#trimming-both-5-and-3-adapters-linked-adapters)

##### ***How to use and tweak data with Spike-in expirements?***  

> An example usage of spike-in libraries and how to add/append spike-in reads of interest to the existing libraries and interpretation is described in the following issues:
>
> - [#27 (comment)](https://github.com/mhalushka/miRge3.0/issues/27#issue-1074611243) 
> - [#48 (comment)](https://github.com/mhalushka/miRge3.0/issues/48#issue-1288527770)


##### ***How to use -dex DESeq2 analysis?*** 

>The documentation for DESeq2 based differentiall expression analysis is available [here](https://mirge3.readthedocs.io/en/master/quick_start.html#performing-differential-expression-analysis)
>
> The following GitHub issues were raised:
> - [#41 (comment)](https://github.com/mhalushka/miRge3.0/issues/41#issue-1228345472)
> - [#33 (comment)](https://github.com/mhalushka/miRge3.0/issues/33#issue-1157209493)


##### ***What is the threshold of the proportion of canonical reads (-ex, --crThreshold)?*** 

>This was answered to an issue on why default value of 0.1 was chosen for --crThreshold in the following issue.
> - [#23 (comment)](https://github.com/mhalushka/miRge3.0/issues/23#issue-1045291352)
> - [#34 (comment)](https://github.com/mhalushka/miRge3.0/issues/34#issuecomment-1063675353)

##### ***How to input paired-end sequencing data?*** 

> miRge3.0 doesn't annotate paired-end data. 
> - [#7 (comment)](https://github.com/mhalushka/miRge3.0/issues/7#issue-737411998)
