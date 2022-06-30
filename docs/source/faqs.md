# Frequently asked question (FAQ) 

## Frequent questions raised on GitHub:

1. How to use Unique Molecular Identifiers (UMIs) 

<div style="background-color: #f8f8f8; border: 1px solid #e1e4e5; margin: 1px 0 24px">
A. A detailed documentation for UMI test run is available [here](https://mirge3.readthedocs.io/en/latest/quick_start.html#running-samples-with-umi) 
  miRge3.0 is designed to process UMIs for Illumina and Qiagen. The parameters to trim UMIs and removing PCR duplicates are different, and also, selecting Qiagen UMI needs an additional parameter. This following issues were raised:

[#32 (comment)](https://github.com/mhalushka/miRge3.0/issues/32#issue-1149944971)<br/>
[#46 (comment)](https://github.com/mhalushka/miRge3.0/issues/46#issue-1273723168)<br/>
[#28 (comment)](https://github.com/mhalushka/miRge3.0/issues/28#issue-1077400071)<br/>

</div>

2. TypeError: Cannot interpret <attribute 'dtype' of 'numpy.generic' objects> as a data type


```
A. I suspect there is a conflict with pandas and numpy in your local machine, I want you to upgrade pandas and try the command again. You can upgrade it as shown (python3.7 if you are using py37) in the following issues:

[#20 (comment)](https://github.com/mhalushka/miRge3.0/issues/20#issuecomment-942408755)
[#47 (comment)](https://github.com/mhalushka/miRge3.0/issues/47#issue-1285035891)

```
3. Is there any way to skip the adaptor trimming process? and how to determine adapter sequence of a Run?

```
A. miRge3.0 allows users to skip the adapter trimming step, and there are several options on how to provide adapter sequences and the following issue provide a list of adapter sequences for various platforms. [Curation date: January 2020]. 

[#20 (comment)](https://github.com/mhalushka/miRge3.0/issues/20#issue-1023411842)
[#20 (looped - comment)](https://github.com/mhalushka/miRge3.0/issues/20#issuecomment-943970154)

Please NOTE: To trim adapter sequences at both ends please follow the documentation [Linked-adapters](https://mirge3.readthedocs.io/en/master/quick_start.html#trimming-both-5-and-3-adapters-linked-adapters)
```

4. How to use and tweak data with Spike-in expirements? 

```
A. An example usage of spike-in libraries and how to add/append spike-in reads of interest to the existing libraries and interpretation is described in the following issues:

[#27 (comment)](https://github.com/mhalushka/miRge3.0/issues/27#issue-1074611243) 
[#48 (comment)](https://github.com/mhalushka/miRge3.0/issues/48#issue-1288527770)
```

5. how to use -dex DESeq2 analysis?

```
A. The documentation for DESeq2 based differentiall expression analysis is available [here](https://mirge3.readthedocs.io/en/master/quick_start.html#performing-differential-expression-analysis)
The following GitHub issues were raised:
[#41 (comment)](https://github.com/mhalushka/miRge3.0/issues/41#issue-1228345472)
[#33 (comment)](https://github.com/mhalushka/miRge3.0/issues/33#issue-1157209493)
```

6. What is the threshold of the proportion of canonical reads (-ex, --crThreshold)?

```
A. This was answered to an issue on why default value of 0.1 was chosen for --crThreshold in the following issue.
[#23 (comment)](https://github.com/mhalushka/miRge3.0/issues/23#issue-1045291352)
[#34 (comment)](https://github.com/mhalushka/miRge3.0/issues/34#issuecomment-1063675353)
```

7. How to input paired-end sequencing data?

```
A. miRge3.0 doesn't annotate paired-end data. 
[#7 (comment)](https://github.com/mhalushka/miRge3.0/issues/7#issue-737411998)
```


Copyright (c) 2020 Arun H. Patil and Marc K. Halushka
