
# Directory Structure
[miRge2](https://github.com/mhalushka/miRge/tree/master/src/mirge)
```
miRge/src/mirge/
```

```bash
├── classes
│   ├── __init__.py 
│   ├── readCluster.py
│   ├── readPrecursor.py
├── forgi
│   ├── __init__.py 
│   ├── graph
│   │   ├── __init__.py
│   │   ├── bulge_graph.py
│   ├── threedee
│   │   ├── __init__.py
│   │   ├── utilities
│   │   │   ├── __init__.py
│   │   │   ├── cytvec.py
│   │   │   ├── mcannotate.py
│   │   │   ├── vector.py
│   ├── utilities
│   │   ├── __init__.py
│   │   ├── debug.py
│   │   ├── stuff.py
├── models
│   ├── human_svc_model.pkl
│   ├── mouse_svc_model.pkl
│   ├── others_svc_model.pkl
│   ├── total_features_namelist.txt
├── rScripts
│   ├── A-to-I_plot.R
├── utils
│   ├── __init__.py
│   ├── cluster_basedon_location.py
│   ├── convert2Fasta.py
│   ├── extractPreMiRName.py
│   ├── fileExist.py
│   ├── filter.py
│   ├── generateReport.py
│   ├── generate_Seq.py
│   ├── generate_featureFiles.py
│   ├── get_precursors.py
│   ├── miRNAmerge.py
│   ├── model_predict.py
│   ├── parseArgument.py
│   ├── preTrimClusteredSeq.py
│   ├── preprocess_featureFiles.py
│   ├── processSam.py
│   ├── quantReads.py
│   ├── renameStrFile.py
│   ├── runAnnotationPipeline.py
│   ├── screen_precusor_candidates.py
│   ├── summarize.py
│   ├── trim_file.py
│   ├── writeDataToCSV.py
│   └── write_novel_report.py
├── __init__.py
└── __main__.py
```
___

**├── utils**

**ANNOTATE**
- [x] ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) extractPreMiRName.py (63)
- [ ] ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) generateReport.py (281)
- [x] ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) filter.py (32)
- [x] ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) miRNAmerge.py (43)
- [x] ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) quantReads.py (46)
- [x] ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) runAnnotationPipeline.py (707)
- [x] ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) trim_file.py (138)
- [ ] ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) summarize.py (66)
- [ ] ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) writeDataToCSV.py (1615)

**PREDICT**
- [ ] ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) cluster_basedon_location.py
- [x] ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) convert2Fasta.py
- [ ] ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) fileExist.py
- [ ] ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) generate_Seq.py
- [ ] ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) generate_featureFiles.py
- [ ] ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) get_precursors.py
- [ ] ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) model_predict.py
- [ ] ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) parseArgument.py
- [ ] ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) preTrimClusteredSeq.py
- [ ] ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) preprocess_featureFiles.py
- [ ] ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) renameStrFile.py
- [ ] ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) screen_precusor_candidates.py
- [ ] ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) write_novel_report.py
- [ ] ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) processSam.py

 ![#1589F0](https://placehold.it/15/1589F0/000000?text=+) `Yet to work`
  
 ![#c5f015](https://placehold.it/15/c5f015/000000?text=+) `Complemented the code - no longer needed`
 ___
This Directory structure github markdown is implemented form [@Aerobatic](https://github.com/aerobatic/markdown-content/blob/master/docs/directory-structure.md)

