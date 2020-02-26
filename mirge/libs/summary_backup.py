#!/usr/bin/env python
import time
import pandas as pd
from pathlib import Path
import numpy as np
import subprocess
from difflib import unified_diff, Differ
"""
THIS SCRIPT CONTAINS LOTS OF PANDAS FUNCTION TO DERIVE THE SUMMARY
IF YOU ARE A DEVELOPER, AND WANT TO UNDERSTAND THIS SCRIPT!! I WOULD RECOMMEND YOU TO BE THOROUGH WITH pandas FUNCTIONS 
THIS FUNCTION RETURNS THE FOLLOWING FILES AS OUTPUT:
    miR.Counts
    miR.RPM
    annotation.report
"""


def mirge_can(can, iso, df, ca_thr, file_name):
    """
    THIS FUNCTION TAKES Exact miRNA & isomiRs FOR EACH SAMPLE and CALCULATES CANNONICAL RATIO
    """
    # MIRGING TAKES PLACE FOR ALL THE ELEMENTS OF ALL exact miRNA LEAVING BEHIND isomiRs WHICH IS NOT IN exact miRNA
    merged_left = pd.merge(left=can,right=iso, how='left', left_on='exact miRNA', right_on='isomiR miRNA') 
    merged_left = merged_left.fillna(0) 
    # IN PANDAS, IF THE COLUMN NAME IS SAME, IT WILL REPLACE THE NAME WITH _x AND _y 
    file_nameX = str(file_name)+"_x" 
    file_nameY = str(file_name)+"_y"
    merged_left[file_nameY] = merged_left[file_nameY].astype(int)
    merged_left.loc[merged_left[file_nameX] < 2, [file_nameX]] = 0 # REPLACE THE exact miRNA COUNT WITH ZERO IF exact miRNA < 2
    merged_left.loc[merged_left[file_nameX] < 2, [file_nameY]] = 0 # REPLACE THE isomiR COUNT WITH ZERO IF exact miRNA < 2
    merged_left.loc[merged_left[file_nameY] > 0, 'ratio'] =  merged_left[file_nameX]/merged_left[file_nameY] #CREATES A COLUMN CALLED RATIO IF DENOMINATOR IS NOT ZERO
    merged_left.loc[merged_left[file_nameY] == 0, 'ratio'] =  merged_left[file_nameX] #ASSIGNS exact miRNA COUNTS AS RATIO IF DENOMINATOR IS ZERO
    cols = [file_nameX, file_nameY]
    merged_left[file_name] = merged_left.loc[merged_left['ratio'] > ca_thr, cols].sum(axis=1)
    merged_left[file_name].fillna(0, inplace=True)
    df = df.join(merged_left[file_name], how='outer')
    del merged_left
    return df


def create_gff(args, pre_mirDict, mirDict, d, filenamegff, cannonical, isomirs, base_names, ref_db, annotation_lib):
    pre_cols1 = ["Sequence","exact miRNA"]
    pre_cols2 = ["Sequence","isomiR miRNA"]
    cols1 = pre_cols1 + base_names
    cols2 = pre_cols2 + base_names
    can_gff_df = pd.DataFrame(cannonical, columns= cols1)
    iso_gff_df = pd.DataFrame(isomirs, columns= cols2)
    canonical_gff = can_gff_df.values.tolist()
    isomir_gff = iso_gff_df.values.tolist()
    gffwrite = open(filenamegff, "w+")
    gffwrite.write("# GFF3 adapted for miRNA sequencing data\n")
    gffwrite.write("## VERSION 0.0.1\n")
    version_db = "miRBase_22" if ref_db == "miRBase" else "MirGeneDB_2.0"
    gffwrite.write("## source-ontology: " + version_db + "\n")
    gffwrite.write("## COLDATA: "+ ",".join(str(nm) for nm in base_names) + "\n")
    # GFF3 adapted for miRNA sequencing data")
    start=end=0
    parent="-"
    filter_var = "Pass"
    strand = "+"
    dots = "."
    precursorSeq = uid_val = ""
    """
    READING ANNOTATION DATA TO GET GENOMIC COORDINATES AND PRECURSOR miRNA 
    """
    pre_cur_name=pre_strand={}
    mature_cor_start=mature_cor_end={}
    with open(annotation_lib) as alib:
        for annlib in alib:
            annlib = annlib.strip()
            annlib_list = annlib.split("\t")
            if ref_db == "MirGeneDB":
                if annlib_list[2] == "pre_miRNA":
                    pre_name = annlib_list[8].split(";")[0]
                    pre_name = pre_name.replace("ID=","")
                    pre_strand[pre_name] = annlib_list[6]
                else:
                    mature_name = annlib_list[8].split(";")[0]
                    mature_name = mature_name.replace("ID=","")
                    pre_cur_name[mature_name] = pre_name
                    mature_cor_start[mature_name] = annlib_list[3]
                    mature_cor_end[mature_name] = annlib_list[4]
            else:
                if annlib_list[2] == "miRNA_primary_transcript":
                    pre_name = annlib_list[8].split(";")[-1]
                    pre_name = pre_name.replace("Name=","")
                    pre_strand[pre_name] = annlib_list[6] 
                else:
                    mature_name = annlib_list[8].split(";")[2]
                    mature_name = mature_name.replace("Name=","")
                    mature_cor_start[mature_name] = annlib_list[3]
                    mature_cor_end[mature_name] = annlib_list[4]

    for cans in canonical_gff:
        seq_m = cans[0]
        if "." in cans[1]:
            seq_master = cans[1].split(".")[0]
        else:
            seq_master = cans[1]
        canonical_expression = ','.join(str(x) for x in cans[2:])
        new_string=""
        try:
            if seq_master not in mirDict:
                seq_master = seq_master.replace("-3p","") 
                seq_master = seq_master.replace("-5p","") 
                seq_master = seq_master.replace("-3p*","") 
                seq_master = seq_master.replace("-5p*","") 
                if ref_db == 'miRBase':
                    parent = seq_master
                else: 
                    parent = seq_master+"_Pre"
            else:
                if ref_db == 'miRBase':
                    seq_master_temp = seq_master.replace("-3p","")
                    seq_master_temp = seq_master.replace("-5p","")
                    parent = seq_master_temp
                else:
                    seq_master_temp = seq_master.replace("-3p","")
                    seq_master_temp = seq_master.replace("-5p","")
                    seq_master_temp = seq_master.replace("-3p*","")
                    seq_master_temp = seq_master.replace("-5p*","")
                    parent = seq_master_temp + "_Pre"
            master_seq = mirDict[seq_master]
            try:
#                new_parent = '-'.join(parent.split('-')[:-1])
#                precursorSeq = pre_mirDict[new_parent]
                parent = parent.lower()
                gffwrite.write(parent+"\n")
                precursorSeq = pre_mirDict[parent]
                pass
            except KeyError:
                print("Key not found: " + str(parent))
                pass

            if seq_m == master_seq:
                type_rna = "ref_miRNA"
                cigar = str(len(seq_m))+"M"
                if precursorSeq != "":
                    start = precursorSeq.find(seq_m) + 1
                    end = start + len(seq_m)
                else:
                    start = 1
                    end = start + len(seq_m)
                mi_var = seq_master+"\t"+version_db+"\t"+type_rna+"\t"+str(start)+"\t"+str(end)+"\t.\t+\t.\tRead "+seq_m+"; UID "+uid_val+"; Name "+ seq_master +"; Parent "+parent+"; Variant NA; Cigar "+cigar+"; Filter Pass; Expression "+ canonical_expression + "\n"
                gffwrite.write(mi_var)

                #hsa-miR-143-3p  miRBase22   ref_miRNA   61  81  .   +   .   
                #Read TGAGATGAAGCACTGTAGCTC; UID 7FTVqJj; Name hsa-miR-143-3p; Parent hsa-mir-143; Variant NA; Cigar 21M; Expression 8210; Filter Pass
            else:
                type_rna = "isomiR"
                result = list(d.compare(master_seq, seq_m))
                re_len = len(master_seq)
                variant="0"
                for idx, bases in enumerate(result):
                    if not bases.startswith("-"):
                        bases = bases.replace(" ","")
                        if bases.startswith("+"):
                            if idx < re_len:
                                variant = "1"
                            new_string += "["+ str(bases) +"]"
                        else:
                            new_string += bases
                    else:
                        try:
                            if bases[idx+1].startswith("+"):
                                pass
                            else:
                                new_string += "[-]"
                        except IndexError:
                            pass
                        print(seq_master +"\t"+ master_seq +"\t"+ seq_m +" "+ new_string+"\t"+variant)
        except KeyError:
            print(seq_m+"\t"+seq_master)
            pass
        pass

    for isos in isomir_gff:
        seq_i = isos[0]
        seq_master_iso = isos[1]
        isomirs_expression = ','.join(str(x) for x in isos[2:])
        pass
    pass

def summarize(args, workDir, ref_db,base_names, pdMapped, sampleReadCounts, trimmedReadCounts, trimmedReadCountsUnique):
    """
    THIS FUNCTION IS CALLED FIRST FROM THE miRge3.0 to summarize the output.  
    """
    ca_thr = float(args.crThreshold)
    mfname = args.organism_name + "_merges_" + ref_db + ".csv"
    mergeFile = Path(args.libraries_path)/args.organism_name/"annotation.Libs"/mfname
    if args.spikeIn:
        col_headers = ['hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','spike-in']
        col_vars = ['hmir','mtrna','pmtrna','snorna','rrna','ncrna','mrna','spikein']
    else:
        col_headers = ['hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA']
        col_vars = ['hmir','mtrna','pmtrna','snorna','rrna','ncrna','mrna']
    empty_list=dict() #Actually this is a dictionary, to collect dictionary of `sample names` as keys and `sum of expression` as values for each element of col_vars. Sorry for naming it _list. 
    for element, col_in in enumerate(col_headers):
        for file_name in base_names:
            if col_vars[element] in empty_list:
                empty_list[col_vars[element]].update({file_name:pdMapped[pdMapped[col_in].astype(bool)][file_name].values.sum()})
            else:
                empty_list[col_vars[element]] = {file_name:pdMapped[pdMapped[col_in].astype(bool)][file_name].values.sum()}

    
    """
    BEGINNING OF WORKING AROUND WITH EXACT miRNA and isomiRs 
    """
    mirMergedNameDic={}
    mirMergedDataframeDic={}
    with open(mergeFile, "r") as merge_file:
        for line in merge_file:
            line_content = line.strip().split(',')
            for item in line_content[1:]:
                mirMergedNameDic.update({item:line_content[0]})
                mirMergedDataframeDic.update({line_content[0]:"1"})
    
    #allSequences = pdMapped.index.shape[0]
    #print(allSequences)

    pdMapped = pdMapped.reset_index(level=['Sequence'])
    subpdMapped = pdMapped[(pdMapped['exact miRNA'].astype(bool) | pdMapped['isomiR miRNA'].astype(bool))]
    cannonical = pdMapped[pdMapped['exact miRNA'].astype(bool)]
    isomirs = pdMapped[pdMapped['isomiR miRNA'].astype(bool)]
    if args.gff_out:
        pre_mirDict = dict()
        mirDict = dict()
        filenamegff = workDir/"sample_miRge3.gff"
        maturefname = args.organism_name + "_mature_" + ref_db + ".fa"
        pre_fname = args.organism_name + "_hairpin_" + ref_db
        fasta_file = Path(args.libraries_path)/args.organism_name/"fasta.Libs"/maturefname
        precursor_file = Path(args.libraries_path)/args.organism_name/"index.Libs"/pre_fname
        annotation_lib = Path(args.libraries_path)/args.organism_name/"annotation.Libs"/args.organism_name+"_"+ref_db+".gff3"

        bwtCommand = Path(args.bowtie_path)/"bowtie-inspect" if args.bowtie_path else "bowtie-inspect"
        bwtExec = bwtCommand+" -a 20000 -e "+ str(precursor_file)
        bowtie = subprocess.run(str(bwtExec), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
        #READING PRECURSOR miRNA SEQUENCES INFORMATION IN A DICTIONARY (pre_mirDict)
        if bowtie.returncode==0:
            bwtOut = bowtie.stdout
            bwtErr = bowtie.stderr
        for srow in bwtOut.split('\n'):
            if '>' in srow:
                srow = srow.replace(">","")
                headmil = srow.split(" ")[0]
                if "MirGeneDB" in ref_db:
                    headmil = headmil.split("_")[0]
            else:
                #headmil = '-'.join(headmil.split('-')[:-1])
                pre_mirDict[headmil] = srow
        #READING MATURE miRNA SEQUENCES INFORMATION IN A DICTIONARY (mirDict)
        with open(fasta_file) as mir:
            for mil in mir:
                mil = mil.strip()
                if '>' in mil:
                    headmil_mi = mil.replace(">","")
                    if "MirGeneDB" in ref_db:
                        headmil_mi = headmil_mi.split("_")[0]
                else:
                    mirDict[headmil_mi] = mil
        d = Differ()
        create_gff(args, pre_mirDict, mirDict, d, filenamegff, cannonical, isomirs, base_names, ref_db, annotation_lib)
        
        #https://stackoverflow.com/questions/35125062/how-do-i-join-2-columns-of-a-pandas-data-frame-by-a-comma
        #ref_db #miRBase | MirGeneDB
        #bowtie-inspect -a 20000 -e human_hairpin_miRBase
        #/Libs/human/index.Libs/human_hairpin_miRBase
        #/Libs/human/fasta.Libs
        #human_mature_miRBase.fa
        pass
    if args.spikeIn:
        cannonical = cannonical.drop(columns=['Sequence','hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','annotFlag','isomiR miRNA','spike-in'])
        isomirs = isomirs.drop(columns=['Sequence','exact miRNA','hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','annotFlag','spike-in'])
        subpdMapped = subpdMapped.drop(columns=['Sequence','hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','annotFlag','spike-in'])
    else:
        cannonical = cannonical.drop(columns=['Sequence','hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','annotFlag','isomiR miRNA'])
        isomirs = isomirs.drop(columns=['Sequence','exact miRNA','hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','annotFlag'])
        subpdMapped = subpdMapped.drop(columns=['Sequence','hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','annotFlag'])

    subpdMapped['miRNA_cbind'] = subpdMapped[['exact miRNA', 'isomiR miRNA']].apply(lambda x: ''.join(x), axis = 1)
    subpdMapped['miRNA_fin'] = subpdMapped['miRNA_cbind'].map(mirMergedNameDic)
    subpdMapped = subpdMapped.fillna(0)
    subpdMapped.loc[subpdMapped.miRNA_fin == 0, 'miRNA_fin'] = subpdMapped.miRNA_cbind
    subpdMapped.set_index('miRNA_cbind',inplace = True)
    cannonical.set_index('exact miRNA',inplace = True)
    isomirs.set_index('isomiR miRNA',inplace = True)
    cann_collapse = cannonical.groupby(['exact miRNA']).sum()[base_names]
    iso_collapse = isomirs.groupby(['isomiR miRNA']).sum()[base_names]
    cann_collapse = cann_collapse.reset_index(level=['exact miRNA'])
    iso_collapse = iso_collapse.reset_index(level=['isomiR miRNA'])
    df = pd.DataFrame(cann_collapse['exact miRNA'].tolist(), columns = ['exact miRNA'])
    for file_name in base_names:
        df = mirge_can(cann_collapse, iso_collapse, df, ca_thr, file_name)
    
    df['miRNA'] = df['exact miRNA'].map(mirMergedNameDic)
    df = df.fillna(0)
    df.loc[df.miRNA == 0, 'miRNA'] = df['exact miRNA']
    df.set_index('miRNA',inplace = True)
    df.drop(columns=['exact miRNA'])
    df = df.groupby(['miRNA']).sum()[base_names]
    #df = df.loc[(df.sum(axis=1) != 0)] # THIS WILL ELEMINATE ROWS ACCROSS SAMPLES WHO'S SUM IS ZERO
    Filtered_miRNA_Reads = df.sum(axis = 0, skipna = True)[base_names]
    Filtered_miRNA_Reads = Filtered_miRNA_Reads.to_dict()
    miR_RPM = (df.div(df.sum(axis=0))*1000000).round(4)
    miRNA_df = subpdMapped.groupby(['miRNA_cbind']).sum()[base_names]
    sumTotal = miRNA_df.sum(axis = 0, skipna = True)
    l_1d = sumTotal.to_dict()
    miRgefileToCSV = Path(workDir)/"miR.Counts.csv"
    miRgeRPMToCSV = Path(workDir)/"miR.RPM.csv"
    indexName  = str(args.organism_name) + '_mirna_' + str(ref_db)
    indexFiles = Path(args.libraries_path)/args.organism_name/"index.Libs"/indexName
    bwtExec = "bowtie-inspect -n " + str(indexFiles)
    #bwtExec = "bowtie-inspect -n /home/arun/repositories/Project_120919/mirge/Libs/human/index.Libs/human_mirna_miRBase"
    bowtie = subprocess.run(str(bwtExec), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
    if bowtie.returncode==0:
        bwtOut = bowtie.stdout
        bwtErr = bowtie.stderr
        lines = bwtOut.strip()
        for srow in lines.split('\n'):
            if srow not in mirMergedNameDic:
                mirMergedDataframeDic.update({srow:"1"})

    mirMerged_df = pd.DataFrame(list(mirMergedDataframeDic.keys()),columns = ['miRNA']) #Contains all the miRNA including those that is not expressed
    mirMerged_df.set_index('miRNA',inplace = True)
    
    mirCounts_completeSet = mirMerged_df.join(df, how='outer').fillna(0)
    mirRPM_completeSet = mirMerged_df.join(miR_RPM, how='outer').fillna(0)
    #df.to_csv(miRgefileToCSV)
    #miR_RPM.to_csv(miRgeRPMToCSV)

    mirCounts_completeSet.to_csv(miRgefileToCSV)
    mirRPM_completeSet.to_csv(miRgeRPMToCSV)

    miRNA_counts={}
    trimmed_counts={}
    for file_name in base_names:
        numOfRows = df.index[df[file_name] > 0].shape[0]
        mapped_rows = pdMapped.index[pdMapped[file_name] > 0].shape[0]
        mirna_dict = {file_name:numOfRows}
        miRNA_counts.update(mirna_dict)

    """
    END OF WORKING AROUND WITH EXACT miRNA and isomiRs 
    """
    trimmed_counts = trimmedReadCountsUnique
    if args.spikeIn:
        pre_summary = {'Total Input Reads':sampleReadCounts,'Trimmed Reads (all)':trimmedReadCounts,'Trimmed Reads (unique)':trimmed_counts,'All miRNA Reads':l_1d,'Filtered miRNA Reads':Filtered_miRNA_Reads,'Unique miRNAs':miRNA_counts, 'Hairpin miRNAs':empty_list[col_vars[0]],'mature tRNA Reads':empty_list[col_vars[1]],'primary tRNA Reads':empty_list[col_vars[2]],'snoRNA Reads':empty_list[col_vars[3]],'rRNA Reads':empty_list[col_vars[4]],'ncRNA others':empty_list[col_vars[5]],'mRNA Reads':empty_list[col_vars[6]],'Spike-in':empty_list[col_vars[7]]}
        col_tosum = ['All miRNA Reads','Hairpin miRNAs','mature tRNA Reads','primary tRNA Reads','snoRNA Reads','rRNA Reads','ncRNA others','mRNA Reads','Spike-in']
        colRearrange = ['Total Input Reads', 'Trimmed Reads (all)','Trimmed Reads (unique)','All miRNA Reads','Filtered miRNA Reads','Unique miRNAs','Hairpin miRNAs','mature tRNA Reads','primary tRNA Reads','snoRNA Reads','rRNA Reads','ncRNA others','mRNA Reads','Spike-in','Remaining Reads']
    else:
        pre_summary = {'Total Input Reads':sampleReadCounts,'Trimmed Reads (all)':trimmedReadCounts,'Trimmed Reads (unique)':trimmed_counts,'All miRNA Reads':l_1d,'Filtered miRNA Reads':Filtered_miRNA_Reads,'Unique miRNAs':miRNA_counts, 'Hairpin miRNAs':empty_list[col_vars[0]],'mature tRNA Reads':empty_list[col_vars[1]],'primary tRNA Reads':empty_list[col_vars[2]],'snoRNA Reads':empty_list[col_vars[3]],'rRNA Reads':empty_list[col_vars[4]],'ncRNA others':empty_list[col_vars[5]],'mRNA Reads':empty_list[col_vars[6]]}
        col_tosum = ['All miRNA Reads','Hairpin miRNAs','mature tRNA Reads','primary tRNA Reads','snoRNA Reads','rRNA Reads','ncRNA others','mRNA Reads']
        colRearrange = ['Total Input Reads', 'Trimmed Reads (all)','Trimmed Reads (unique)','All miRNA Reads','Filtered miRNA Reads','Unique miRNAs','Hairpin miRNAs','mature tRNA Reads','primary tRNA Reads','snoRNA Reads','rRNA Reads','ncRNA others','mRNA Reads','Remaining Reads']

    summary = pd.DataFrame.from_dict(pre_summary).astype(int)
    summary['Remaining Reads'] = summary['Trimmed Reads (all)'] - (summary[col_tosum].sum(axis=1))
    summary = summary.reindex(columns=colRearrange)
    summary.index.name = "Sample name(s)"
    report = Path(workDir)/"annotation.report.csv"
    report_html = Path(workDir)/"annotation.report.html"
    summary.to_csv(report)
    summary = summary.reset_index(level=['Sample name(s)'])
    summary.index += 1

    data_in_html = summary.to_html(index=False)
    import re
    table = '<table style="color="black";font-size:15px; text-align:center; border:0.2px solid black; border-collapse:collapse; table-layout:fixed; height="550"; text-align:center">'
    th = '<th style ="background-color: #3f51b5; color:#ffffff; text-align:center">'
    td ='<td style="vertical-align: middle;background-color: #edf6ff;font-size: 14px;font-family: Arial;font-weight: normal;color: #000000;text-align:center; height="250";border:0; ">'
    data_in_html = data_in_html.replace("<table>", th)
    data_in_html = data_in_html.replace("<th>", th)
    data_in_html = data_in_html.replace("<td>", td)
    with open(report_html,'w') as f:
        f.write(data_in_html)
