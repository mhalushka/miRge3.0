#!/usr/bin/env python
import time
import pandas as pd
from pathlib import Path
import numpy as np

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



def summarize(args, workDir, ref_db,base_names, pdMapped, sampleReadCounts, trimmedReadCounts, trimmedReadCountsUnique):
    """
    THIS FUNCTION IS CALLED FIRST FROM THE miRge3.0 to summarize the output.  
    """
    mirMergedNameDic={}
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
    with open(mergeFile, "r") as merge_file:
        for line in merge_file:
            line_content = line.strip().split(',')
            for item in line_content[1:]:
                mirMergedNameDic.update({item:line_content[0]})
    pdMapped = pdMapped.reset_index(level=['Sequence'])
    subpdMapped = pdMapped[(pdMapped['exact miRNA'].astype(bool) | pdMapped['isomiR miRNA'].astype(bool))]
    cannonical = pdMapped[pdMapped['exact miRNA'].astype(bool)]
    isomirs = pdMapped[pdMapped['isomiR miRNA'].astype(bool)]
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
    df = df.loc[(df.sum(axis=1) != 0)] # THIS WILL ELEMINATE ROWS ACCROSS SAMPLES WHO'S SUM IS ZERO
    Filtered_miRNA_Reads = df.sum(axis = 0, skipna = True)[base_names]
    Filtered_miRNA_Reads = Filtered_miRNA_Reads.to_dict()
    miR_RPM = (df.div(df.sum(axis=0))*1000000).round(4)
    miRNA_df = subpdMapped.groupby(['miRNA_cbind']).sum()[base_names]
    sumTotal = miRNA_df.sum(axis = 0, skipna = True)
    l_1d = sumTotal.to_dict()
    miRgefileToCSV = Path(workDir)/"miR.Counts.csv"
    miRgeRPMToCSV = Path(workDir)/"miR.RPM.csv"
    df.to_csv(miRgefileToCSV)
    miR_RPM.to_csv(miRgeRPMToCSV)
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
