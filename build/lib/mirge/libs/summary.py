#!/usr/bin/env python
import time
import math
import re
import pandas as pd
from pathlib import Path
import numpy as np
import subprocess
from difflib import unified_diff, Differ
from mirge.libs.miRgeEssential import UID
from mirge.libs.bamFmt import sam_header, bow2bam, createBAM
from mirge.libs.mirge2_tRF_a2i import trna_deliverables, a2i_editing
import os, sys
from mirge.classes.exportHTML import FormatJS
"""
THIS SCRIPT CONTAINS LOTS OF PANDAS FUNCTION TO DERIVE THE SUMMARY (EXCEPT FOR GFF-FUNCTION)
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


def create_gff(args, pre_mirDict, mirDict, d, filenamegff, cannonical, isomirs, base_names, ref_db, annotation_lib, workDir, mirRPM_completeSet):
    cols1 = ["Sequence","exact miRNA"] + base_names 
    cols2 = ["Sequence","isomiR miRNA"] + base_names
    canonical_gff = pd.DataFrame(cannonical, columns= cols1).values.tolist() # Gives list of list containg Sequence, miRNA name, expression values for the samples - ref miRNA
    isomir_gff = pd.DataFrame(isomirs, columns= cols2).values.tolist() # Gives list of list containg Sequence, miRNA name, expression values for the samples - isomiR 
    #can_gff_df = pd.DataFrame(cannonical, columns= cols1) # Gives list of list containg Sequence, miRNA name, expression values for the samples - ref miRNA
    #iso_gff_df = pd.DataFrame(isomirs, columns= cols2) # Gives list of list containg Sequence, miRNA name, expression values for the samples - isomiR 
    #canonical_gff = can_gff_df.values.tolist() 
    #isomir_gff = iso_gff_df.values.tolist()
    canonical_gff.extend(isomir_gff)  # APPENDING THE LIST OF ISOMIRS TO CANONICAL # Making one big list to get coordinates and anntations for GFF3 format of miRTop
    gffwrite = open(filenamegff, "w+") # creating file to write the gff output # sample_miRge3.gff
    gffwrite.write("# GFF3 adapted for miRNA sequencing data\n")
    gffwrite.write("## VERSION 0.0.1\n")
    version_db = "miRBase22" if ref_db == "miRBase" else "MirGeneDB2.0"
    gffwrite.write("## source-ontology: " + version_db + "\n")
    gffwrite.write("## COLDATA: "+ ",".join(str(nm) for nm in base_names) + "\n")
    #forJSgffData = open(str(Path(workDir)/"gff_data.csv"), "w+")
    JS_hmap_iso_miRcounts = dict()
    JS_hmap_iso_miRVar = dict()
    JS_hmap_list_ref = dict()

    JS_filenames = ",".join(str(nm) for nm in base_names)
    #JS_headername = 'Seq,miR,var,'+JS_filenames+',Total\n'
    #forJSgffData.write(JS_headername)
    # GFF3 adapted for miRNA sequencing data")
    start=0
    end=0
    parent="-"
    filter_var = "Pass"
    strand = "+"
    dots = "."
    precursorSeq = uid_val = ""
    """
    READING ANNOTATION DATA TO GET GENOMIC COORDINATES AND PRECURSOR miRNA 
    """

    pre_cur_name={}
    pre_strand={}
    mature_chromosome={}
    mature_cor_start={}
    mature_cor_end={}
    mature_strand = {}
    with open(annotation_lib) as alib: # Reading annotations GTF from miRBase or miRGeneDB based on user and gather coordinates and name and sequence of the precursor miRNA 
        for annlib in alib:
            annlib = annlib.strip()
            annlib_list = annlib.split("\t")
            try:
                if ref_db == "MirGeneDB":
                    if annlib_list[2] == "pre_miRNA":
                        pre_name = annlib_list[8].split(";")[0]
                        pre_name = pre_name.replace("ID=","")
                        pre_strand[pre_name] = annlib_list[6]
                        #print(pre_name)
                    else:
                        mature_name = annlib_list[8].split(";")[0]
                        mature_name = mature_name.replace("ID=","")
                        #print(mature_name)
                        if mature_name not in pre_cur_name:
                            pre_cur_name[mature_name] = pre_name
                            mature_chromosome[mature_name] = annlib_list[0]
                            mature_cor_start[mature_name] = annlib_list[3]
                            mature_cor_end[mature_name] = annlib_list[4]
                            mature_strand[mature_name] =  annlib_list[6] # Genomic strand 
                else:
                    if annlib_list[2] == "miRNA_primary_transcript":
                        pre_name = annlib_list[8].split(";")[-1]
                        pre_name = pre_name.replace("Name=","")
                        pre_strand[pre_name] = annlib_list[6] 
                    else:
                        mature_name = annlib_list[8].split(";")[2]
                        mature_name = mature_name.replace("Name=","")
                        if mature_name not in pre_cur_name:
                            pre_cur_name[mature_name] = pre_name
                            mature_chromosome[mature_name] = annlib_list[0] # Chromosome location
                            mature_cor_start[mature_name] = annlib_list[3] # Genomic coordinates and not miRNA seq to precursor sequence
                            mature_cor_end[mature_name] = annlib_list[4] # Genomic coordinates of miRNA and not its position w.r.t precursor sequence 
                            mature_strand[mature_name] =  annlib_list[6] # Genomic strand 
            except IndexError:
                pass
    #print(pre_cur_name)
    # bam_can_dict={}
    # bam_expression_dict={}
    JS_variantType_dataDict = dict()
    for cans in canonical_gff:
        gen_start=0
        gen_end=0
        seq_m = cans[0] # Sequence from datasest/pandas 
         
        if "." in cans[1]:
            seq_master = cans[1].split(".")[0] # miRNA name with .SNP extension
        else:
            seq_master = cans[1] # miRNA name 
        canonical_expression = ','.join(str(x) for x in cans[2:]) # Expression/collapsed counts for each sample - joining by ','
        JS_exprn_total = str(sum(cans[2:]))
        # bam_expression_dict[seq_m] = [int(x) for x in cans[2:]]
        new_string=""
        #print(seq_master)
        try:
            if seq_master in pre_cur_name:
                master_seq = mirDict[seq_master] # Fetch mature miRNA sequence 
                req_precursor_name = pre_cur_name[seq_master] # Fetch name of the corresponding precursor miRNA name
                gen_chr = mature_chromosome[seq_master]
                gen_start = int(mature_cor_start[seq_master])
                gen_end = int(mature_cor_end[seq_master])
                gen_strand = mature_strand[seq_master] 
                #print(master_seq, req_precursor_name, gen_chr, gen_start, gen_end, gen_strand)
            else:
                seq_master = seq_master.replace("-3p","")
                seq_master = seq_master.replace("-5p","")
                seq_master = seq_master.replace("-3p*","")
                seq_master = seq_master.replace("-5p*","")
                #print(seq_master)
                master_seq = mirDict[seq_master] # Fetch mature miRNA sequence
                req_precursor_name = pre_cur_name[seq_master] # Fetch name of the corresponding precursor miRNA name
                #print(req_precursor_name)
                gen_chr = mature_chromosome[seq_master]
                #print(gen_chr)
                gen_start = int(mature_cor_start[seq_master])
                gen_end = int(mature_cor_end[seq_master])
                gen_strand = mature_strand[seq_master] 
                #print(gen_start, gen_end, gen_strand)
                #print(master_seq, req_precursor_name, gen_chr, gen_start, gen_end, gen_strand)
            #print(req_precursor_name)
            precursorSeq = pre_mirDict[req_precursor_name] # # Fetch sequence of the corresponding precursor miRNA 
            #print(precursorSeq)

            if precursorSeq != "":
                start = precursorSeq.find(master_seq) + 1 # This is to calculate coordinates of the mature miRNA seq wrt precursor 
                end = start + len(master_seq) - 1
            else:
                start = 1 # Well, if the coordinates are not availble, I will put them as 1 to length of seq. Although this case is rare to none 
                end = start + len(master_seq) - 1
            #print()
            #print(precursorSeq+"\n"+master_seq+"\n")
            #print("2: "+ seq_m + " "+ str(start)+ " " + str(gen_start))
            if seq_m == master_seq: # If mature miRNA sequence is same as query FASTQ sequence, then it is reference miRNA
                type_rna = "ref_miRNA"
                if "N" not in seq_m: 
                    uid_val = UID(seq_m, "ref") # It is the Unique identifier, this function is present miRgeEssential.py 
                else:
                    uid_val = "." # If the sequence contains ambiguous bases (N), then UID is replaced by dot (.)
                cigar = str(len(seq_m))+"M" # CIGAR format for ref_miRNA is complete match, i.e., length of miRNA and M for example 25M (meaning 25 bases match)
                # Finally to write GFF output for ref_miRNA 

                mi_var = seq_master+"\t"+version_db+"\t"+type_rna+"\t"+str(start)+"\t"+str(end)+"\t.\t+\t.\tRead="+seq_m+"; UID="+uid_val+"; Name="+ seq_master +"; Parent="+req_precursor_name+"; Variant=NA; Cigar="+cigar+"; Expression="+canonical_expression +"; Filter=Pass; Hits="+ canonical_expression + "\n"
                #print(mi_var)
                gffwrite.write(mi_var)
                #JS_hmap_list_ref.append([seq_master, "ref", canonical_expression])
                #JS_hmap_list_ref.append([seq_master, "ref"] + canonical_expression.split(","))
                try:
                    JS_hmap_list_ref[seq_master].append(canonical_expression.split(","))
                except KeyError:
                    JS_hmap_list_ref[seq_master] = canonical_expression.split(",")
                #forJSgffData.write(str(seq_m)+","+seq_master+",ref,"+str(canonical_expression)+","+JS_exprn_total+"\n")
                # bam_can_dict[seq_m] = str(gen_chr) +"\t"+ str(gen_start) +"\t"+ cigar + "\t"+ gen_strand 
                    #bow2bam
                """
                FOR SAM/BAM FILE FORMAT: gen_chr, gen_start, gen_end
                """
            else: # If mature miRNA sequence is same as query FASTQ sequence, then it is reference miRNA else it is an isomiR
                type_rna = "isomiR"
                if "N" not in seq_m:
                    uid_val = UID(seq_m, "iso")
                else:
                    uid_val = "."
                result = list(d.compare(master_seq, seq_m)) # Python function difflib - Differ to detect changes between two strings
                re_len = len(master_seq)
                variant=""
                master_variant = []
                #print()
                #print("--**START**--")
                master_seq_bc = list(master_seq)
                result_seq_bc = result
                for mdx, mbases in enumerate(master_seq_bc):
                    if result_seq_bc[mdx].startswith("-"):
                        pass # Thought of something and left it as placeholder 
                    elif result_seq_bc[mdx].startswith("+"):
                        master_seq_bc.insert(mdx,'-') # Inserting '-' in the reference sequence if there is a change is base or insertion in the sequence 
                result_seq_bc = [bc.replace(" ", "") for bc in result_seq_bc if not bc.startswith("-")] # Cleaning the results and removing extra spaces obtained from difflib - Differ
                #print("--**MID**--")
                sub=[]
                for idx, bases in enumerate(result): # Creating an arrays if the reference now starts with '-', append "_" at that index position, etc and making two arrays of same lenght with variations
                    if bases.startswith("-"):
                        sub.append("_")
                    elif bases.startswith("+"):
                        sub.append(bases.replace(" ", ""))
                    else:
                        sub.append(bases.replace(" ",""))
                
                diff = len(sub) - len(master_seq_bc)
                for x in range(diff):
                    master_seq_bc.append("-")
                # print(master_seq_bc)
                # print(sub)
                for yidx, ys in enumerate(master_seq_bc): # <FORWARD> For upto 2 bases, find variants/base changes as shown in below comment A>T,T>G,T>G and basically delete '-' and "_" from two arrays 
                    if yidx > 0: 
                        try:
                            if ys == "-" and sub[yidx-1] == "_":
                                if yidx-2 > 0 and sub[yidx-2] == "_" and master_seq_bc[yidx+1] == "-":
                                    del master_seq_bc[yidx:yidx+2]
                                    del sub[yidx-2:yidx]
                #['A', 'A', 'A', 'C', 'C', 'G', 'T', 'T', 'A',  '-', 'C', 'C', 'A', 'T', 'T', 'A', 'C', 'T', 'G', 'A', 'G', 'T', 'T', '-',   '-']
                #['A', 'A', 'A', 'C', 'C', 'G', 'T', 'T', '_', '+T', 'C', 'C', 'A', 'T', 'T', 'A', 'C', 'T', 'G', '_', 'G', '_', '_', '+G', '+G']
                                else:
                                    temp_1 = master_seq_bc.pop(yidx)
                                    temp_2 = sub.pop(yidx-1)
                                #AAACCGTTTCCATTACTGGGG - This is the output of the loop
                #['A', 'A', 'A', 'C', 'C', 'G', 'T', 'T',  'A', 'C', 'C', 'A', 'T', 'T', 'A', 'C', 'T', 'G', 'A', 'G',  'T',  'T']
                #['A', 'A', 'A', 'C', 'C', 'G', 'T', 'T', '+T', 'C', 'C', 'A', 'T', 'T', 'A', 'C', 'T', 'G', '_', 'G', '+G', '+G']
                        except IndexError:
                            pass

                for yidx, ys in enumerate(master_seq_bc): # <REVERSE> For upto 2 bases, find variants/base changes as shown in below comment C>T,A>_ and basically delete '-' and "_" from two arrays 
                    if yidx > 0: 
                        try:
                            if ys == "-" and sub[yidx+1] == "_":
                                if yidx+2 <= len(sub) and sub[yidx+2] == "_" and master_seq_bc[yidx+1] == "-":
                                    del master_seq_bc[yidx:yidx+2]
                                    del sub[yidx:yidx+2]
                #['-', 'A', 'A', 'C', 'G', 'G', 'C', 'A', 'A', 'T', 'G', 'A', 'C', 'T', 'T', 'T', 'T', 'G', 'T', 'A', 'C', '-', 'C', 'A']
                #['+A', 'A', 'A', 'C', 'G', 'G', 'C', 'A', 'A', 'T', 'G', 'A', 'C', 'T', 'T', 'T', 'T', 'G', 'T', 'A', 'C', '+T', '_', '_']
                                else:
                                    temp_1 = master_seq_bc.pop(yidx)
                                    temp_2 = sub.pop(yidx+1)
                                # hsa-miR-548al   AACGGCAATGACTTTTGTACCA  AAACGGCAATGACTTTTGTACT
                #['-', 'A', 'A', 'C', 'G', 'G', 'C', 'A', 'A', 'T', 'G', 'A', 'C', 'T', 'T', 'T', 'T', 'G', 'T', 'A', 'C', 'C', 'A']
                #['+A', 'A', 'A', 'C', 'G', 'G', 'C', 'A', 'A', 'T', 'G', 'A', 'C', 'T', 'T', 'T', 'T', 'G', 'T', 'A', 'C', '+T', '_']
                        except IndexError:
                            pass
                #print(master_seq_bc)
                #print(sub)
                """
                    ['T', 'T', 'T', 'T', 'T', 'C', 'A', 'T', 'T', 'A', 'T', 'T', 'G', 'C', '-', 'T', 'C', 'C', 'T', 'G', 'A', 'C', '-', 'C'] =>  "-" in this line means insertion (Ref)
                    ['T', 'T', 'T', 'T', 'T', 'C', 'A', 'T', 'T', 'A', 'T', 'T', 'G', '_', '+G', 'T', 'C', 'C', 'T', 'G', '_', 'C', '+T', 'C'] => "_" in this line means deletion (Query)
                    hsa-miR-335-3p  TTTTTCATTATTGCTCCTGACC  TTTTTCATTATTGGTCCTGCTC
                """
                #print(seq_master +"\t"+ master_seq +"\t"+ seq_m +" "+ new_string+"\t"+variant)
                iso_add={} # Insertions
                iso_del={} # Deletions
                iso_sub={} # Substitutions
                iso_5p_add=iso_5p_del=""
                iso_3p_add=iso_3p_del=""
                #print("2"+ seq_m + " "+ str(start)+ " " + str(gen_start))
                for pidx, v in enumerate(master_seq_bc):
                    if v == "-": # Insertions
                        iso_add[pidx] = sub[pidx]
                    elif sub[pidx] == "_": # Deletions
                        iso_del[pidx] = v
                    elif v != sub[pidx]: # Substitutions
                        iso_sub[pidx] = sub[pidx]
                limit_master = len(master_seq_bc)
                ## Loop to detect 5p changes ##
                for i in range(limit_master):
                    if i in iso_add:
                        iso_5p_add += iso_add[i]
                    elif i in iso_del:
                        iso_5p_del += iso_del[i]
                    else: 
                        break 
                ## Loop to detect 3p changes ##
                for i in range(limit_master, -1, -1):
                    if i-1 in iso_add:
                        iso_3p_add += iso_add[i-1]
                    elif i-1 in iso_del:
                        iso_3p_del += iso_del[i-1]
                    else: 
                        break 
                ## Loop to detect internal changes ##
                ## Find and trim all the 5p and 3p changes to retain only the internal variants ##
                #cigar5padd = cigar5pdel = len3padd = cigar3pdel ="" # variant=iso_3p:-1; Cigar=21M; 
                variant = ""
                #print(precursorSeq)
                #print(precursorSeq[start-1:end])
                #print(seq_m)
                if iso_5p_add != "":
                    a5p = iso_5p_add
                    a5p = a5p.replace("+","")
                    len5padd = len(a5p)
                    pre_5pAdd = list(precursorSeq[start-len5padd-1:start-1])
                    try: 
                        template5p=non_template5p=0
                        for e5pidx, each5ps in enumerate(a5p):
                            if each5ps == pre_5pAdd[e5pidx]:
                                template5p += 1
                            else:
                                non_template5p += 1
                        if template5p != 0:
                            variant+= "iso_5p:+"+str(template5p)+","
                        if non_template5p != 0:
                            variant+= "iso_add5p:+"+str(non_template5p)+","
                    except IndexError:
                        variant+= "iso_5p:-"+str(len5padd)+","
                    start = start - len5padd
                    if gen_strand != "-":
                        #gen_start = gen_start + len5pdel
                        gen_start = gen_start - len5padd
                    # If the addition is a SNV w.r.t precursor, then it will be => iso_add5p:N. Number of non-template nucleotides added at 3p. 
                    #del sub[0:len5padd]
                    #del master_seq_bc[0:len5padd]
                    #print("5p_add:" + a5p + "Len"+ str(len5padd))
                    #cigar5padd = str(len5padd)+"I"
                if iso_5p_del != "":
                    d5p = iso_5p_del
                    len5pdel = len(d5p)
                    variant+= "iso_5p:+"+str(len5pdel)+","
                    start = start + len5pdel
                    if gen_strand != "-":
                        gen_start = gen_start + len5pdel
                    #del sub[0:len5pdel]
                    #del master_seq_bc[0:len5pdel]
                    #print("5p_del:" + d5p +"Len"+ str(len5pdel))
                    #cigar5pdel = str(len5pdel)+"D"
                if iso_3p_add != "":
                    a3p = "".join(iso_3p_add[::-1])
                    a3p = a3p.replace("+","")
                    len3padd = len(a3p)
                    #variant+= "iso_3p:+"+str(len3padd)+","
                    pre_3pAdd = list(precursorSeq[end:end+len3padd])
                    try: 
                        template3p=non_template3p=0
                        for e3pidx, each3ps in enumerate(a3p):
                            if each3ps == pre_3pAdd[e3pidx]:
                                template3p += 1
                            else:
                                non_template3p += 1
                        if template3p != 0:
                            variant+= "iso_3p:+"+str(template3p)+","
                        if non_template3p != 0:
                            variant+= "iso_add3p:+"+str(non_template3p)+","
                    except IndexError:
                        variant+= "iso_3p:+"+str(len3padd)+","
                    #iso_add3p
                    end = end + len3padd
                    gen_end = gen_end + len3padd
                    if gen_strand == "-":
                        gen_start-=len3padd
                    # If the addition is a SNV w.r.t precursor, then it will be => iso_add3p:N. Number of non-template nucleotides added at 3p. 
                    #del sub[-len3padd:]
                    #del master_seq_bc[-len3padd:]
                    #print("3p_add:" + a3p + "Len"+ str(len3padd))
                    #cigar3padd = str(len3padd)+"I"
                if iso_3p_del != "":
                    d3p = "".join(iso_3p_del[::-1])
                    len3pdel = len(d3p)
                    variant+= "iso_3p:-"+str(len3pdel)+","
                    end = end - len3pdel
                    gen_end = gen_end - len3pdel
                    if gen_strand == "-":
                        gen_start+=len3pdel
                    #del sub[-len3pdel:]
                    #del master_seq_bc[-len3pdel:]
                    #print("3p_del:" + d3p +"Len"+ str(len3pdel))
                    #cigar3pdel = str(len3pdel)+"D"
                # Now, these array's for reference and query doesn't have changes at the 5' or 3' ends. So, any variant correspond to internal changes
                #print(precursorSeq[start-1:end])
                new_var_type={}
                if iso_sub:
                    for xs in iso_sub.keys():
                        if xs == 7:
                            new_var_type["iso_snv_central_offset,"] = "1"
                        elif xs >= 1 and xs <= 6:
                            new_var_type["iso_snv_seed,"] = "1"
                        elif xs >= 8 and xs <= 12:
                            new_var_type["iso_snv_central,"] = "1"
                        elif xs >= 13 and xs <= 17:
                            new_var_type["iso_snv_central_supp,"] = "1"
                        else:
                            new_var_type["iso_snv,"] = "1"
                    variant+= "".join(new_var_type.keys())
                    #print(new_var_type)
                if variant.endswith(","):
                    variant = re.sub(',$','',variant)
                #print(variant)
                """
                # PREPARING CIGAR BODY 
                """
                #print("Arun:")
                #print(master_seq_bc)
                #print("Patil:")
                match_case = ""
                for snv_id, snv in enumerate(master_seq_bc):
                    if snv == sub[snv_id]:
                        match_case+= "M"
                    elif snv == "-":
                        match_case+= "M"
                        #match_case+= "I"
                    elif sub[snv_id] == "_":
                        #match_case+= "D"
                        match_case+= "M"
                    else: 
                        match_case+=master_seq_bc[snv_id] # 11MA7M to indicates there is a mismatch at position 12, where A is the reference nucleotide.
                        #match_case+=sub[snv_id]
                ## CREATING THE CIGAR FORMAT HERE ##
                match_case = match_case.replace("+","")
                #print(seq_master, match_case, seq_m)
                count_4cigar=0
                iso_cigar="" # This varialbe is actually CIGAR variable which collects CIGAR information
                for isx, ist in enumerate(match_case):
                    if isx != 0:
                        if ist == match_case[isx-1]:
                            count_4cigar +=1
                        else:
                            if count_4cigar != 1:
                                iso_cigar += str(count_4cigar)+match_case[isx-1]
                                count_4cigar =1
                            else: 
                                iso_cigar += match_case[isx-1]
                                count_4cigar =1
                    else:
                        count_4cigar +=1
                if count_4cigar != 1:
                    iso_cigar += str(count_4cigar)+ist
                else: 
                    iso_cigar += ist

                if "A" not in match_case and "T" not in match_case and "G" not in match_case and "C" not in match_case:
                    iso_cigar = str(len(seq_m))+"M"
                else:
                    pass
                    #print(seq_m, iso_cigar)
                if variant == "":
                    variant = "iso_snv"
                iso_mi_var = seq_master+"\t"+version_db+"\t"+type_rna+"\t"+str(start)+"\t"+str(end)+"\t.\t+\t.\tRead="+seq_m+"; UID="+uid_val+"; Name="+ seq_master +"; Parent="+req_precursor_name+"; Variant="+variant+"; Cigar="+iso_cigar+"; Expression="+canonical_expression +"; Filter=Pass; Hits="+ canonical_expression + "\n"
                gffwrite.write(iso_mi_var)
                iovariant = re.sub(',',';',variant)
                iovarlist = iovariant.split(";")
                for iv in iovarlist:
                    if iv != "":
                        if ":" in iv:
                            iv = iv.split(":")[0]
                        try:
                            JS_variantType_dataDict[str(iv)] += int(JS_exprn_total)
                        except KeyError:
                            JS_variantType_dataDict[str(iv)] = int(JS_exprn_total)

                #forJSgffData.write(str(seq_m)+","+seq_master+","+iovariant+","+str(canonical_expression)+","+JS_exprn_total+"\n")
                valStr = ','.join([str(elem) for elem in iovarlist]) +"#"+str(canonical_expression)
                try:
                    JS_hmap_iso_miRVar[seq_master].append(valStr)
                except KeyError:
                    JS_hmap_iso_miRVar[seq_master] = [valStr]

                try:
                    JS_hmap_iso_miRcounts[seq_master].append(canonical_expression.split(","))
                except KeyError:
                    JS_hmap_iso_miRcounts[seq_master] = [canonical_expression.split(",")]
                """
                FOR SAM/BAM FILE FORMAT: gen_chr, gen_start, gen_end
                """
                # bam_can_dict[seq_m] = str(gen_chr) +"\t"+ str(gen_start) +"\t"+ iso_cigar + "\t"+ gen_strand
                #print(seq_m+"\t"+ str(gen_chr) +"\t"+ str(gen_start) +"\t"+ iso_cigar+"\n")

                #print("--**END**--")
        except KeyError:
            pass
            #print(seq_m+"\t"+seq_master)
            #ACTGGCCTTGGAGTCAGAAGGC  hsa-miR-378g
    html_data.openDoChartJSD()
    for xy, xz in JS_variantType_dataDict.items():
        html_data.donutChartJSD(xy, xz)
    html_data.closeDoChartJSD()
    #print(JS_hmap_list_ref)
    sum_JS_hmap_iso_miRcounts=dict()
    #print(JS_hmap_iso_miRcounts)
    for ji, jj in JS_hmap_iso_miRcounts.items():
        new_list = [list(map(int, lst)) for lst in jj]
        sum_JS_hmap_iso_miRcounts[ji] = [sum(i) for i in zip(*new_list)]

    #print(sum_JS_hmap_iso_miRcounts)
    #print(JS_hmap_iso_miRVar)
    df_rpm = mirRPM_completeSet.reset_index()
    for indx, name in enumerate(base_names):
        # FIRST PICK TOP 40 isomiRs FROM EACH SAMPLE
        tempDict = dict()
        for td1, td2 in sum_JS_hmap_iso_miRcounts.items():
            tempDict[td1] = td2[indx]

        req_mir_names = sorted(df_rpm[['miRNA', name]].values.tolist(), key=lambda x: x[1], reverse=True)
        #req_mir_names = sorted(df_rpm[['miRNA', name]].values.tolist(), key=lambda x: x[1], reverse=True)[:20]
        only_AbundantmiRs = [ re.sub('/.*', '', item[0]) for item in req_mir_names] 
        abundant_miRsJS = []
        for chkExts in only_AbundantmiRs:
            if len(abundant_miRsJS) == 20:
                break
            try:
                if JS_hmap_iso_miRVar[chkExts]:
                    abundant_miRsJS.append(chkExts)
            except KeyError:
                pass
        
        req_mir_names = sorted(abundant_miRsJS)
        #req_mir_names = sorted(only_AbundantmiRs)
        #req_mir_names = sorted(tempDict, key=tempDict.get, reverse=True)[:20] # returns keys for top 40 miRNA expression (Only keys)
        html_data.openisoHmapTop(name, req_mir_names)
        for kindx, km in enumerate(req_mir_names):
            #print(km, JS_hmap_iso_miRVar[km])
            vary = dict()
            var_list_items = ['iso_3p:-1', 'iso_3p:-2', 'iso_5p:+2','iso_5p:+1','iso_3p:+1','iso_add3p:+1','iso_3p:+2','iso_add3p:+2','iso_3p:+3', 'iso_add3p:+3', 'iso_3p:+4', 'iso_add3p:+4', 'iso_5p:-1', 'iso_add5p:+1', 'iso_5p:-2', 'iso_add5p:+2', 'iso_5p:-3', 'iso_add5p:+3', 'iso_5p:-4', 'iso_add5p:+4', 'iso_snv_seed', 'iso_snv_central_offset', 'iso_snv_central', 'iso_snv_central_supp','iso_snv'] 
            for k_var in var_list_items: # Initializing dictionaries to zero to avoid try catch except block in the later
                vary[k_var] = 0

            for valsy in JS_hmap_iso_miRVar[km]:
                lhs = valsy.split("#")
                #print(valsy, lhs)
                reqExprval = int(lhs[1].split(",")[indx])
                reqVaritems = lhs[0].split(",")
                for varkeys in reqVaritems:
                    #print(km, varkeys, reqExprval, valsy)
                    try:
                        vary[varkeys]+= int(reqExprval)
                    except KeyError:
                        vary[varkeys] = int(reqExprval)

            iso_3pm1 = "%.1f" %math.log2(vary['iso_3p:-1']) if vary['iso_3p:-1'] > 0 else 0
            iso_3pm2 = "%.1f" %math.log2(vary['iso_3p:-2']) if vary['iso_3p:-2'] > 0 else 0
            iso_5pm1 = "%.1f" %math.log2(vary['iso_5p:+1']) if vary['iso_5p:+1'] > 0 else 0
            iso_5pm2 = "%.1f" %math.log2(vary['iso_5p:+2']) if vary['iso_5p:+2'] > 0 else 0

            iso_3pp1 = "%.1f" %math.log2(vary['iso_3p:+1'] + vary['iso_add3p:+1']) if (vary['iso_3p:+1'] + vary['iso_add3p:+1']) > 0 else 0
            iso_3pp2 = "%.1f" %math.log2(vary['iso_3p:+2'] + vary['iso_add3p:+2']) if (vary['iso_3p:+2'] + vary['iso_add3p:+2']) > 0 else 0
            iso_3pp3 = "%.1f" %math.log2(vary['iso_3p:+3'] + vary['iso_add3p:+3']) if (vary['iso_3p:+3'] + vary['iso_add3p:+3']) > 0 else 0
            iso_3pp4 = "%.1f" %math.log2(vary['iso_3p:+4'] + vary['iso_add3p:+4']) if (vary['iso_3p:+4'] + vary['iso_add3p:+4']) > 0 else 0
            iso_5pp1 = "%.1f" %math.log2(vary['iso_5p:-1'] + vary['iso_add5p:+1']) if (vary['iso_5p:-1'] + vary['iso_add5p:+1']) > 0 else 0
            iso_5pp2 = "%.1f" %math.log2(vary['iso_5p:-2'] + vary['iso_add5p:+2']) if (vary['iso_5p:-2'] + vary['iso_add5p:+2']) > 0 else 0
            iso_5pp3 = "%.1f" %math.log2(vary['iso_5p:-3'] + vary['iso_add5p:+3']) if (vary['iso_5p:-3'] + vary['iso_add5p:+3']) > 0 else 0
            iso_5pp4 = "%.1f" %math.log2(vary['iso_5p:-4'] + vary['iso_add5p:+4']) if (vary['iso_5p:-4'] + vary['iso_add5p:+4']) > 0 else 0
            
            iso_3pp5 = 0
            iso_5pp5 = 0
            for item_3p in [x for x in vary.keys() if "iso_3p:+" in x]:
                if 'iso_3p:+1' not in item_3p and 'iso_3p:+2' not in item_3p and 'iso_3p:+3' not in item_3p and 'iso_3p:+4' not in item_3p:
                    iso_3pp5 += int(vary[item_3p])
            for item_add3p in [ xa for xa in vary.keys() if "iso_add3p:+" in xa]:
                if 'iso_add3p:+1' not in item_add3p and 'iso_add3p:+2' not in item_add3p and 'iso_add3p:+3' not in item_add3p and 'iso_add3p:+4' not in item_add3p:
                    iso_3pp5 += int(vary[item_add3p])
            for item_5p in [x for x in vary.keys() if "iso_5p:-" in x]:
                if 'iso_5p:-1' not in item_5p and 'iso_5p:-2' not in item_5p and 'iso_5p:-3' not in item_5p and 'iso_5p:-4' not in item_5p:
                    iso_5pp5 += int(vary[item_5p])
            for item_add5p in [ xa for xa in vary.keys() if "iso_add5p:+" in xa]:
                if 'iso_add5p:+1' not in item_add5p and 'iso_add5p:+2' not in item_add5p and 'iso_add5p:+3' not in item_add5p and 'iso_add5p:+4' not in item_add5p:
                    iso_5pp5 += int(vary[item_add5p])

            iso_3pp5 = "%.1f" %math.log2(iso_3pp5) if iso_3pp5 > 0 else 0
            iso_5pp5 = "%.1f" %math.log2(iso_5pp5) if iso_5pp5 > 0 else 0
            try:
                ref_val = "%.1f" %math.log2(int(list(JS_hmap_list_ref[km])[indx])) if int(list(JS_hmap_list_ref[km])[indx]) > 0 else 0
            except KeyError:
                ref_val = 0 
            snv_1 = "%.1f" %math.log2(vary['iso_snv_seed']) if vary['iso_snv_seed'] > 0 else 0
            snv_2 = "%.1f" %math.log2(vary['iso_snv_central_offset']) if vary['iso_snv_central_offset'] > 0 else 0
            snv_3 = "%.1f" %math.log2(vary['iso_snv_central']) if vary['iso_snv_central'] > 0 else 0
            snv_4 = "%.1f" %math.log2(vary['iso_snv_central_supp']) if vary['iso_snv_central_supp'] > 0 else 0
            snv_5 = "%.1f" %math.log2(vary['iso_snv']) if vary['iso_snv'] > 0 else 0

            iso_data_js = "["+ str(kindx) + ",0," + str(iso_3pp5) + "],[" + str(kindx) + ",1," + str(iso_3pp4) + "],[" + str(kindx) + ",2," + str(iso_3pp3) + "],[" + str(kindx) + ",3," + str(iso_3pp2) + "],[" + str(kindx) + ",4," + str(iso_3pp1) + "],[" + str(kindx) + ",5," + str(ref_val) + "],[" + str(kindx) + ",6," + str(iso_3pm1) + "],[" + str(kindx) + ",7," + str(iso_3pm2) + "],[" + str(kindx) + ",8," + str(snv_1) + "],[" + str(kindx) + ",9," + str(snv_2) + "],[" + str(kindx) + ",10," + str(snv_3) + "],[" + str(kindx) + ",11," + str(snv_4) + "],[" + str(kindx) + ",12," + str(snv_5) + "],[" + str(kindx) + ",13," + str(iso_5pm2) + "],[" + str(kindx) + ",14," + str(iso_5pm1) + "],[" + str(kindx) + ",15," + str(iso_5pp1) + "],[" + str(kindx) + ",16," + str(iso_5pp2) + "],[" + str(kindx) + ",17," + str(iso_5pp3) + "],[" + str(kindx) + ",18," + str(iso_5pp4) + "],[" + str(kindx) + ",19," + str(iso_5pp5) + "]"
            html_data.isoHmapData(iso_data_js)
            #print(km, vary, iso_3pp1)
            #print(iso_data_js)
        html_data.closeisoHmapBottom()

    """
    if args.bam_out:
        mirna_samFile = Path(workDir)/"miRge3_miRNA.sam"
        genC=0
        genS=0
        cig=0
        with open(mirna_samFile) as miSam:
            for mi_idx, mi_sam in enumerate(miSam):
                mi_sam = mi_sam.strip()
                mi_sam_list = mi_sam.split("\t")
                try:
                #if bam_can_dict[mi_sam_list[0]]:
                    sam_exprn_list = bam_expression_dict[mi_sam_list[0]]
                    for ex_idx, exprn in enumerate(sam_exprn_list):
                        (genC, genS, cig, strand) = bam_can_dict[mi_sam_list[0]].split("\t")
                        if exprn >= 1:
                            file_sam_name = str(base_names[ex_idx]) +".sam"
                            sam_name = Path(workDir)/file_sam_name
                            xbam = open(sam_name, "a+")
                            for numexp in range(exprn):
                                #print()
                                #readname = "r"+str(mi_idx) + "_" + str(numexp)
                                readname = mi_sam_list[0] + "_" + str(numexp)
                                phredQual = "I"*len(mi_sam_list[0])
                                cigar = str(len(mi_sam_list[0]))+"M"
                                #xbamout = readname+"\t"+mi_sam_list[1]+"\t"+genC+"\t"+str(genS)+"\t"+mi_sam_list[4]+"\t"+mi_sam_list[5]+"\t"+mi_sam_list[6]+"\t"+mi_sam_list[7]+"\t"+mi_sam_list[8]+"\t"+mi_sam_list[9]+"\t"+mi_sam_list[10]+"\n"
                                if strand == "+":
                                    xbamout = readname+"\t"+mi_sam_list[1]+"\t"+genC+"\t"+str(genS)+"\t"+mi_sam_list[4]+"\t"+cigar+"\t"+mi_sam_list[6]+"\t"+mi_sam_list[7]+"\t"+mi_sam_list[8]+"\t"+mi_sam_list[0]+"\t"+phredQual+"\n"
                                else:
                                    xbamout = readname+"\t"+mi_sam_list[1]+"\t"+genC+"\t"+str(genS)+"\t"+mi_sam_list[4]+"\t"+cigar+"\t"+mi_sam_list[6]+"\t"+mi_sam_list[7]+"\t"+mi_sam_list[8]+"\t"+mi_sam_list[0][::-1]+"\t"+phredQual+"\n"

                                xbam.write(xbamout)
                            xbam.close()
                except KeyError:
                    pass
    """

def addDashNew(seq, totalLength, start, end):
    newSeq = '-'*(start-1)+seq+'-'*(totalLength-end)
    return newSeq

def trfTypes(seq, tRNAName, start, trnaStruDic):
    if 'pre_' not in tRNAName:
        tRNASeq = trnaStruDic[tRNAName]['seq']
        tRNAStru = trnaStruDic[tRNAName]['stru']
        tRNASeqLen = len(tRNASeq)
        # anticodonStart and anticodonEnd is 1-based position, so change it into 0-based
        anticodonStart = trnaStruDic[tRNAName]['anticodonStart']-1
        anticodonEnd = trnaStruDic[tRNAName]['anticodonEnd']-1
        if start == 0:
            if start+len(seq) == tRNASeqLen:
                trfType = 'tRF-whole'
            elif start+len(seq)-1 >= anticodonStart-2 and start+len(seq)-1 <= anticodonStart+1:
                trfType = "5'-half"
            else:
                trfType = "5'-tRF"
        else:
            if start+len(seq)-1 >= tRNASeqLen-1-2 and start+len(seq)-1 <= tRNASeqLen-1:
                if start >= anticodonStart-1 and start <= anticodonStart+2:
                    trfType = "3'-half"
                else:
                    trfType = "3'-tRF"
            else:
                trfType = 'i-tRF'
    else:
        trfType = 'tRF-1'
    return trfType


def summarize(args, workDir, ref_db,base_names, pdMapped, sampleReadCounts, trimmedReadCounts, trimmedReadCountsUnique):
    """
    THIS FUNCTION IS CALLED FIRST FROM THE miRge3.0 to summarize the output.  
    """
    global html_data
    html_data = FormatJS(workDir) 
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
    try:
        with open(mergeFile, "r") as merge_file:
            for line in merge_file:
                line_content = line.strip().split(',')
                for item in line_content[1:]:
                    mirMergedNameDic.update({item:line_content[0]})
                    mirMergedDataframeDic.update({line_content[0]:"1"})
    except FileNotFoundError:
        pass
    
    #allSequences = pdMapped.index.shape[0]
    #print(allSequences)

    pdMapped = pdMapped.reset_index(level=['Sequence'])
    subpdMapped = pdMapped[(pdMapped['exact miRNA'].astype(bool) | pdMapped['isomiR miRNA'].astype(bool))]
    cannonical = pdMapped[pdMapped['exact miRNA'].astype(bool)]
    isomirs = pdMapped[pdMapped['isomiR miRNA'].astype(bool)]
    cannonical_4bam = pdMapped[pdMapped['exact miRNA'].astype(bool)]
    isomirs_4bam = pdMapped[pdMapped['isomiR miRNA'].astype(bool)]
    cannonical_4ie = cannonical
    isomirs_4ie = isomirs
    cannonical_4gff = cannonical
    isomirs_4gff = isomirs
    #### MOVED TWO IF CONDITIONS args.gff_out or args.bam_out: and args.bam_out: from next line 
    
    if args.spikeIn:
        cannonical = cannonical.drop(columns=['Sequence','hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','annotFlag','isomiR miRNA','spike-in'])
        isomirs = isomirs.drop(columns=['Sequence','exact miRNA','hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','annotFlag','spike-in'])
        cannonical_4ie = cannonical_4ie.drop(columns=['hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','annotFlag','isomiR miRNA','spike-in'])
        isomirs_4ie = isomirs_4ie.drop(columns=['exact miRNA','hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','annotFlag','spike-in'])
        subpdMapped = subpdMapped.drop(columns=['Sequence','hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','annotFlag','spike-in'])
    else:
        cannonical = cannonical.drop(columns=['Sequence','hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','annotFlag','isomiR miRNA'])
        isomirs = isomirs.drop(columns=['Sequence','exact miRNA','hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','annotFlag'])
        cannonical_4ie = cannonical_4ie.drop(columns=['hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','annotFlag','isomiR miRNA'])
        isomirs_4ie = isomirs_4ie.drop(columns=['exact miRNA','hairpin miRNA','mature tRNA','primary tRNA','snoRNA','rRNA','ncrna others','mRNA','annotFlag'])
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
    bwtCommand = Path(args.bowtie_path)/"bowtie-inspect" if args.bowtie_path else "bowtie-inspect"
    bwtExec = str(bwtCommand) + " -n " + str(indexFiles)
    #bwtExec = "bowtie-inspect -n /home/arun/repositories/Project_120919/mirge/Libs/human/index.Libs/human_mirna_miRBase"
    bowtie = subprocess.run(str(bwtExec), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
    if bowtie.returncode==0:
        bwtOut = bowtie.stdout
        bwtErr = bowtie.stderr
        lines = bwtOut.strip()
        for srow in lines.split('\n'):
            if "segs:" in srow:
                srow = srow.split(" ")[0]
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

    if args.gff_out:
    #if args.gff_out or args.bam_out:
        pre_mirDict = dict()
        mirDict = dict()
        filenamegff = workDir/"sample_miRge3.gff"
        maturefname = args.organism_name + "_mature_" + ref_db + ".fa"
        pre_fname = args.organism_name + "_hairpin_" + ref_db
        fasta_file = Path(args.libraries_path)/args.organism_name/"fasta.Libs"/maturefname
        precursor_file = Path(args.libraries_path)/args.organism_name/"index.Libs"/pre_fname
        annotation_pre_fname = args.organism_name+"_"+ref_db+".gff3"
        annotation_lib = Path(args.libraries_path)/args.organism_name/"annotation.Libs"/annotation_pre_fname

        bwtCommand = Path(args.bowtie_path)/"bowtie-inspect" if args.bowtie_path else "bowtie-inspect"
        bwtExec = str(bwtCommand) + " -a 20000 -e "+ str(precursor_file)
        bowtie = subprocess.run(str(bwtExec), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
        #READING PRECURSOR miRNA SEQUENCES INFORMATION IN A DICTIONARY (pre_mirDict)
        if bowtie.returncode==0:
            bwtOut = bowtie.stdout
            bwtErr = bowtie.stderr
        for srow in bwtOut.split('\n'):
            if '>' in srow:
                srow = srow.replace(">","")
                headmil = srow.split(" ")[0]
                #if "MirGeneDB" in ref_db:
                #    headmil = headmil.split("_")[0]
            else:
                #headmil = '-'.join(headmil.split('-')[:-1])
                pre_mirDict[headmil] = srow
        #READING MATURE miRNA SEQUENCES INFORMATION IN A DICTIONARY (mirDict)
        with open(fasta_file) as mir:
            for mil in mir:
                mil = mil.strip()
                if '>' in mil:
                    headmil_mi = mil.replace(">","")
                    #if "MirGeneDB" in ref_db:
                    #    headmil_mi = headmil_mi.split("_")[0]
                else:
                    mirDict[headmil_mi] = mil
        d = Differ()
        create_gff(args, pre_mirDict, mirDict, d, filenamegff, cannonical_4gff, isomirs_4gff, base_names, ref_db, annotation_lib, workDir, mirRPM_completeSet)

    if args.bam_out:
    #if args.bam_out:
        header = sam_header(args)
        for names in base_names:
            file_sam_nameH = str(names) +".sam"
            sam_nameH = Path(workDir)/file_sam_nameH
            with open(sam_nameH,"w+") as samH:
                #samH.write("@HD\tVN:3.0\tSO:coordinate\n")
                samH.write(header)
        cols1 = ["Sequence","exact miRNA"] + base_names
        cols2 = ["Sequence","isomiR miRNA"] + base_names
        canonical_gff = pd.DataFrame(cannonical_4bam, columns= cols1).values.tolist() # Gives list of list containg Sequence, miRNA name, expression values for the samples - ref miRNA
        isomir_gff = pd.DataFrame(isomirs_4bam, columns= cols2).values.tolist() # Gives list of list containg Sequence, miRNA name, expression values for the samples - isomiR
        canonical_gff.extend(isomir_gff)  # APPENDING THE LIST OF ISOMIRS TO CANONICAL # Making one big list to get coordinates and anntations for GFF3 format of miRTop

        pd_frame = ['snoRNA','rRNA','ncrna others','mRNA']
        bwt_idx_prefname = ['snorna','rrna','ncrna_others','mrna']
        for igv_idx, igv_name in enumerate(pd_frame):
            dfRNA2sam = pdMapped[pdMapped[igv_name].astype(bool)]
            pre_cols_birth = ["Sequence", igv_name]
            cols1 = pre_cols_birth + base_names
            df_sam_out = pd.DataFrame(dfRNA2sam, columns= cols1) # Gives list of list containg Sequence, RNA type, expression values for the samples 
            df_expr_list = df_sam_out.values.tolist()
            rna_type = args.organism_name + "_" + bwt_idx_prefname[igv_idx]
            index_file_name = Path(args.libraries_path)/args.organism_name/"index.Libs"/rna_type
            bow2bam(args, workDir, ref_db, df_expr_list, base_names, index_file_name, rna_type, bwt_idx_prefname[igv_idx])
        #https://stackoverflow.com/questions/35125062/how-do-i-join-2-columns-of-a-pandas-data-frame-by-a-comma
        rna_type = args.organism_name +"_mirna_"+ref_db
        index_file_name = Path(args.libraries_path)/args.organism_name/"index.Libs"/rna_type
        bow2bam(args, workDir, ref_db, canonical_gff, base_names, index_file_name, rna_type, "miRNA" )
        #bow2bam(args, workDir, ref_db, canonical_gff, base_names, index_file_name, rna_type, "mirna" )
        # Create BAM for hairpin miRNA 
        rna_type = args.organism_name +"_hairpin_"+ref_db
        index_file_name = Path(args.libraries_path)/args.organism_name/"index.Libs"/rna_type
        hpin_exprn = pdMapped[pdMapped["hairpin miRNA"].astype(bool)]
        cols_hpin = ["Sequence","hairpin miRNA"] + base_names
        hpin_exlist = pd.DataFrame(hpin_exprn, columns=cols_hpin).values.tolist()
        bow2bam(args, workDir, ref_db, hpin_exlist, base_names, index_file_name, rna_type, "hairpin_miRNA" )
        # hairpin_miRNA
        createBAM(args, workDir, base_names)

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
    
    """
    Calcuate isomir entropy
    """
    def calcEntropy(inputList):
        sum1 = sum(inputList)
        entropy = 0
        for i in range(len(inputList)):
            if inputList[i] > 1:
                freq = float(inputList[i])/sum1
                entropy = entropy + -1*freq*math.log(freq, 2)
        return entropy

    def create_ie(args, cannonical, isomirs, base_names, workDir, Filtered_miRNA_Reads):
        isomirFile = Path(workDir)/"isomirs.csv"
        isomirSampleFile = Path(workDir)/"isomirs.samples.csv"
        outf1 = open(isomirFile, 'w')
        outf2 = open(isomirSampleFile, 'w')
        outf1.write('miRNA,sequence')
        outf2.write('miRNA')
        for i in range(len(base_names)):
            outf1.write(','+base_names[i])
            outf2.write(','+base_names[i]+' isomir+miRNA Entropy')
            outf2.write(','+base_names[i]+' Canonical Sequence')
            outf2.write(','+base_names[i]+' Canonical RPM')
            outf2.write(','+base_names[i]+' Top Isomir RPM')
        outf1.write(',Entropy\n')
        outf2.write('\n')
        pre_cols1 = ["Sequence","exact miRNA"] 
        pre_cols2 = ["Sequence","isomiR miRNA"]
        cols1 = pre_cols1 + base_names
        cols2 = pre_cols2 + base_names
        can_gff_df = pd.DataFrame(cannonical, columns= cols1) # Gives list of list containg Sequence, miRNA name, expression values for the samples - ref miRNA
        iso_gff_df = pd.DataFrame(isomirs, columns= cols2) # Gives list of list containg Sequence, miRNA name, expression values for the samples - isomiR 
        canonical_gff = can_gff_df.values.tolist() 
        isomir_gff = iso_gff_df.values.tolist()
        freq_list=[]
        for fname in base_names:
            try:
                freq_list.append(1000000/Filtered_miRNA_Reads[fname])
            except ZeroDivisionError:
                freq_list.append(0)
        maxEntropy = math.log(len(base_names), 2)

        """
        Collecting miRNA values across each samples into an array
        """
        miR_can = {}
        for each_can in canonical_gff:
            canValScore = each_can[2:]
            if ".SNP" in each_can[1]:
                each_can[1] = each_can[1].split('.')[0]
            try:
                miR_can[each_can[1]].append(canValScore)
            except KeyError:
                miR_can[each_can[1]] = [canValScore]

        """
        Collecting miRNA values across each samples into an array - Here the values for each sample is summed 
        """
        for key_mir, val_mir in miR_can.items():
            res = [sum(i) for i in zip(*val_mir)]
            miR_can[key_mir] = res
        
        miR_iso = {}
        for each_isoSeq in isomir_gff:
            valueScore = each_isoSeq[2:]
            entropy = calcEntropy(valueScore)
            if maxEntropy == 0:
                entropy= "NA"
            else:
                entropy = str(entropy/maxEntropy)
            emptyListEntropy = []
            #topIsomir = []
            #isomirSum = []
            for idxn, ival in enumerate(valueScore):
                emptyListEntropy.append(str(ival*freq_list[idxn]))
                #topIsomir.append(str(max(valueScore)*rpmFactor))
                #isomirSum.append(str(sum(sampleIsomirs[sampleLane])*rpmFactor))
            samplesEntropy = "\t".join(emptyListEntropy)
            if ".SNP" in each_isoSeq[1]:
                each_isoSeq[1] = each_isoSeq[1].split('.')[0]
            try:
                miR_iso[each_isoSeq[1]].append(valueScore)
            except KeyError:
                miR_iso[each_isoSeq[1]] = [valueScore]
            outf1.write(each_isoSeq[1]+"\t"+each_isoSeq[0]+"\t"+ samplesEntropy +"\t"+ entropy + "\n")

        for isokey, isoval in miR_iso.items():
            #print(list(zip(*isoval)))
            res = [i for i in zip(*isoval)]
            isomirOut = [isokey]
            for xn, x in enumerate(res):
                iso_vals_asList =list(x)
                topIsomir = max(iso_vals_asList)*freq_list[xn]
                isomirSum = sum(iso_vals_asList)*freq_list[xn]
                if isokey in miR_can:
                    iso_can_vals_list = iso_vals_asList + [miR_can[isokey][xn]]
                    miRNARPM = miR_can[isokey][xn] * freq_list[xn] 
                    sampleEntropyWithmiRNA = calcEntropy(iso_can_vals_list)
                    maxEntropy = len(iso_vals_asList)
                    if maxEntropy > 1:
                        sampleEntropyWithmiRNA = str(sampleEntropyWithmiRNA/(math.log(maxEntropy,2)))
                    else:
                        sampleEntropyWithmiRNA = 'NA'
                    isomirOut.append(sampleEntropyWithmiRNA)
                    combined = miRNARPM + isomirSum
                    if combined >0:
                        isomirOut.append(str(100.0*miRNARPM/combined))
                    else:
                        isomirOut.append('NA')
                    isomirOut.append(str(miRNARPM))
                    isomirOut.append(str(topIsomir))
                    
                else:
                    pass
                    #print(list(x))
            if len(isomirOut) > 1:
                outf2.write(','.join(isomirOut))
                outf2.write('\n')

        outf1.close()
        outf2.close()
            #print(isomirOut)
            #print(isokey, isoval)
            #print(res)
            #print(each_isoSeq)
            # ['AAAAAACTCTAAACAA', 'hsa-miR-3145-5p', 0, 1]

    if args.isoform_entropy:        
        create_ie(args, cannonical_4ie, isomirs_4ie, base_names, workDir, Filtered_miRNA_Reads)
    
    if args.AtoI:
        reqCols = ['miRNA']+base_names
        mirCounts_completeSet = mirCounts_completeSet.reset_index(level=['miRNA'])
        mirCC = pd.DataFrame(mirCounts_completeSet, columns= reqCols).values.tolist() 
        mirDic={}
        for mC in mirCC:
            mirDic[mC[0]] = mC[1:]
        pre_cols1 = ["Sequence"] 
        cols1 = pre_cols1 + base_names
        cols2 = pre_cols1 + base_names
        can_ai_df = pd.DataFrame(cannonical_4ie, columns= cols1) # Gives list of list containg Sequence, miRNA name, expression values for the samples - ref miRNA
        iso_ai_df = pd.DataFrame(isomirs_4ie, columns= cols2) # Gives list of list containg Sequence, miRNA name, expression values for the samples - isomiR 
        canonical_ai = can_ai_df.values.tolist() 
        onlyCannon = canonical_ai
        onlyCanmiRNA={}
        for oC in onlyCannon:
            onlyCanmiRNA[oC[0]] = oC[1:]
        isomir_ai = iso_ai_df.values.tolist()
        canonical_ai.extend(isomir_ai) 
        seqDic={}
        for sD in canonical_ai:
            seqDic[sD[0]] = sD[1:]
        #print(seqDic)
        a2i_editing(args, cannonical_4ie, isomirs_4ie, base_names, workDir, Filtered_miRNA_Reads, mirMergedNameDic, mirDic, ref_db, seqDic, onlyCanmiRNA)
        pass
    
    if args.tRNA_frag:
        m_trna_pre = pdMapped[pdMapped['mature tRNA'].astype(bool)]
        p_trna_pre = pdMapped[pdMapped['primary tRNA'].astype(bool)]
        m_trna_cols1 = ["Sequence","mature tRNA"] + base_names
        p_trna_cols2 = ["Sequence","primary tRNA"] + base_names
        m_trna = pd.DataFrame(m_trna_pre, columns= m_trna_cols1).values.tolist() # Gives list of list containg Sequence, mature tRNA, expression values for the samples - mature tRNA 
        p_trna = pd.DataFrame(p_trna_pre, columns= p_trna_cols2).values.tolist() # Gives list of list containg Sequence, primary tRNA, expression values for the samples - primary tRNA
        trnaStruDic={}
        fname = args.organism_name+'_trna.str'
        trna_stru_file = Path(args.libraries_path)/args.organism_name/"annotation.Libs"/fname
        allFiles_inPlace = 1
        try:
            with open(trna_stru_file, 'r') as inf:
                a=0
                lines = inf.readlines()
                for idx, xline in enumerate(lines):
                    xline=xline.strip()
                    if xline.startswith(">"):
                        trnaName = xline.replace(">","")
                        a+=1
                    elif a == 1:
                        a+=1
                        trnaSeq = xline
                    elif a == 2:
                        a=0
                        trnaStru = xline
                        anticodonStart = trnaStru.index('XXX')+1
                        anticodonEnd = anticodonStart+2
                        trnaStruDic.update({trnaName:{'seq':trnaSeq, 'stru':trnaStru, 'anticodonStart':anticodonStart, 'anticodonEnd':anticodonEnd}})
        except IOError:
            allFiles_inPlace = 0
            print(f"File {trna_stru_file} does not exist!!\nProceeding the annotation with out -trf\n")

        fname2 = args.organism_name+'_trna_aminoacid_anticodon.csv'
        trna_aa_anticodon_file = Path(args.libraries_path)/args.organism_name/"annotation.Libs"/fname2
        trnaAAanticodonDic = {}
        try:
            with open(trna_aa_anticodon_file, 'r') as inf:
                for line in inf:
                    contentTmp = line.strip().split(',')
                    trnaAAanticodonDic.update({contentTmp[0]:{'aaType':contentTmp[1], 'anticodon':contentTmp[2]}})
        except IOError:
            allFiles_inPlace = 0
            print(f"File {trna_aa_anticodon_file} does not exist!!\nProceeding the annotation with out -trf\n")
        
        fname3 = args.organism_name+'_trna_deduplicated_list.csv'
        trna_duplicated_list_file = Path(args.libraries_path)/args.organism_name/"annotation.Libs"/fname3
        duptRNA2UniqueDic = {}
        try:
            with open(trna_duplicated_list_file, 'r') as inf:
                line = inf.readline()
                line = inf.readline()
                while line != '':
                    contentTmp = line.strip().split(',')
                    for item in contentTmp[1].split('/'):
                        duptRNA2UniqueDic.update({item.strip():contentTmp[0].strip()})
                    line = inf.readline()
        except IOError:
            allFiles_inPlace = 0
            print(f"File {trna_duplicated_list_file} does not exist!!\nProceeding the annotation with out -trf\n")
        
        fname4 = args.organism_name+'_tRF_infor.csv'
        tRF_infor_file = Path(args.libraries_path)/args.organism_name/"annotation.Libs"/fname4
        tRNAtrfDic = {}
        try:
            with open(tRF_infor_file, 'r') as  inf:
                line = inf.readline()
                line = inf.readline()
                while line != '':
                    content = line.strip().split(',')
                    tRNAName = content[0].split('_Cluster')[0]
                    tRNAClusterName = content[0]
                    seq = content[4]
                    tRNAlength = len(content[5])
                    start = int(content[3].split('-')[0])
                    end = int(content[3].split('-')[1])
                    if tRNAName not in tRNAtrfDic.keys():
                        tRNAtrfDic.update({tRNAName:{}})
                    tRNAtrfDic[tRNAName].update({addDashNew(seq, tRNAlength, start, end):tRNAClusterName})
                    line = inf.readline()
        except IOError:
            allFiles_inPlace = 0
            print(f"File {tRF_infor_file} does not exist!!\nProceeding the annotation with out -trf\n")

        # Load predifined tRF merged file
        fname5 = args.organism_name+"_tRF_merges.csv"
        tRF_merge_file = Path(args.libraries_path)/args.organism_name/"annotation.Libs"/fname5
        trfMergedNameDic = {}
        trfMergedList = []
        try:
            with open(tRF_merge_file, 'r') as inf:
                for line in inf:
                    tmp = line.strip().split(',')
                    mergedName = tmp[0]
                    trfMergedList.append(mergedName)
                    for item in tmp[1].split('/'):
                        trfMergedNameDic.update({item:mergedName})
        except IOError: 
            allFiles_inPlace = 0
            print(f"File {tRF_merge_file} does not exist!!\nProceeding the annotation with out -trf\n")
        
        pretrnaNameSeqDic = {}
        file_pre_tRNA = args.organism_name+'_pre_trna'
        indexFiles = Path(args.libraries_path)/args.organism_name/"index.Libs"/file_pre_tRNA
        bwtCommand = Path(args.bowtie_path)/"bowtie-inspect" if args.bowtie_path else "bowtie-inspect"
        bwtExec = str(bwtCommand) +" -a 20000 -e "+ str(indexFiles)
        bowtie = subprocess.run(str(bwtExec), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
        #READING PRECURSOR miRNA SEQUENCES INFORMATION IN A DICTIONARY (pre_mirDict)
        if bowtie.returncode==0:
            bwtOut2 = bowtie.stdout
            for srow in bwtOut2.split('\n'):
                if srow != "":
                    if '>' in srow:
                        header = srow.replace(">","")
                    else:
                        pretrnaNameSeqDic.update({header:str(srow)})
        else:
            allFiles_inPlace = 0
        # Deal with the alignment of mature tRNA
        #alignmentResult[content[0]].append((content[2], content[1], content[3], content[5]))
        if allFiles_inPlace == 1:
            m_trna.extend(p_trna)
            trfContentDic = {}
            for mtrna_item in m_trna:
                trfContentDic.update({mtrna_item[0]:{'count':mtrna_item[2:]}})
                if "N" not in mtrna_item[0]:
                    trfContentDic[mtrna_item[0]]['uid'] = UID(mtrna_item[0], "tRF") # It is the Unique identifier, this function is present miRgeEssential.py
                else:
                    trfContentDic[mtrna_item[0]]['uid'] = "."
            # Open sam file for the mapped mature tRNA 
            matureMappedtRNA = workDir/"miRge3_tRNA.sam"
            with open(matureMappedtRNA, "r") as minf:
                for mline in minf:
                    mline=mline.strip()
                    #AAAACATCAGATTGTGAGTC    0       trnaMT_HisGTG_MT_+_12138_12206  18      255     20M     *       0       0       AAAACATCAGATTGTGAGTC    IIIIIIIIIIIIIIIIIIII    XA:i:1  MD:Z:17A2       NM:i:1  XM:i:2
                    item = mline.split("\t")
                    startTmp = int(item[3])-1
                    trfContentDic[item[0]][item[2]] = {}
                    trfContentDic[item[0]][item[2]]['start'] = startTmp
                    trfContentDic[item[0]][item[2]]['end'] = startTmp+len(item[0])-1
                    trfContentDic[item[0]][item[2]]['cigar'] = 'undifined'
                    trfContentDic[item[0]][item[2]]['tRFType'] = trfTypes(item[0], item[2], startTmp, trnaStruDic)
            
            primaryMappedtRNA = workDir/"miRge3_pre_tRNA.sam"
            with open(primaryMappedtRNA, "r") as minf:
                for mline in minf:
                    mline=mline.strip()
                    item = mline.split("\t")
                    startTmp = int(item[3])-1
                    trfContentDic[item[0]][item[2]] = {}
                    trfContentDic[item[0]][item[2]]['start'] = startTmp
                    trfContentDic[item[0]][item[2]]['cigar'] = 'undifined'
                    trfContentDic[item[0]][item[2]]['tRFType'] = trfTypes(item[0], item[2], startTmp, trnaStruDic)
                    # Only end coordinate will change
                    cutRemainderSeqLen = re.search('T{3,}$', item[0]).span(0)[0]
                    lenPostTrimming = len(item[0])-cutRemainderSeqLen ## Length after trimming the sequences at end for more than 3 TTT's
                    trfContentDic[item[0]][item[2]]['end'] = startTmp+len(item[0])-1-lenPostTrimming 
            ## CALLING EXTERNAL FUNCTION FROM miRge2 TO OUTPUT THE tRNF RESULT FILES 
            mature_tRNA_Reads_values = list(empty_list[col_vars[1]].values())
            primary_tRNA_Reads_values = list(empty_list[col_vars[2]].values())
            trna_deliverables(args, workDir, pretrnaNameSeqDic, trfContentDic, mature_tRNA_Reads_values, primary_tRNA_Reads_values, trnaAAanticodonDic, base_names, trnaStruDic, duptRNA2UniqueDic, trfMergedList, tRNAtrfDic, trfMergedNameDic)

            #pretrnaNameSeqDic
    summary = pd.DataFrame.from_dict(pre_summary).fillna(0).astype(int)
    summary['Remaining Reads'] = summary['Trimmed Reads (all)'] - (summary[col_tosum].sum(axis=1))
    readDistSample = str(list(pre_summary['Total Input Reads'].keys()))
    readDistGraph = """
        [ { name: 'mature miRNA', data: """+  str(list(pre_summary['Filtered miRNA Reads'].values())) + """},
          { name: 'Hairpin miRNA', data: """+  str(list(pre_summary['Hairpin miRNAs'].values())) + """},  
          { name: 'primary tRNA', data: """+  str(list(pre_summary['primary tRNA Reads'].values())) + """}, 
          { name: 'mature tRNA', data: """+  str(list(pre_summary['mature tRNA Reads'].values())) + """}, 
          { name: 'snoRNA', data: """+  str(list(pre_summary['snoRNA Reads'].values())) + """}, 
          { name: 'rRNA', data: """+  str(list(pre_summary['rRNA Reads'].values())) + """}, 
          { name: 'ncRNA', data: """+  str(list(pre_summary['ncRNA others'].values())) + """}, 
          { name: 'mRNA', data: """+  str(list(pre_summary['mRNA Reads'].values())) + """},  
          { name: 'remaining reads', data: """+  str(summary["Remaining Reads"].tolist()) + """}
        ]
    """
    html_data.readDist(readDistSample, readDistGraph)
    df_rpm = mirRPM_completeSet.reset_index()

    def fmtExprn(name, x, y, val):
        honeyName = ":".join(name.split("-")[1:])
        honeyName = honeyName.split(".")[0]
        honeyName = honeyName.split("/")[0]
        var_item = """
        {
            'hc-a2': '""" + honeyName + """',
            name: '"""+ name + """',
            x: """ + str(x) + """,
            y: """ + str(y) + """,
            value: """+ str(val) + """
        }"""
        return var_item 
    idname =1
    for nme in base_names:
        exprnDivID = "exprnDivID_" + str(idname)
        honey_dataTmp = sorted(df_rpm[['miRNA', nme]].values.tolist(), key=lambda x: x[1], reverse=True)[:40]
        honey_data = sorted(honey_dataTmp)
        tilemapArray = []
        idn = 0
        for ix in range(5):
            for iy in range(8):
                name = honey_data[idn][0]
                val = honey_data[idn][1]
                dataVal = fmtExprn(name, ix, iy, val)
                tilemapArray.append(dataVal)
                idn += 1
        dataSeries = ",".join(tilemapArray)
        html_data.topExprnChart(exprnDivID, nme, dataSeries)
        idname += 1
        #print(dataSeries)
        #print(sorted(honey_data))
        #dict(sorted(mirRPM_completeSet.values.tolist(), key=int, reverse=True)[:40])
    #print(list(summary['Remaining Reads'].values()))
    summary = summary.reindex(columns=colRearrange)
    summary.index.name = "Sample name(s)"
    report = Path(workDir)/"annotation.report.csv"
    report_html = Path(workDir)/"annotation.report.html"
    summary.to_csv(report)
    summary = summary.reset_index(level=['Sample name(s)'])
    summary.index += 1

    data_in_html = summary.to_html(index=False)
    table = '<table style="color="black";font-size:15px; text-align:center; border:0.2px solid black; border-collapse:collapse; table-layout:fixed; height="550"; text-align:center">'
    th = '<th style ="background-color: #3f51b5; color:#ffffff; text-align:center">'
    td ='<td style="vertical-align: middle;background-color: #edf6ff;font-size: 14px;font-family: Arial;font-weight: normal;color: #000000;text-align:center; height="250";border:0; ">'
    data_in_html = data_in_html.replace("<table>", th)
    data_in_html = data_in_html.replace("<th>", th)
    data_in_html = data_in_html.replace("<td>", td)
    with open(report_html,'w') as f:
        f.write(data_in_html)
