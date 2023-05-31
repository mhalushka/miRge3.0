import subprocess
from pathlib import Path
import pandas as pd
import time
import os
import re
import concurrent.futures
from collections import defaultdict
import numpy as np
from mirge.libs.miRgeEssential import UID


def alignPlusParse(bwtExec, iter_number, pdDataFrame, args, workDir):
    """
    ALIGN TO BOWTIE, PARSE SAM FILE AND UPDATE THE DATAFRAME
    """
    #indexNames = ['_mirna_', '_hairpin_', '_mature_trna', '_pre_trna', '_snorna', '_rrna', '_ncrna_others', '_mrna', '_mirna_', '_spike-in']
    colnames = list(pdDataFrame.columns)
    colToAct = 1 + int(iter_number)
    runlogFile = Path(workDir)/"run.log"
    outlog = open(str(runlogFile),"a+")
    if not args.quiet:
        print("Alignment starting for " + colnames[colToAct] + "...")
    outlog.write("Alignment starting for " + colnames[colToAct] + "..." + "\n")
    outlog.close()

    #Generate a file for isomiR output
    exact_alignment = Path(workDir)/"exact_miRNAs.csv"
    isomir_alignment = Path(workDir)/"isomir_alignment.csv"

    bowtie = subprocess.run(str(bwtExec), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
    if args.bam_out:
        if iter_number == 0 or iter_number == 8:
            bwtoutput = Path(workDir)/"miRge3_miRNA.sam"
            bwto = open(bwtoutput, "a+")
        elif iter_number == 1:
            bwtoutput = Path(workDir)/"miRge3_hairpin_miRNA.sam"
            bwto = open(bwtoutput, "a+")
        elif iter_number == 4:
            bwtoutput = Path(workDir)/"miRge3_snorna.sam"
            bwto = open(bwtoutput, "a+")
        elif iter_number == 5:
            bwtoutput = Path(workDir)/"miRge3_rrna.sam"
            bwto = open(bwtoutput, "a+")
        elif iter_number == 6:
            bwtoutput = Path(workDir)/"miRge3_ncrna_others.sam"
            bwto = open(bwtoutput, "a+")
        elif iter_number == 7:
            bwtoutput = Path(workDir)/"miRge3_mrna.sam"
            bwto = open(bwtoutput, "a+")
    if args.tRNA_frag:
        if iter_number == 2:
            bwtoutput = Path(workDir)/"miRge3_tRNA.sam"
            bwto = open(bwtoutput, "a+")
        elif iter_number == 3:
            bwtoutput = Path(workDir)/"miRge3_pre_tRNA.sam"
            bwto = open(bwtoutput, "a+")
    if bowtie.returncode==0:
        bwtOut = bowtie.stdout
        bwtErr = bowtie.stderr
    align_store = defaultdict(list)
    exact_store = defaultdict(list)
    isomir_store = defaultdict(list)
    for srow in bwtOut.split('\n'):
        if not srow.startswith('@'):
            sam_line = srow.split('\t')
            if sam_line != ['']:
                if sam_line[2] != "*":
                    sub_dict = defaultdict(list)
                    if sam_line[0] in align_store:
                        align_store[sam_line[0]][sam_line[2]] = sam_line
                    else:
                        sub_dict[sam_line[2]] = sam_line
                        align_store[sam_line[0]] = sub_dict
                    #Code that needs to be kept in samline iteration
                    elif iter_number != 2 and iter_number != 3 and iter_number != 9: #tRNA and pre_tRNA and spike-ins are ignored
                        if args.bam_out:
                            bwto.write(srow+"\n")
                    elif args.tRNA_frag:
                        if iter_number != 9:
                            bwto.write(srow+"\n")


    #Functionality to handle ambiguous alignments
    mol_count = 0
    multi_match_count = 0
    single_count =0
    name = colnames[colToAct]
    name = name.split(" ")
    name = "_".join(name)
    fname = name + "_secondary_matches.csv"
    full_second_match_file = Path(workDir)/fname
    isomir_alignment_info = Path(workDir)/"isomiR_alignment_info.tsv"
    for keys, values in align_store.items():
        #len values = number of subkeys
        if len(values) != 1:
            multi_match_count +=1
            mol_count += 1
            #Then we have multi-mappers!
            #Append the first alignment to the dataframe for count purposes and mark it as aligned
            top_score = 0
            winner = "null"
            with open(full_second_match_file, "a+") as f:
                for k,v in values.items():
                    #Output results to file
                    v_format = "\t".join(v)
                    f.write(v_format + "\n")
                    #Check if the alignment is the top
                    score = int(v[4])
                    if score > top_score:
                        top_score = score
                        winner = k
                    winning_entry = values[winner]

            #map_name = winning_entry[2]+"*"
            map_name = winning_entry[2]
            pdDataFrame.at[winning_entry[0], colnames[colToAct]] = map_name
            pdDataFrame.at[winning_entry[0], colnames[0]] = 1

        else:
            single_count +=1
            mol_count += 1
            #There will only be a single key and value but without the key (name of matching sequence) this is the only way to access it
            for k,v in values.items():
                pdDataFrame.at[v[0], colnames[colToAct]] = v[2]
                pdDataFrame.at[v[0], colnames[0]] = 1


    if not args.quiet:
        print("Number of molecules with valid alternate alignments for " + colnames[colToAct] + ":" + str(multi_match_count))
    outlog = open(str(runlogFile),"a+")
    outlog.write("Number of molecules with valid alternate alignments for " + colnames[colToAct] + ":" + str(multi_match_count) + "\n")
    outlog.close()
    if multi_match_count !=0 and not args.quiet:
        multi_prop = round(multi_match_count/mol_count * 100, 1)
        print(" This accounts for " + str(multi_prop) + "% of all aligned molecules")
        print(" See " + str(full_second_match_file) + " for alignment details")
        outlog = open(str(runlogFile),"a+")
        outlog.write(" This accounts for " + str(multi_prop) + "% of all aligned molecules" + "\n")
        outlog.write(" See " + str(full_second_match_file) + " for alignment details" + "\n")
        outlog.close()

    if iter_number ==0:
        pd_exact_file  = Path(workDir)/"exact_miRNA_counts.csv"
        pd_exact = pdDataFrame.copy(deep=True)
        pd_slice = pd_exact[pdDataFrame.annotFlag != 0]
        pd_exact_anno = pd_slice.iloc[:,[1]]
        pd_exact_count = pd_slice.iloc[:,11:]
        pd_exact_combine = pd.concat([pd_exact_anno, pd_exact_count], axis=1, join='inner')
        pd_exact_combine.to_csv(pd_exact_file)
    elif iter_number == 8:
        pd_isomir = Path(workDir)/"isomir_miRNA_counts.csv"
        #print(pdDataFrame.iloc[:,[9]])
        pd_iso_change = pdDataFrame.copy(deep=True)
        pd_iso_change['isomiR miRNA'].replace('', np.nan, inplace=True)
        pd_iso_change.dropna(subset=['isomiR miRNA'], inplace=True)
        pd_iso_anno = pd_iso_change.iloc[:,[9]]
        pd_iso_count = pd_iso_change.iloc[:,11:]
        pd_iso_combine = pd.concat([pd_iso_anno, pd_iso_count], axis=1, join='inner')
        pd_iso_combine.to_csv(pd_isomir)


    return pdDataFrame



def bwtAlign(args,pdDataFrame,workDir,ref_db):
    """
    THIS FUNCTION COLLECTS DATAFRAME AND USER ARGUMENTS TO MAP TO VARIOUS DATABASES USING BOWTIE. CALLED FIRST AND ONCE.
    """
    global threads
    threads = args.threads
    begningTime = time.perf_counter()
    bwtCommand = str(Path(args.bowtie_path)/"bowtie ") if args.bowtie_path else "bowtie "
    if args.bowtieVersion == "False": # That is if version is v1.3.0
        bwtCommand += " -x "
    bwtInput = Path(workDir)/"bwtInput.fasta"
    runlogFile = Path(workDir)/"run.log"
    outlog = open(str(runlogFile),"a+")
    if not args.quiet:
        print("Alignment in progress ...")
    outlog.write("Alignment in progress ...\n")
    outlog.close()
    indexNames = ['_mirna_', '_hairpin_', '_mature_trna', '_pre_trna', '_snorna', '_rrna', '_ncrna_others', '_mrna', '_mirna_', '_spike-in']
    parameters = [' -n 0 -f --norc -S --threads ', ' -n 1 -f --norc -S --threads ', ' -v 1 -f -a --best --strata --norc -S --threads ', ' -v 0 -f -a --best --strata --norc -S --threads ', ' -n 1 -f --norc -S --threads ', ' -n 1 -f --norc -S --threads ', ' -n 1 -f --norc -S --threads ', ' -n 0 -f --norc -S --threads ', ' -5 1 -3 2 -v 2 -f --norc --best -S --threads ', ' -n 0 -f --norc -S --threads ']
    if args.spikeIn:
        iterations = 10
    else:
        iterations = 9

    for bwt_iter in range(iterations):
        if bwt_iter == 0:
            with open(bwtInput, 'w') as wseq:
                for sequences in (pdDataFrame.index[pdDataFrame.index.str.len() < 26]):
                    wseq.write(">"+str(sequences)+"\n")
                    wseq.write(str(sequences)+"\n")

            indexName  = str(args.organism_name) + str(indexNames[bwt_iter]) + str(ref_db)
            indexFiles = Path(args.libraries_path)/args.organism_name/"index.Libs"/indexName
            bwtExec = str(bwtCommand) + str(indexFiles) + str(parameters[bwt_iter]) + str(args.threads) + " " + str(bwtInput)
            alignPlusParse(bwtExec, bwt_iter, pdDataFrame, args, workDir)

        elif bwt_iter == 1:
            with open(bwtInput, 'w') as wseq:
                for sequences in (pdDataFrame.index[pdDataFrame.index.str.len() > 25]):
                    wseq.write(">"+str(sequences)+"\n")
                    wseq.write(str(sequences)+"\n")

            indexName  = str(args.organism_name) + str(indexNames[bwt_iter]) + str(ref_db)
            indexFiles = Path(args.libraries_path)/args.organism_name/"index.Libs"/indexName
            bwtExec = str(bwtCommand) + str(indexFiles) + str(parameters[bwt_iter]) + str(args.threads) + " " + str(bwtInput)
            alignPlusParse(bwtExec, bwt_iter, pdDataFrame, args, workDir)

        else:
            if bwt_iter == 8:
                indexName  = str(args.organism_name) + str(indexNames[bwt_iter]) + str(ref_db)
            else:
                indexName  = str(args.organism_name) + str(indexNames[bwt_iter])
            if bwt_iter == 3:
                with open(bwtInput, 'w') as wseq:
                    for sequences in (pdDataFrame.index[pdDataFrame.annotFlag.eq(0)]):
                        try:
                            footer = sequences[:(re.search('T{3,}$', sequences).span(0)[0])]
                            wseq.write(">"+str(sequences)+"\n")
                            wseq.write(str(footer)+"\n")
                        except AttributeError:
                            pass
            else:
                with open(bwtInput, 'w') as wseq:
                    for sequences in (pdDataFrame.index[pdDataFrame.annotFlag.eq(0)]):
                        wseq.write(">"+str(sequences)+"\n")
                        wseq.write(str(sequences)+"\n")

            indexFiles = Path(args.libraries_path)/args.organism_name/"index.Libs"/indexName
            bwtExec = str(bwtCommand) + str(indexFiles) + str(parameters[bwt_iter]) + str(args.threads) + " " + str(bwtInput)
            alignPlusParse(bwtExec, bwt_iter, pdDataFrame, args, workDir)
    finish = time.perf_counter()
    if not args.spikeIn:
        pdDataFrame = pdDataFrame.drop(columns=['spike-in'])

    os.remove(bwtInput)
    pdDataFrame = pdDataFrame.fillna('')
    outlog = open(str(runlogFile),"a+")
    if not args.quiet:
        print(f'Alignment completed in {round(finish-begningTime, 4)} second(s)\n')
    outlog.write(f'Alignment completed in {round(finish-begningTime, 4)} second(s)\n')
    outlog.close()

    return pdDataFrame

