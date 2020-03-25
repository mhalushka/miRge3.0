import time
import re
import pandas as pd
from pathlib import Path
import subprocess
import os, sys

def sam_header(args):
    if args.organism_name == "human":
        header = "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:248956422\n@SQ\tSN:chr2\tLN:242193529\n@SQ\tSN:chr3\tLN:198295559\n@SQ\tSN:chr4\tLN:190214555\n@SQ\tSN:chr5\tLN:181538259\n@SQ\tSN:chr6\tLN:170805979\n@SQ\tSN:chr7\tLN:159345973\n@SQ\tSN:chr8\tLN:145138636\n@SQ\tSN:chr9\tLN:138394717\n@SQ\tSN:chr10\tLN:133797422\n@SQ\tSN:chr11\tLN:135086622\n@SQ\tSN:chr12\tLN:133275309\n@SQ\tSN:chr13\tLN:114364328\n@SQ\tSN:chr14\tLN:107043718\n@SQ\tSN:chr15\tLN:101991189\n@SQ\tSN:chr16\tLN:90338345\n@SQ\tSN:chr17\tLN:83257441\n@SQ\tSN:chr18\tLN:80373285\n@SQ\tSN:chr19\tLN:58617616\n@SQ\tSN:chr20\tLN:64444167\n@SQ\tSN:chr21\tLN:46709983\n@SQ\tSN:chr22\tLN:50818468\n@SQ\tSN:chrX\tLN:156040895\n@SQ\tSN:chrY\tLN:57227415\n@SQ\tSN:chrM\tLN:16569\n@PG\tID:Bowtie\tVN:1.2.3\tCL:\"bowtie-1.2.3-linux-x86_64/bowtie-align-s --wrapper basic-0 sample_index sample.fastq -S -p 12\"\n"
#        >ENST00000516122.1 ncrna chromosome:GRCh38:12:59450673:59450772:-1 gene:ENSG00000251931.1 gene_biotype:snRNA transcript_biotype:snRNA gene_symbol:RNU6-871P description:RNA, U6 small nuclear 871, pseudogene [Source:HGNC Symbol;Acc:HGNC:47834]
        pass
    elif args.organism_name == "mouse":
        header = "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:195471971\n@SQ\tSN:chr2\tLN:182113224\n@SQ\tSN:chr3\tLN:160039680\n@SQ\tSN:chr4\tLN:156508116\n@SQ\tSN:chr5\tLN:151834684\n@SQ\tSN:chr6\tLN:149736546\n@SQ\tSN:chr7\tLN:145441459\n@SQ\tSN:chr8\tLN:129401213\n@SQ\tSN:chr9\tLN:124595110\n@SQ\tSN:chr10\tLN:130694993\n@SQ\tSN:chr11\tLN:122082543\n@SQ\tSN:chr12\tLN:120129022\n@SQ\tSN:chr13\tLN:120421639\n@SQ\tSN:chr14\tLN:124902244\n@SQ\tSN:chr15\tLN:104043685\n@SQ\tSN:chr16\tLN:98207768\n@SQ\tSN:chr17\tLN:94987271\n@SQ\tSN:chr18\tLN:90702639\n@SQ\tSN:chr19\tLN:61431566\n@SQ\tSN:chrX\tLN:171031299\n@SQ\tSN:chrY\tLN:91744698\n@SQ\tSN:chrM\tLN:16299\n@SQ\tSN:GL456210.1\tLN:169725\n@SQ\tSN:GL456211.1\tLN:241735\n@SQ\tSN:GL456212.1\tLN:153618\n@SQ\tSN:GL456213.1\tLN:39340\n@SQ\tSN:GL456216.1\tLN:66673\n@SQ\tSN:GL456219.1\tLN:175968\n@SQ\tSN:GL456221.1\tLN:206961\n@SQ\tSN:GL456233.1\tLN:336933\n@SQ\tSN:GL456239.1\tLN:40056\n@SQ\tSN:GL456350.1\tLN:227966\n@SQ\tSN:GL456354.1\tLN:195993\n@SQ\tSN:GL456359.1\tLN:22974\n@SQ\tSN:GL456360.1\tLN:31704\n@SQ\tSN:GL456366.1\tLN:47073\n@SQ\tSN:GL456367.1\tLN:42057\n@SQ\tSN:GL456368.1\tLN:20208\n@SQ\tSN:GL456370.1\tLN:26764\n@SQ\tSN:GL456372.1\tLN:28664\n@SQ\tSN:GL456378.1\tLN:31602\n@SQ\tSN:GL456379.1\tLN:72385\n@SQ\tSN:GL456381.1\tLN:25871\n@SQ\tSN:GL456382.1\tLN:23158\n@SQ\tSN:GL456383.1\tLN:38659\n@SQ\tSN:GL456385.1\tLN:35240\n@SQ\tSN:GL456387.1\tLN:24685\n@SQ\tSN:GL456389.1\tLN:28772\n@SQ\tSN:GL456390.1\tLN:24668\n@SQ\tSN:GL456392.1\tLN:23629\n@SQ\tSN:GL456393.1\tLN:55711\n@SQ\tSN:GL456394.1\tLN:24323\n@SQ\tSN:GL456396.1\tLN:21240\n@SQ\tSN:JH584292.1\tLN:14945\n@SQ\tSN:JH584293.1\tLN:207968\n@SQ\tSN:JH584294.1\tLN:191905\n@SQ\tSN:JH584295.1\tLN:1976\n@SQ\tSN:JH584296.1\tLN:199368\n@SQ\tSN:JH584297.1\tLN:205776\n@SQ\tSN:JH584298.1\tLN:184189\n@SQ\tSN:JH584299.1\tLN:953012\n@SQ\tSN:JH584300.1\tLN:182347\n@SQ\tSN:JH584301.1\tLN:259875\n@SQ\tSN:JH584302.1\tLN:155838\n@SQ\tSN:JH584303.1\tLN:158099\n@SQ\tSN:JH584304.1\tLN:114452\n@SQ\tSN:KQ030490.1\tLN:120154\n@SQ\tSN:KB469738.3\tLN:210641\n@SQ\tSN:JH792830.1\tLN:246751\n@SQ\tSN:KK082442.1\tLN:154766\n@SQ\tSN:KQ030495.1\tLN:490000\n@SQ\tSN:KQ030485.1\tLN:185548\n@SQ\tSN:KQ030486.1\tLN:399265\n@SQ\tSN:KQ030487.1\tLN:316842\n@SQ\tSN:KB469739.1\tLN:331111\n@SQ\tSN:KB469741.1\tLN:267241\n@SQ\tSN:JH792829.1\tLN:352455\n@SQ\tSN:JH792826.1\tLN:368286\n@SQ\tSN:JH792828.1\tLN:473738\n@SQ\tSN:KB469740.1\tLN:316140\n@SQ\tSN:JH792832.1\tLN:544189\n@SQ\tSN:JH792833.1\tLN:221588\n@SQ\tSN:JH792827.1\tLN:205713\n@SQ\tSN:JH792834.1\tLN:331480\n@SQ\tSN:JH792831.1\tLN:182256\n@SQ\tSN:KB469742.1\tLN:1059955\n@SQ\tSN:KQ030484.1\tLN:214957\n@SQ\tSN:KQ030494.1\tLN:506812\n@SQ\tSN:KQ030496.1\tLN:260447\n@SQ\tSN:KQ030497.1\tLN:219721\n@SQ\tSN:KQ030492.1\tLN:166012\n@SQ\tSN:KQ030493.1\tLN:269476\n@SQ\tSN:KQ030491.1\tLN:129865\n@SQ\tSN:KQ030488.1\tLN:45901\n@SQ\tSN:KQ030489.1\tLN:188269\n@SQ\tSN:KK082441.1\tLN:456798\n@PG\tID:Bowtie\tVN:1.2.3\tCL:\"bowtie-1.2.3-linux-x86_64/bowtie-align-s --wrapper basic-0 sample_index sample.fastq -S -p 12\"\n"
        pass
    else:
        pass
    return header

def fetchGenCor(args, index_file_name, dict_gen_coordinates):
    bwtCommand = Path(args.bowtie_path)/"bowtie-inspect" if args.bowtie_path else "bowtie-inspect"
    bwtExec = bwtCommand+" -n "+ str(index_file_name)
    bowtie = subprocess.run(str(bwtExec), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
    if bowtie.returncode==0:
        bwtOut = bowtie.stdout
        bwtErr = bowtie.stderr
        for srow in bwtOut.split('\n'):
            if not "PATCH" in srow:
                header_line = srow.split(' ')
                try:
                    dict_gen_coordinates[header_line[0]] = header_line[2]
                except IndexError:
                    pass
    return dict_gen_coordinates
    #ENST00000516122.1 ncrna chromosome:GRCh38:12:59450673:59450772:-1 gene:ENSG00000251931.1 gene_biotype:snRNA transcript_biotype:snRNA gene_symbol:RNU6-871P description:RNA, U6 small nuclear 871, pseudogene [Source:HGNC Symbol;Acc:HGNC:47834]

def bow2bam(args, workDir, ref_db, df_list, base_names, index_file_name, rna_type, samName):
    samtoolsCommandPre = Path(args.samtools_path)/"samtools " if args.samtools_path else "samtools "
    dict_genCors = {}
    dict_genCors = fetchGenCor(args, index_file_name, dict_genCors) # Key: ENST00000516122.1 and value: chromosome:GRCh38:12:59450673:59450772:-1

    bam_expression_dict={}
    for seq in df_list:
        bam_expression_dict[seq[0]] = [int(x) for x in seq[2:]]
    
    req_sam = "miRge3_"+samName+".sam"
    otrna_samFile = Path(workDir)/req_sam
    with open(otrna_samFile) as miSam:
        for mi_idx, mi_sam in enumerate(miSam):
            mi_sam = mi_sam.strip()
            mi_sam_list = mi_sam.split("\t")
            try:
                sam_exprn_list = bam_expression_dict[mi_sam_list[0]] # list of expression values 
                try:
                    chromo = "chr"+ str(dict_genCors[mi_sam_list[2]].split(":")[2])
                    if chromo == "chrMT": 
                        chromo = chromo.replace("chrMT", "chrM")
                    start = int(dict_genCors[mi_sam_list[2]].split(":")[3]) + int(mi_sam_list[3]) - 1
                    strand = dict_genCors[mi_sam_list[2]].split(":")[-1]
                    for ex_idx, exprn in enumerate(sam_exprn_list):
                        if exprn >= 1:
                            file_sam_name = str(base_names[ex_idx]) +".sam"
                            sam_name = Path(workDir)/file_sam_name
                            xbam = open(sam_name, "a+")
                            for numexp in range(exprn):
                                readname = mi_sam_list[0] + "_" + str(numexp)
                                phredQual = "I"*len(mi_sam_list[0])
                                #xbamout = readname+"\t"+mi_sam_list[1]+"\t"+genC+"\t"+str(genS)+"\t"+mi_sam_list[4]+"\t"+mi_sam_list[5]+"\t"+mi_sam_list[6]+"\t"+mi_sam_list[7]+"\t"+mi_sam_list[8]+"\t"+mi_sam_list[9]+"\t"+mi_sam_list[10]+"\n"
                                if strand == "1":
                                    xbamout = readname+"\t"+mi_sam_list[1]+"\t"+chromo+"\t"+str(start)+"\t"+mi_sam_list[4]+"\t"+mi_sam_list[5]+"\t"+mi_sam_list[6]+"\t"+mi_sam_list[7]+"\t"+mi_sam_list[8]+"\t"+mi_sam_list[0]+"\t"+phredQual+"\n"
                                else:
                                    xbamout = readname+"\t"+mi_sam_list[1]+"\t"+chromo+"\t"+str(start)+"\t"+mi_sam_list[4]+"\t"+mi_sam_list[5]+"\t"+mi_sam_list[6]+"\t"+mi_sam_list[7]+"\t"+mi_sam_list[8]+"\t"+mi_sam_list[0][::-1]+"\t"+phredQual+"\n"
                                xbam.write(xbamout)
                            xbam.close()
                except IndexError:                                                                                                                                                                                                    pass
            except KeyError:
                pass


def createBAM(args, workDir, base_names):
    samtoolsCommandPre = Path(args.samtools_path)/"samtools " if args.samtools_path else "samtools "
    # SAM to BAM conversion
    for files_sam in base_names:
        file_sam_name = str(files_sam) +".sam"
        file_bam_name = str(files_sam) +".bam"
        file_Sortbam_name = str(files_sam) +"_sorted.bam"
        file_Sortbam_idx = str(files_sam) +"_sorted.bai"
        sam_name = Path(workDir)/file_sam_name
        bam_name = Path(workDir)/file_bam_name
        bam_sortname = Path(workDir)/file_Sortbam_name
        bam_sortidx = Path(workDir)/file_Sortbam_idx
        samCom2bam = str(samtoolsCommandPre) + "view -bS " + str(sam_name) + " > " + str(bam_name)
        samcreation = subprocess.run(str(samCom2bam), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
        if samcreation.returncode !=0:
            print("Error in creating BAM file!. Ignoring this step\n")
            pass
        else:
            #os.remove(sam_name)
            pass
        bamsort = str(samtoolsCommandPre) + "sort -@ " + str(args.threads) + " " + str(bam_name) + " -o " + str(bam_sortname)
        bamsorting = subprocess.run(str(bamsort), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
        if bamsorting.returncode !=0:
            print("Error in Sorting BAM file!. Ignoring this step\n")
            pass
        else:
            #os.remove(bam_name)
            pass
        bamindex = str(samtoolsCommandPre) + "index " + str(bam_sortname) + " " + str(bam_sortidx)
        bamindexing = subprocess.run(str(bamindex), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
        if bamindexing.returncode !=0:
            print("Error indexing BAM file!. Ignoring this step\n")
            pass
        #mfname = args.organism_name + "_merges_" + ref_db + ".csv"
    #mergeFile = Path(args.libraries_path)/args.organism_name/"annotation.Libs"/mfname
