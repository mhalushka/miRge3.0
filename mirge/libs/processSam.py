from Bio import SeqIO

def split_fasta_from_sam(infTmp1, infTmp2, imperfect_FASTA):
    # The arguments of this function are as follows:
    # mapped_mirna_SRR944031_vs_representative_seq_tmp1.sam  mapped_mirna_SRR944031.fa
    filterReadNameList = []
    with open(infTmp1, 'r') as inf1:
        for line in inf1:
            if line[0] != '@':
                temp = line.strip().split('\t')
                if temp[1] in ['0', '16']:
                    readName = line.strip().split('\t')[0]
                    if readName not in filterReadNameList:
                        filterReadNameList.append(readName)
    outf = open(imperfect_FASTA,'w')
    for record in SeqIO.parse(infTmp2, "fasta"):
        if record.id not in filterReadNameList:
            outf.write('>'+record.id+'\n')
            outf.write(str(record.seq)+'\n')
    outf.close()

def combineSam(infTmp1, infTmp2, combined_Sam):
    # The parameters of this function are as follows:
    # mapped_mirna_SRR944031_vs_representative_seq_tmp1.sam  mapped_mirna_SRR944031_vs_representative_seq_tmp2.sam
    outf = open(combined_Sam, 'w')
    with open(infTmp1,'r') as inf1:
        for line in inf1:
            if line[0] == '@':
                outf.write(line)
            else:
                temp = line.strip().split('\t')
                if temp[1] in ['0', '16']:
                    outf.write(line)
                else:
                    pass
    with open(infTmp2, 'r') as inf2:
        for line in inf2:
            if line[0] != '@':
                outf.write(line)
            else:
                pass
    outf.close()

def decorateSam(infTmp1, infTmp2, outfile_modifiedSam, *arg):
    # The parameters of this function are as follows:
    # unmapped_mirna_SRR944031_vs_representative_seq.sam unmapped_mirna_SRR944031.fa mapped_mirna_SRR944031_vs_genome_sorted_clusters_overlap14_trimmed_orig.fa
    seqDic = {}
    for record in SeqIO.parse(infTmp2, "fasta"):
        seqDic.update({record.id:str(record.seq)})
    outf = open(outfile_modifiedSam,"w")
    if len(arg) == 1:
        clusterSeqDic = {}
        for record in SeqIO.parse(arg[0], "fasta"):
            clusterSeqDic.update({record.id:str(record.seq)})
        with open(infTmp1,"r") as inf:
            for line in inf:
                if line[0] == '@':
                    outf.write(line)
                else:
                    miRNAName = line.strip().split('\t')[0]
                    clusterRNAName = line.strip().split('\t')[2]
                    if clusterRNAName in clusterSeqDic.keys():
                        outf.write('\t'.join([miRNAName, miRNAName.split('_')[-1], seqDic[miRNAName], clusterSeqDic[clusterRNAName]])+'\t'+'\t'.join(line.strip().split('\t')[1:])+'\n')
                    else:
                        outf.write('\t'.join([miRNAName, miRNAName.split('_')[-1], seqDic[miRNAName], '*'])+'\t'+'\t'.join(line.strip().split('\t')[1:])+'\n')
    else:
        with open(infTmp1,"r") as inf:
            for line in inf:
                if line[0] == '@':
                    outf.write(line)
                else:
                    miRNAName = line.strip().split('\t')[0]
                    outf.write('\t'.join([miRNAName, miRNAName.split('_')[-1], seqDic[miRNAName]])+'\t'+'\t'.join(line.strip().split('\t')[1:])+'\n') 
    outf.close()

def parse_refine_sam(infTmp, SelectTSV, RevKeptTSV):
    # The argument of this function is as follows:
    # unmapped_mirna_SRR944031_vs_representative_seq_bowtie2_modified.sam
    outf1 = open(SelectTSV,"w")
    outf2 = open(RevKeptTSV, "w")
    with open(infTmp,"r") as inf:
        for line in inf:
            if line[0] != "@":
                content = line.strip().split("\t")
                readName= content[0]
                readCounts = int(content[1])
                readSeq = content[2]
                clusterdSeq = content[3]
                # If the read in unmapped, the clusterdSeq will be "*"
                samFlags = content[4]
                clusterdName = content[5]
                alignedSeq = content[12]
                if samFlags in ["0", "256"]:
                    outf1.write(line)
                    outf2.write(line)
                elif samFlags in ["16", "272"]:
                    outf2.write(line)
                else:
                    pass
    outf1.close()
    outf2.close()
