#Build a class of a read cluster including memeber reads and the clustered seq.
#Define it's attributes and methods

import os
import sys
from pathlib import Path
#import cPickle
import pickle
from Bio.Seq import Seq
from Bio import SeqIO
from mirge.classes.readCluster import ReadCluster


def interset(l1,l2):
    outlist = []
    for item in l1:
        if item in l2:
            outlist.append(item)
    return outlist

def distance(tup1, tup2):
    (s1, e1) = tup1
    (s2, e2) = tup2
    # In default, s1 < s2
    d = s2 - e1-1
    return d

def headDashCount(clusterSeq):
    # clusterSeq is: CTAATACTGCCTGGTAATGATGACGG-
    headDashCount = 0
    for symbol in clusterSeq:
        if symbol == '-':
            headDashCount = headDashCount + 1
        else:
            break
    return headDashCount

def tailDashCount(clusterSeq):
    # clusterSeq is: CTAATACTGCCTGGTAATGATGACGG-
    tailDashCount = 0
    for symbol in clusterSeq[::-1]:
        if symbol == '-':
            tailDashCount = tailDashCount + 1
        else:
            break
    return tailDashCount

def removeDash(seq):
    # clusterSeq is: --CTAATACTGCCTGGTAATGATGACGG-
    outSeq = ''
    for symbol in seq:
        if symbol != '-':
            outSeq = outSeq + symbol
    return outSeq

def generate_featureFiles(outputdir2, infFile, chrSeqDic, chrSeqLenDic, miRNAchrCoordivateDic, exactmiRNASeqDic, readNameIndex=1, readCountIndex=2, readSeqIndex=3, clusterNameIndex=6, clusterSeqIndex=4):
    #def generate_featureFiles(infFile, readNameIndex=1, readCountIndex=2, readSeqIndex=3, clusterNameIndex=6, clusterSeqIndex=4):
    # *.py unmapped_mirna_SRR944034_vs_representative_seq_modified_selected_sorted.tsv hg38.pckl hsa_miRNA.gff3 hsa_mature.fa readNameIndex(1) readCountIndex(2) readSeqIndex(3) clusterNameIndex(6) clusterSeqIndex(4)"
    files_in = str(Path(outputdir2)/(infFile+"_modified_selected_sorted.tsv"))
    readCountLimit = 10
    seqCountLimit = 3
    stableClusterSeqLenLimit = 16-6-3
    closestDistance = 9
    fathestDistance = 44
    # closest distance from miRNA to the terminal of chromosome is set to be 20
    distance2Terminal = 20
    fileName = infFile
    #fileName = os.path.basename(infFile) ### NEEDS CHANGE 
    readNameIndex = readNameIndex-1
    readCountIndex = readCountIndex-1
    readSeqIndex = readSeqIndex-1
    clusterNameIndex = clusterNameIndex-1
    clusterSeqIndex = clusterSeqIndex-1
    chrContentDic = {}
    chrList = []
    chrContentDetailedDic = {}
    with open(files_in, 'r') as inf:
        for line in inf:
            lineContent = line.strip().split('\t')
            readName = lineContent[readNameIndex]
            readCount = int(lineContent[readCountIndex])
            readSeq = lineContent[readSeqIndex]
            clusterName = lineContent[clusterNameIndex]
            clusterSeq = lineContent[clusterSeqIndex]
            chr = clusterName.split(':')[2].strip()
            startPos = int(clusterName.split(':')[-1][:-1].split('_')[0].strip())
            endPos = int(clusterName.split(':')[-1][:-1].split('_')[1].strip())
            if chr not in chrContentDic.keys():
                chrList.append(chr)
                chrContentDic.update({chr:[]})
                chrContentDetailedDic.update({chr:[]})
                chrContentDic[chr].append([startPos, endPos, clusterName, clusterSeq, [readName], [readSeq], [readCount]])
            else:
                if clusterName == chrContentDic[chr][-1][2]:
                    chrContentDic[chr][-1][4].append(readName)
                    chrContentDic[chr][-1][5].append(readSeq)
                    chrContentDic[chr][-1][6].append(readCount)
                else:
                    chrContentDic[chr].append([startPos, endPos, clusterName, clusterSeq, [readName], [readSeq], [readCount]])
    chrList.sort()
    for chr in chrList:
        chrContentDic[chr].sort()

    for chr in chrList:
        for t in range(len(chrContentDic[chr])):
            cluster = chrContentDic[chr][t]
            startPos = cluster[0]
            endPos = cluster[1]
            clusterName = cluster[2]
            clusterSeq = cluster[3]
            readNameList = cluster[4]
            readSeqList = cluster[5]
            readCountList = cluster[6]
            strand = clusterName[-1]
            # Only retain the clusters
            # where the count of reads is no less than readCountLimit,
            # the count of seqs is no less than seqCountLimit
            if sum(readCountList) >= readCountLimit and len(readSeqList) >= seqCountLimit:
                if startPos > distance2Terminal and endPos < chrSeqLenDic[chr]-distance2Terminal:
                    clusterInstance = ReadCluster(chrSeqDic, chr, strand, startPos, endPos, clusterName, clusterSeq, readNameList, readSeqList, readCountList)
                    alignSeqList, exacMatch, clusterSeq, clusterSeqRatioList, majorSeq, majorSeqRationList, headStartPosition, tailStartPosition = clusterInstance.locateStartPosition()
                    if (headStartPosition is None) or (tailStartPosition is None):
                        pass
                    else:
                        chrContentDetailedDic[chr].append((clusterInstance, alignSeqList, exacMatch, clusterSeq, clusterSeqRatioList, majorSeq, majorSeqRationList, headStartPosition, tailStartPosition))

    outf5 = open(str(Path(outputdir2)/(infFile+'_cluster.txt')),'w')
    with open(str(Path(outputdir2)/(infFile+'_features.tsv')),'w') as outf:
        outf.write('realMicRNA\trealMicRNAName\tchr\tstartPos\tendPos\tclusterName\tclusterSeq\tmajoritySeq\tstableClusterSeq\talignedClusterSeq\tadjustedClusterSeq\tclusterSecondSeq\ttemplateSeq\tseqCount\treadCountSum\texactMatchRatio\theadUnstableLength\ttailUnstableLength\t')
        s = 0
        for chr in chrList:
            for t in range(len(chrContentDetailedDic[chr])):
                clusterInstance, alignSeqList, exacMatch, clusterSeq, clusterSeqRatioList, majorSeq, majorSeqRationList, headStartPosition, tailStartPosition = chrContentDetailedDic[chr][t]
                startPos = clusterInstance.startPos
                endPos = clusterInstance.endPos
                clusterName = clusterInstance.clusterName
                readNameList = clusterInstance.readNameList
                readSeqList = clusterInstance.readSeqList
                readCountList = clusterInstance.readCountList
                strand = clusterName[-1]
                alignedClusterSeq = alignSeqList[0]
                
                # To do.... Calculate the neigbors, three categories: good(distance (5+4)-(40+4) and ++ or --), bad(distance <=(8+4) or overlap or +- ), null(distance > (40+4)) 
                if len(chrContentDetailedDic[chr]) >=3:
                    if t == 0:
                        loopLength = distance((startPos+headStartPosition-headDashCount(clusterSeq), endPos+1+tailStartPosition+tailDashCount(clusterSeq)), (chrContentDetailedDic[chr][t+1][0].startPos+chrContentDetailedDic[chr][t+1][7]-headDashCount(chrContentDetailedDic[chr][t+1][3]), chrContentDetailedDic[chr][t+1][0].endPos+1+chrContentDetailedDic[chr][t+1][8]+tailDashCount(chrContentDetailedDic[chr][t+1][3])))
                        if loopLength >=closestDistance and loopLength <= fathestDistance and strand == chrContentDetailedDic[chr][t+1][0].clusterName[-1]:
                            neighborState = 'Good'
                            #neighborState = 'Good:'+str(loopLength)
                        elif loopLength > fathestDistance:
                            neighborState = 'Null'
                            #neighborState = 'Null:'+str(loopLength)
                        else:
                            neighborState = 'Bad'
                            #neighborState = 'Bad:'+str(loopLength)
                        upstreamDistance = None
                        downstreamDistance = loopLength
                    elif t == len(chrContentDetailedDic[chr])-1:
                        loopLength = distance((chrContentDetailedDic[chr][t-1][0].startPos+chrContentDetailedDic[chr][t-1][7]-headDashCount(chrContentDetailedDic[chr][t-1][3]), chrContentDetailedDic[chr][t-1][0].endPos+1+chrContentDetailedDic[chr][t-1][8]+tailDashCount(chrContentDetailedDic[chr][t-1][3])), (startPos+headStartPosition-headDashCount(clusterSeq), endPos+1+tailStartPosition+tailDashCount(clusterSeq)))
                        if loopLength >=closestDistance and loopLength <= fathestDistance and chrContentDetailedDic[chr][t-1][0].clusterName[-1] == strand:
                            neighborState = 'Good'
                            #neighborState = 'Good:'+str(loopLength)
                        elif loopLength > fathestDistance:
                            neighborState = 'Null'
                            #neighborState = 'Null:'+str(loopLength)
                        else:
                            neighborState = 'Bad'
                            #neighborState = 'Bad:'+str(loopLength)
                        upstreamDistance = loopLength
                        downstreamDistance = None
                    else:
                        loopLength1 = distance((chrContentDetailedDic[chr][t-1][0].startPos+chrContentDetailedDic[chr][t-1][7]-headDashCount(chrContentDetailedDic[chr][t-1][3]), chrContentDetailedDic[chr][t-1][0].endPos+1+chrContentDetailedDic[chr][t-1][8]+tailDashCount(chrContentDetailedDic[chr][t-1][3])),(startPos+headStartPosition-headDashCount(clusterSeq), endPos+1+tailStartPosition+tailDashCount(clusterSeq)))
                        loopLength2 = distance((startPos+headStartPosition-headDashCount(clusterSeq), endPos+1+tailStartPosition+tailDashCount(clusterSeq)), (chrContentDetailedDic[chr][t+1][0].startPos+chrContentDetailedDic[chr][t+1][7]-headDashCount(chrContentDetailedDic[chr][t+1][3]), chrContentDetailedDic[chr][t+1][0].endPos+1+chrContentDetailedDic[chr][t+1][8]+tailDashCount(chrContentDetailedDic[chr][t+1][3])))
                        if (loopLength1 >=closestDistance and loopLength1 <= fathestDistance and chrContentDetailedDic[chr][t-1][0].clusterName[-1] == strand) or (loopLength2 >=closestDistance and loopLength2 <= fathestDistance and strand == chrContentDetailedDic[chr][t+1][0].clusterName[-1]):
                            neighborState = 'Good'
                        elif (loopLength1 > fathestDistance) and (loopLength2 > fathestDistance):
                            neighborState = 'Null'
                        else:
                            neighborState = 'Bad'
                        upstreamDistance = loopLength1
                        downstreamDistance = loopLength2
                elif len(chrContentDetailedDic[chr]) ==2:
                    if t == 0:
                        loopLength = distance((startPos+headStartPosition-headDashCount(clusterSeq), endPos+1+tailStartPosition+tailDashCount(clusterSeq)),(chrContentDetailedDic[chr][t+1][0].startPos+chrContentDetailedDic[chr][t+1][7]-headDashCount(chrContentDetailedDic[chr][t+1][3]), chrContentDetailedDic[chr][t+1][0].endPos+1+chrContentDetailedDic[chr][t+1][8]+tailDashCount(chrContentDetailedDic[chr][t+1][3])))
                        if loopLength >=closestDistance and loopLength <= fathestDistance and strand == chrContentDetailedDic[chr][t+1][0].clusterName[-1]:
                            neighborState = 'Good'
                            #neighborState = 'Good:'+str(loopLength)
                        elif loopLength > fathestDistance:
                            neighborState = 'Null'
                            #neighborState = 'Null:'+str(loopLength)
                        else:
                            neighborState = 'Bad'
                            #neighborState = 'Bad:'+str(loopLength)
                        upstreamDistance = None
                        downstreamDistance = loopLength
                    else:
                        loopLength = distance((chrContentDetailedDic[chr][t-1][0].startPos+chrContentDetailedDic[chr][t-1][7]-headDashCount(chrContentDetailedDic[chr][t-1][3]), chrContentDetailedDic[chr][t-1][0].endPos+1+chrContentDetailedDic[chr][t-1][8]+tailDashCount(chrContentDetailedDic[chr][t-1][3])),(startPos+headStartPosition-headDashCount(clusterSeq), endPos+1+tailStartPosition+tailDashCount(clusterSeq)))
                        if loopLength >=closestDistance and loopLength <= fathestDistance and chrContentDetailedDic[chr][t-1][0].clusterName[-1] == strand:
                            neighborState = 'Good'
                            #neighborState = 'Good:'+str(loopLength)
                        elif loopLength > fathestDistance:
                            neighborState = 'Null'
                            #neighborState = 'Null:'+str(loopLength)
                        else:
                            neighborState = 'Bad'
                            #neighborState = 'Bad:'+str(loopLength)
                        upstreamDistance = loopLength
                        downstreamDistance = None
                else:
                    neighborState = 'Null'
                    upstreamDistance = None
                    downstreamDistance = None
                    
                flag = 'Null'
                miRNA = 'Null'
                
                outf5.write("Cluster Name: %s\n"%(clusterName))
                outf5.write("%s (%d - %d%s) cluster %d:\n"%(chr, startPos, endPos, clusterName[-1], t))
                outf5.write("%s: %s\n"%(clusterSeq, '\t'.join(map(lambda x: str(x), clusterSeqRatioList))))
                outf5.write("%s: %s\n"%(majorSeq, '\t'.join(map(lambda x: str(x), majorSeqRationList))))
                for i in range(len(alignSeqList)):
                    if i == 0:
                        outf5.write("%s\n"%(alignSeqList[i]))
                    else:
                        outf5.write('%s\t%d\n'%(alignSeqList[i], readCountList[i-1]))
                features = clusterInstance.calculateFeature()
                keptState = features.pop()
                # Write the header 
                if s == 0:
                    outf.write('\t'.join(features[10])+'\tneighborState\tupstreamDistance\tdownstreamDistance\n')
                else:
                    pass
                s = 1
                # If headUnstableLength or tailUnstableLength is None Or stableClusterSeqLen < stableClusterSeqLenLimit,
                # the cluster sequence will not be written to output file.
                if keptState:
                    stableClusterSeqLen = len(features[2])-features[8]-features[9]
                    if stableClusterSeqLen >= stableClusterSeqLenLimit:
                        features.append(neighborState)
                        features.append(upstreamDistance)
                        features.append(downstreamDistance)
                        #return (self.clusterName, self.clusterSeq, adjustedClusterSeq, templateSeq, seqCount, readCountSum, exactMatchRatio, neclueotidePositionList, neclueotideCountList)
                        # To do...
                        # Output three more columns: stable cluster sequence, most abundant seqence, adjusted stable cluster seuqence.
                        # Test which sequence could reflect the real miRNA.
                        # .....
                        listTmp = []
                        for i in range(len(alignSeqList)-1):
                            listTmp.append([readCountList[i], alignSeqList[i+1]])
                        listTmp.sort(reverse=True)
                        majoritySeq = removeDash(listTmp[0][1])
                        seqTmp = alignSeqList[0]
                        if tailStartPosition == -1:
                            stableClusterSeqTmp = seqTmp[headStartPosition:]
                        else:
                            stableClusterSeqTmp = seqTmp[headStartPosition: tailStartPosition+1]
                        stableClusterSeq = removeDash(stableClusterSeqTmp)
                        
                        majoritySeq, stableClusterSeq
                        outf.write('%s\t%s\t'%(flag, miRNA))
                        outf.write('\t'.join([chr, str(startPos), str(endPos), features[0], features[1], majoritySeq, stableClusterSeq, alignedClusterSeq, str(features[2]), str(features[3]), str(features[4]), str(features[5]), str(features[6]), str(features[7]), str(features[8]), str(features[9])]))
                        outf.write('\t')
                        for i in range(len(features[11])):
                            outf.write(str(features[11][i])+'\t')
                        outf.write(features[12]+'\t')
                        if upstreamDistance is None:
                            outf.write('None'+'\t')
                        else:
                            outf.write(str(features[13])+'\t')
                        if downstreamDistance is None:
                            outf.write('None'+'\n')
                        else:
                            outf.write(str(features[14])+'\n')

