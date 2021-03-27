#!/usr/bin/env python
import time
import math
import random
import re
import pandas as pd
import subprocess
from pathlib import Path
import numpy as np
import os, sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
# from Bio.Alphabet import IUPAC, Gapped
from scipy import stats

def addDashNew(seq, totalLength, start, end):
    newSeq = '-'*(start-1)+seq+'-'*(totalLength-end)
    return newSeq


def assign_cluster(dashedSeq, tRNAName, tRNAtrfDic):
    disCutoff = 8.0
    disList = []
    if tRNAName in tRNAtrfDic.keys():
        for seq2Search in tRNAtrfDic[tRNAName].keys():
            disList.append((getDistance2(dashedSeq, seq2Search), tRNAtrfDic[tRNAName][seq2Search]))
        disList.sort()
        if disList[0][0] <= disCutoff:
            assignedName = disList[0][1]
            nearestCluster = disList[0][1]
        else:
            assignedName = 'Undef'
            nearestCluster = disList[0][1]
        dist = disList[0][0]
    else:
        assignedName = 'Dele'
        dist = 100
        nearestCluster = 'Null'
    return (assignedName, dist, nearestCluster)


def getDistance2(seq1, seq2):
    '''
    Calculate the editing distance between two seqs.
    seq1: GCATTGTTGGTTCAGTGGTAGAATTCTCGCCTGC----------------------------------------
    seq2: ---------------TGGTAGAATGCTCGCCTGCCACGCGGGA-------------------------------
    the distance between seq1 and seq2 is calculated by:
    weight1*(number of front insertion)+weight2*(number of rear insertion)+weight3*(number of substitution)
    '''
    weight1, weight2, weight3 = 1, 1, 1
    coordinate1 = coordinate(seq1)
    coordinate2 = coordinate(seq2)
    substiCount = 0
    for k in range(len(seq1)):
        try:
            if seq1[k] != '-' and seq2[k] != '-' and seq1[k] != seq2[k]:
                substiCount = substiCount + 1
        except IndexError:
            substiCount = substiCount + 1
    distance = weight1*abs(coordinate1[0]-coordinate2[0])+weight2*abs(coordinate1[1]-coordinate2[1])+weight3*substiCount
    return distance


def coordinate(dashedSeq):
    startPos = 1
    for i in dashedSeq:
        if i != '-':
            break
        else:
            startPos = startPos + 1
    endPosTmp = 0
    for i in dashedSeq[::-1]:
        if i != '-':
            break
        else:
            endPosTmp = endPosTmp + 1
    endPos = len(dashedSeq) - endPosTmp
    return (startPos, endPos)


def locationOfStartEnd(dashedSeq):
    leftCount = 0
    for i in dashedSeq:
        if i == '-':
            leftCount = leftCount + 1
        else:
            break
    rightCount = 0
    for i in dashedSeq[::-1]:
        if i == '-':
            rightCount = rightCount + 1
        else:
            break
    return (leftCount, rightCount)


def load_data_new(in_file):
    readInforList = []
    with open(in_file, 'r') as inf:
        content = inf.readlines()
    indexLabel = []
    indexLabel.append(-1)
    for index, item in enumerate(content):
        if "mature tRNA" in item or "primary tRNA trailer" in item:
            indexLabel.append(index)
    if len(indexLabel) > 0:
        for i in range(len(indexLabel)-1):
            readInforDic = {}
            j = 1
            for line in content[indexLabel[i]+1+1:indexLabel[i+1]+1-1]:
                data = line.strip().split('\t')
                readInforDic[j] = {'allignedSeq':data[0], 'type':data[1], 'count':int(data[2]), 'RPM':float(data[3])}
                j = j + 1
            if len(readInforDic) != 0:
                readInforList.append(readInforDic)
    else:
        pass
    return readInforList


def local_density(distanceDic, readInforDic, max_id, dc, guass=True, cutoff=False):
    '''
    Compute all points' local density
    Args:
        distanceDic : distance dict
        readInforDic: readInfor dict
        max_id      : max continues id
        dc          : threshold of distance
        gauss       : use guass func or not(can't use together with cutoff)
        cutoff      : use cutoff func or not(can't use together with guass)
    
    Returns:
        local density vector that index is the point index that start from 1
    '''
    assert guass ^ cutoff
    guass_func = lambda dij, dc: math.exp(- (dij / dc) ** 2)
    cutoff_func = lambda dij, dc: 1 if dij < dc else 0
    func = guass and guass_func or cutoff_func
    rho = [-1] + [0] * max_id
    for i in range(1, max_id):
        for j in range(i + 1, max_id + 1):
            rho[i] += func(distanceDic[(i, j)], dc)*readInforDic[j]['RPM']
            rho[j] += func(distanceDic[(i, j)], dc)*readInforDic[i]['RPM']
    # For each point, the number of the read count(RPM) must to added to the rho value
    for i in range(1, max_id+1):
        rho[i] = rho[i]+readInforDic[i]['RPM']
    return np.array(rho, np.float32)


def getDistance(readInforDic):
    '''
    Recently, we don't allow deletion in Bowtie.
    seq1: GCATTGTTGGTTCAGTGGTAGAATTCTCGCCTGC----------------------------------------
    seq2: ---------------TGGTAGAATGCTCGCCTGCCACGCGGGA-------------------------------
    the distance between seq1 and seq2 is calculated by:
    weight1*(number of front insertion)+weight2*(number of rear insertion)+weight3*(number of substitution)
    '''
    weight1, weight2, weight3 = 1.0, 1.0, 1.0
    distanceDic = {}
    min_dis, max_dis, max_id = sys.float_info.max, 0.0, 0
    idList = list(readInforDic.keys())
    for i in range(0, len(idList)):
        for j in range(i+1, len(idList)):
            x1 = idList[i]
            x2 = idList[j]
            seq1 = readInforDic[x1]['allignedSeq']
            seq2 = readInforDic[x2]['allignedSeq']
            coordinate1 = coordinate(seq1)
            coordinate2 = coordinate(seq2)
            substiCount = 0
            for k in range(len(seq1)):
                if seq1[k] != '-' and seq2[k] != '-' and seq1[k] != seq2[k]:
                    substiCount = substiCount + 1
            distance = weight1*abs(coordinate1[0]-coordinate2[0])+weight2*abs(coordinate1[1]-coordinate2[1])+weight3*substiCount
            min_dis, max_dis = min(min_dis, distance), max(max_dis, distance)
            distanceDic[(x1, x2)] = distance
            distanceDic[(x2, x1)] = distance
    max_id = max(idList)
    for i in range(1, max_id+1):
        distanceDic[(i, i)] = 0.0
    return (distanceDic, max_dis, min_dis, max_id)


def min_distance(distanceDic, max_dis, max_id, rho):
    '''
    Compute all points' min distance to the higher local density point(which is the nearest neighbor)
    Args:
        max_id    : max continues id
        max_dis   : max distance for all points
        distanceDic : distance dict
        rho       : local density vector that index is the point index that start from 1
    Returns:
        min_distance vector, nearest neighbor vector
    '''
    sort_rho_idx = np.argsort(-rho)
    delta, nneigh, nneighDistance = [0.0] + [float(max_dis)]*(len(rho)-1), [0]*len(rho), [0.0]*len(rho)
    delta[sort_rho_idx[0]] = -1.0
    for i in range(1, max_id):
        for j in range(0, i):
            old_i, old_j = sort_rho_idx[i], sort_rho_idx[j]
            if distanceDic[(old_i, old_j)] <= delta[old_i]:
                delta[old_i] = distanceDic[(old_i, old_j)]
                nneigh[old_i] = old_j
                nneighDistance[old_i] = distanceDic[(old_i, old_j)]
    delta[sort_rho_idx[0]] = max(delta)
    return np.array(delta, np.float32), np.array(nneigh, np.int32), np.array(nneighDistance, np.float32), sort_rho_idx.astype(np.int32)


def removeDash(seq):
    headDashCount, tailDashCount = DashCount(seq)
    return seq[headDashCount:len(seq)-tailDashCount]


def detectMismach(target_seq, template_seq, position_tmp):
    mimatchState = 'N'
    mismatchPos = []
    startPos, endPos = int(position_tmp.split(':')[0])-1, int(position_tmp.split(':')[1])-1
    seqTmp = template_seq[startPos:endPos+1]
    for i in range(len(target_seq)):
        try:
            if target_seq[i] != seqTmp[i]:
                mimatchState = 'Y'
                mismatchPos.append(str(startPos+1+i))
        except IndexError:
            pass
    return (mimatchState, ','.join(mismatchPos))


def DashCount(seq):
    headDashCount = 0
    tailDashCount = 0
    for i in seq:
        if i == '-':
            headDashCount = headDashCount + 1
        else:
            break
    for i in range(len(seq)-1,-1,-1):
        if seq[i] == '-':
            tailDashCount = tailDashCount + 1
        else:
            break
    return (headDashCount, tailDashCount)


def align2TargetSeq(targetSeq, seqList):
    #align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
    alignSeqList = []
    stateList = []
    for i in range(len(seqList)):
        # Perform pairwise local alignment. And increase the penalty of gap opening and extending into -10 and -10,
        # so that there will be no gap in the alignned result.
        seq = seqList[i]
        alignTmp = pairwise2.align.localms(targetSeq, seq, 2, -1, -20, -20)
        targetSeqTmpNew = alignTmp[0][0]
        seqTmp = alignTmp[0][1]
        
        state = judgeAllign(targetSeqTmpNew, seqTmp)
        stateList.append(state)

        if len(alignSeqList) == 0:
            alignSeqList.append(targetSeqTmpNew)
            alignSeqList.append(seqTmp)
        else:
            if targetSeqTmpNew == alignSeqList[0]:
                alignSeqList.append(seqTmp)
            else:
                headAdd1 = targetSeqTmpNew.index(targetSeq)
                tailAdd1 = len(targetSeqTmpNew)-headAdd1-len(targetSeq)
                headAdd2 = alignSeqList[0].index(targetSeq)
                tailAdd2 = len(alignSeqList[0])-headAdd2-len(targetSeq)
                if headAdd1 >= headAdd2 and tailAdd1 >= tailAdd2:
                    for j in range(len(alignSeqList)):
                        alignSeqList[j] = (headAdd1 - headAdd2)*'-'+alignSeqList[j]+(tailAdd1 - tailAdd2)*'-'
                    alignSeqList.append(seqTmp)
                elif headAdd1 >= headAdd2 and tailAdd1 < tailAdd2:
                    for j in range(len(alignSeqList)):
                        alignSeqList[j] = (headAdd1 - headAdd2)*'-'+alignSeqList[j]
                    alignSeqList.append(seqTmp+(tailAdd2 - tailAdd1)*'-')
                elif headAdd1 < headAdd2 and tailAdd1 >= tailAdd2:
                    for j in range(len(alignSeqList)):
                        alignSeqList[j] = alignSeqList[j]+(tailAdd1 - tailAdd2)*'-'
                    alignSeqList.append((headAdd2 - headAdd1)*'-'+seqTmp)
                elif headAdd1 < headAdd2 and tailAdd1 < tailAdd2:
                    alignSeqList.append((headAdd2 - headAdd1)*'-'+seqTmp+(tailAdd2 - tailAdd1)*'-')
                else:
                    pass
    return (alignSeqList, stateList)


def startAndEnd(seq1, seq2):
    # seq1 is the substring of seq2
    start = seq2.index(seq1)
    end = start+len(seq1)-1
    return (start, end)
    

def judgeAllign(targetSeqTmpNew, seqTmp):
    mismatchCount = 0
    matchCount = 0
    mismatchLimit = 1
    headShift = 1
    tailShift = 3
    state = True
    headDash_target, tailDash_target = DashCount(targetSeqTmpNew)
    headDash_seqTmp, tailDash_seqTmp = DashCount(seqTmp)
    len1 = len(targetSeqTmpNew) - headDash_target - tailDash_target
    matchLimit = len1 - tailShift - mismatchLimit
    #len2 = len(seqTmp) - headDash_seqTmp - tailDash_seqTmp
    start_pos1 = headDash_target
    end_pos1 = len(targetSeqTmpNew)-headDash_target-1-tailShift
    start_pos2 = headDash_seqTmp
    end_pos2 = len(seqTmp)-tailDash_seqTmp-1
    if start_pos2-start_pos1 > headShift:
    	state = False
    else:
    	for pos in range(headDash_target, min([end_pos1, end_pos2])+1):
    		if seqTmp[pos] == '-':
    			pass
    		elif targetSeqTmpNew[pos] != seqTmp[pos]:
    			mismatchCount = mismatchCount + 1
    		else:
    			matchCount = matchCount + 1
    	if start_pos2-start_pos1 == headShift:
    		matchLimitNew = matchLimit - headShift
    	else:
    		matchLimitNew = matchLimit
    	if mismatchCount > mismatchLimit:
    		state = False
    	if matchCount < matchLimitNew:
    		state = False
    return state


def A2IEditing(targetSeq, seqList, countList, mirName, outfTmp, retainedSeqDic, startBase, endBase):
    p_mismatch = 0.001
    tailShift = 5
    alignSeqListTmp, stateListTmp = align2TargetSeq(targetSeq, seqList)
    startLabel, endLabel = startAndEnd(targetSeq, alignSeqListTmp[0])
    targetSeqAligned = alignSeqListTmp[0]
    targetSeqAlignedHeadCount = DashCount(targetSeqAligned)[0]
    #count = 0
    canonicalSeqCount = 0
    countSumTrue = 0
    seqCountTrue = 0
    a2IPositionCountDic = {}
    positionList = []
    alignSeqListKeptNew = []

    for j, seq in enumerate(alignSeqListTmp[1:]):
    	if stateListTmp[j] and removeDash(seq) in retainedSeqDic:
    	#if stateListTmp[j]:
    		if seqList[j] in targetSeq:
    			canonicalSeqCount = canonicalSeqCount + countList[j]
    		alignSeqListKeptNew.append(seq)
    		seqCountTrue = seqCountTrue + 1
    		countSumTrue = countSumTrue + countList[j]
    		for i in range(startLabel, endLabel+1-tailShift):
    			try:
    				if targetSeqAligned[i] == startBase and alignSeqListTmp[j+1][i] == endBase:
    					positionNew = i+1-targetSeqAlignedHeadCount
    					if positionNew not in positionList:
    						positionList.append(positionNew)
    						a2IPositionCountDic.update({positionNew:countList[j]})
    					else:
    						a2IPositionCountDic[positionNew] = a2IPositionCountDic[positionNew] + countList[j]
    					#count = count + countList[j]
    					#break
    			except IndexError:
    				pass
    #print 'Canonical_Seq of %s: %s'%(mirName, targetSeq)
    #print 'seqList size is: %d, %d'%(len(seqList), len(countList))
    outfTmp.write('Canonical_Seq of %s: %s\n'%(mirName, targetSeq))
    outfTmp.write('seqList size is: %d, %d\n'%(len(seqList), len(countList)))
    for index, itemTmp in enumerate(alignSeqListTmp):
    	if index != 0:
    		#print '\t'.join([itemTmp, str(countList[index-1]), str(stateListTmp[index-1])])
    		outfTmp.write('\t'.join([itemTmp, str(countList[index-1]), str(stateListTmp[index-1])]))
    		outfTmp.write('\n')
    #print '****************'
    #print 'retained seqList size is: %d'%(seqCountTrue)
    outfTmp.write('****************\n')
    outfTmp.write('retained seqList size is: %d\n'%(seqCountTrue))
    alignSeqListKept = []
    for index, itemTmp in enumerate(alignSeqListTmp):
    	if index != 0:
    		if stateListTmp[index-1]:
    			alignSeqListKept.append(itemTmp)
    			#print '\t'.join([itemTmp, str(countList[index-1]), str(stateListTmp[index-1])])
    			outfTmp.write('\t'.join([itemTmp, str(countList[index-1]), str(stateListTmp[index-1])]))
    			outfTmp.write('\n')
    outfTmp.write('****************\n')
    outfTmp.write('retained sequences after filering are:\n')
    alignSeqListKept2 = []
    for index, itemTmp in enumerate(alignSeqListTmp):
    	if index != 0:
    		if stateListTmp[index-1] and removeDash(itemTmp) in retainedSeqDic:
    			alignSeqListKept2.append(itemTmp)
    			#print '\t'.join([itemTmp, str(countList[index-1]), str(stateListTmp[index-1])])
    			outfTmp.write('\t'.join([itemTmp, str(countList[index-1]), str(stateListTmp[index-1])]))
    			outfTmp.write('\n')

    a2IPositionRatioDic = {}
    a2IPositionPvalueDic = {}
    for position in positionList:
    	try:
    		a2IPositionRatioDic.update({position:float(a2IPositionCountDic[position])/countSumTrue})
    	except ZeroDivisionError:
    		a2IPositionRatioDic.update({position:0})
    	if countSumTrue - a2IPositionCountDic[position] >= 0:
    		p_value = stats.binom.cdf(countSumTrue - a2IPositionCountDic[position], countSumTrue, 1-p_mismatch)
    	else:
    		p_value = 1.0
    	a2IPositionPvalueDic.update({position:p_value})
    #try:
    #	aIValue = float(count)/countSumTrue
    #except ZeroDivisionError:
    #	aIValue = 0
    return (alignSeqListKeptNew, positionList, a2IPositionCountDic, a2IPositionRatioDic, a2IPositionPvalueDic, countSumTrue, seqCountTrue, canonicalSeqCount)
    #return (aIValue, countSumTrue, seqCountTrue)
    #return (float(count)/sum(countList), count, seqCountTrue)


def mismatchCountAnalysis(targetSeq, seqList, countList, retainedSeqDic):
    tailShift = 5
    alignSeqListTmp, stateListTmp = align2TargetSeq(targetSeq, seqList)
    #print 'Canonical_Seq of %s: %s'%(mirName, targetSeq)
    #print 'seqList size is: %d, %d'%(len(seqList), len(countList))
    #outfTmp.write('Canonical_Seq of %s: %s\n'%(mirName, targetSeq))
    #outfTmp.write('seqList size is: %d, %d\n'%(len(seqList), len(countList)))
    #for index, itemTmp in enumerate(alignSeqListTmp):
    #	if index != 0:
    #		#print '\t'.join([itemTmp, str(countList[index-1]), str(stateListTmp[index-1])])
    #		outfTmp.write('\t'.join([itemTmp, str(countList[index-1]), str(stateListTmp[index-1])]))
    #		outfTmp.write('\n')
    startLabel, endLabel = startAndEnd(targetSeq, alignSeqListTmp[0])
    targetSeqAligned = alignSeqListTmp[0]
    targetSeqAlignedHeadCount = DashCount(targetSeqAligned)[0]
    basePairMismachCountList = []
    for basePair in [('A', 'G'),('A', 'C'),('A', 'T'),('T', 'G'),('T', 'A'),('T', 'C'),('C', 'G'),('C', 'A'),('C', 'T'),('G', 'A'),('G', 'C'),('G', 'T')]:
    	startBaseTmp = basePair[0]
    	endBaseTmp = basePair[1]
    	#startBaseTmp = 'A'
    	#endBaseTmp = 'G'
    	#count = 0
    	countSum = 0
    	seqCount = 0
    	countSumTrue = 0
    	seqCountTrue = 0
    	countSumTrue2 = 0
    	seqCountTrue2 = 0
    	
    	mismatchCountDic_raw = {}
    	positionList_raw = [] 
    	mismatchCountDic = {}
    	positionList = []
    	mismatchCountDic_filter = {}
    	positionList_filter = []

    	for j, seq in enumerate(alignSeqListTmp[1:]):
    		seqCount = seqCount + 1
    		countSum = countSum + countList[j]
    		for i in range(startLabel, endLabel+1-tailShift):
    			try:
    				if targetSeqAligned[i] == startBaseTmp and alignSeqListTmp[j+1][i] == endBaseTmp:
    					positionNew = i+1-targetSeqAlignedHeadCount
    					if positionNew not in positionList_raw:
    						positionList_raw.append(positionNew)
    						mismatchCountDic_raw.update({positionNew:countList[j]})
    					else:
    						mismatchCountDic_raw[positionNew] = mismatchCountDic_raw[positionNew] + countList[j]
    					#count = count + countList[j]
    					#break
    			except IndexError:
    				pass
    		if stateListTmp[j]:
    			seqCountTrue = seqCountTrue + 1
    			countSumTrue = countSumTrue + countList[j]
    			for i in range(startLabel, endLabel+1-tailShift):
    				try:
    					if targetSeqAligned[i] == startBaseTmp and alignSeqListTmp[j+1][i] == endBaseTmp:
    						positionNew = i+1-targetSeqAlignedHeadCount
    						if positionNew not in positionList:
    							positionList.append(positionNew)
    							mismatchCountDic.update({positionNew:countList[j]})
    						else:
    							mismatchCountDic[positionNew] = mismatchCountDic[positionNew] + countList[j]
    						#count = count + countList[j]
    						#break
    				except IndexError:
    					pass
    		if stateListTmp[j] and removeDash(seq) in retainedSeqDic:
    			seqCountTrue2 = seqCountTrue2 + 1
    			countSumTrue2 = countSumTrue2 + countList[j]
    			for i in range(startLabel, endLabel+1-tailShift):
    				try:
    					if targetSeqAligned[i] == startBaseTmp and alignSeqListTmp[j+1][i] == endBaseTmp:
    						positionNew = i+1-targetSeqAlignedHeadCount
    						if positionNew not in positionList_filter:
    							positionList_filter.append(positionNew)
    							mismatchCountDic_filter.update({positionNew:countList[j]})
    						else:
    							mismatchCountDic_filter[positionNew] = mismatchCountDic_filter[positionNew] + countList[j]
    						#count = count + countList[j]
    						#break
    				except IndexError:
    					pass
    	basePairMismachCount_raw = 0
    	basePairMismachCount = 0
    	basePairMismachCount_filter = 0
    	for key in mismatchCountDic_raw.keys():
    		basePairMismachCount_raw = basePairMismachCount_raw + mismatchCountDic_raw[key]
    	for key in mismatchCountDic.keys():
    		basePairMismachCount = basePairMismachCount + mismatchCountDic[key]
    	for key in mismatchCountDic_filter.keys():
    		basePairMismachCount_filter = basePairMismachCount_filter + mismatchCountDic_filter[key]
    	basePairMismachCountList.append((basePairMismachCount_raw, basePairMismachCount, basePairMismachCount_filter))
    return basePairMismachCountList



def trna_deliverables(args, workDir, pretrnaNameSeqDic, trfContentDic, mature_tRNA_Reads_values, primary_tRNA_Reads_values, trnaAAanticodonDic, sampleList, trnaStruDic, duptRNA2UniqueDic, trfMergedList, tRNAtrfDic, trfMergedNameDic):
    # update RPM value for tRFs. Here since tRFs account for around 1.6~2.7% of total reads. 
    # We choose Reads Per 100K mapped reads to normalize the read count.
    for trf in trfContentDic.keys():
        RPMList = []
        readCountListTmp = trfContentDic[trf]['count']
        for i in range(len(sampleList)):
            try:
                RPMList.append((100000.0*readCountListTmp[i])/(mature_tRNA_Reads_values[i]+primary_tRNA_Reads_values[i]))
            except ZeroDivisionError:
                RPMList.append(0.0)
        trfContentDic[trf]['RPM'] = RPMList
    
    trfFile = Path(workDir)/"tRFs.aligned.report.tsv"
    with open(trfFile, 'w') as outf:
        outf.write('read sequence\tuid\tread count(%s)\tRP100K (%s)\tamino acid all hits\tamino acid-anticodon all hits\ttRF information all hits\tamino acid all deduplicated hits\tamino acid-anticodon all deduplicated hits\ttRF information all deduplicated hits\tamino acid one hit\tamino acid-anticodon one hit\ttRF information one hit\n'%(';'.join(sampleList), ';'.join(sampleList)))
        for seqKey in trfContentDic.keys():
            tRFInforList = []
            aaAnticodonList = []
            aaTypeList = []
            candidatetRNAlist = []
            for subkey in trfContentDic[seqKey].keys():
                if subkey not in ['uid', 'RPM', 'count']:
                    tRFInfor = ':'.join([subkey, trfContentDic[seqKey][subkey]['tRFType'], str(trfContentDic[seqKey][subkey]['start']+1), str(trfContentDic[seqKey][subkey]['end']+1)])
                    tRFInforList.append(tRFInfor)
                    candidatetRNAlist.append(subkey)
                    anticodon = trnaAAanticodonDic[subkey]['anticodon']
                    aaTypes = trnaAAanticodonDic[subkey]['aaType']
                    if 'pre_' in subkey:
                        aaTypes = 'pre:'+aaTypes
                    if aaTypes not in ('Und', 'pre:Und'):
                        if aaTypes+'-'+anticodon not in aaAnticodonList:
                            aaAnticodonList.append(aaTypes+'-'+anticodon)
                        if aaTypes not in aaTypeList:
                            aaTypeList.append(aaTypes)
            candidatetRNAUniquelist = []
            for candidatetRNA in candidatetRNAlist:
                try:
                    candidatetRNAUniquelist.append(duptRNA2UniqueDic[candidatetRNA])
                except KeyError:
                    candidatetRNAUniquelist.append(candidatetRNA)
            candidatetRNAUniquelist = list(set(candidatetRNAUniquelist))
            # select a tRNA type from candidatetRNAUniquelist randomly
            # remove the redundant tRNA names and keep the selected one.
            if len(candidatetRNAUniquelist) > 0:
                slectedtRNA = random.choice(candidatetRNAUniquelist)
                outf.write('\t'.join([seqKey, trfContentDic[seqKey]['uid'], ';'.join([str(s) for s in trfContentDic[seqKey]['count']]), ';'.join(['%.3f'%(round(s, 3)) for s in trfContentDic[seqKey]['RPM']])]))
                outf.write('\t')
                outf.write('\t'.join([','.join(aaTypeList), ','.join(aaAnticodonList), ','.join(tRFInforList)]))
                outf.write('\t')
                aaTypeDeduplicatedList = []
                aaAnticodonDeduplicatedList = []
                tRFInforDeduplicatedList = []
                for tRNAtmp in candidatetRNAUniquelist:
                    anticodonTmp = trnaAAanticodonDic[tRNAtmp]['anticodon']
                    aaTypeTmp = trnaAAanticodonDic[slectedtRNA]['aaType']
                    if 'pre_' in tRNAtmp:
                        aaTypeTmp = 'pre:'+aaTypeTmp
                    anticodonTmp = aaTypeTmp + '-' + anticodonTmp
                    aaTypeDeduplicatedList.append(aaTypeTmp)
                    aaAnticodonDeduplicatedList.append(anticodonTmp)
                    tRFInforTmp = ':'.join([tRNAtmp, trfContentDic[seqKey][tRNAtmp]['tRFType'], str(trfContentDic[seqKey][tRNAtmp]['start']+1), str(trfContentDic[seqKey][tRNAtmp]['end']+1)])
                    tRFInforDeduplicatedList.append(tRFInforTmp)
                outf.write('\t'.join([','.join(aaTypeDeduplicatedList), ','.join(aaAnticodonDeduplicatedList), ','.join(tRFInforDeduplicatedList)]))
                outf.write('\t')
                anticodonOneHit = trnaAAanticodonDic[slectedtRNA]['anticodon']
                aaTypeOneHit = trnaAAanticodonDic[slectedtRNA]['aaType']
                if 'pre_' in slectedtRNA:
                    aaTypeOneHit = 'pre:'+aaTypeOneHit
                anticodonOneHit = aaTypeOneHit + '-' + anticodonOneHit
                outf.write('\t'.join([aaTypeOneHit, anticodonOneHit]))
                outf.write('\t')
                outf.write(':'.join([slectedtRNA, trfContentDic[seqKey][slectedtRNA]['tRFType'], str(trfContentDic[seqKey][slectedtRNA]['start']+1), str(trfContentDic[seqKey][slectedtRNA]['end']+1)]))
                outf.write('\n')
                for subitem2 in list(trfContentDic[seqKey].keys()):
                    if subitem2 not in ('uid', 'count', 'RPM', slectedtRNA):
                        del trfContentDic[seqKey][subitem2]
            else:
                del trfContentDic[seqKey]
    #outf3.close()
    #outf4.close()
    # Generate the report table of read count and RP100K values for tRF entities of all the input samples.
    sampletRFEntityDic = {}
    sampletRFEntitySummaryDic = {}
    # Initialize sampletRFEntityDic and sampletRFEntitySummaryDic
    for sample in sampleList:
        sampletRFEntityDic.update({sample:{}})
        for trfMerged in trfMergedList:
            sampletRFEntityDic[sample].update({trfMerged:[0, 0.0]})
            # in [0, 0.0], the first one stores read count and the second on stores RP100K
        sampletRFEntitySummaryDic.update({sample:[0, 0]})
        # in [0, 0], the first one stores the sum of discarded reads count and the second one stores the sum of the total reads (which aligned to tRNA libraries)
    with open(trfFile, 'r') as inf:
        line = inf.readline()
        line = inf.readline()
        while line != '':
            content = line.strip().split('\t')
            seq = content[0]
            tRFinforTmp = content[-1].split(':')
            tRNAName = tRFinforTmp[0]
            start = int(tRFinforTmp[2])
            end = int(tRFinforTmp[3])
            if 'pre' not in tRNAName:
                tRNALength = len(trnaStruDic[tRNAName]['seq'])
            else:
                tRNALength = len(pretrnaNameSeqDic[tRNAName])
            dashedSeq = addDashNew(seq, tRNALength, start, end)
            readcountListTmp = [int(tmp) for tmp in content[2].split(';')]
            RP100KListTmp = [float(tmp) for tmp in content[3].split(';')]
            assignedName, nearestDist, nearestCluster = assign_cluster(dashedSeq, tRNAName, tRNAtrfDic)
            for i, sample in enumerate(sampleList):
                sampletRFEntitySummaryDic[sample][1] += readcountListTmp[i]
                if assignedName not in ['Undef', 'Dele']:
                    mergedName = trfMergedNameDic[assignedName]
                    sampletRFEntityDic[sample][mergedName][0] += readcountListTmp[i]
                    sampletRFEntityDic[sample][mergedName][1] += RP100KListTmp[i]
                else:
                    sampletRFEntitySummaryDic[sample][0] += readcountListTmp[i]
            line = inf.readline()
    discard_fname = Path(workDir)/"discarded.reads.summary.assigningtRFs.csv"
    with open(discard_fname, 'w') as outf:
        outf.write('sample name,percentage of discarded reads,details\n')
        for sample in sampleList:
            try:
                outf.write(sample+',%.2f%%,%d\\%d\n'%(round((float(sampletRFEntitySummaryDic[sample][0])/sampletRFEntitySummaryDic[sample][1])*100.0, 2), sampletRFEntitySummaryDic[sample][0], sampletRFEntitySummaryDic[sample][1]))
            except ZeroDivisionError:
                outf.write(sample+',0.00%%,%d\\%d\n'%(sampletRFEntitySummaryDic[sample][0], sampletRFEntitySummaryDic[sample][1]))
    outf3_1 = open((Path(workDir)/'tRF.Counts.csv'), 'w')
    outf3_2 = open((Path(workDir)/'tRF.RP100K.csv'), 'w')
    outf3_1.write('entry name,'+','.join(sampleList)+'\n')
    outf3_2.write('entry name,'+','.join(sampleList)+'\n')
    for trfMerged in trfMergedList:
        outf3_1.write(trfMerged+','+','.join([str(sampletRFEntityDic[sample][trfMerged][0]) for sample in sampleList])+'\n')
        outf3_2.write(trfMerged+','+','.join(['%.2f'%(round(sampletRFEntityDic[sample][trfMerged][1], 2)) for sample in sampleList])+'\n')
    outf3_1.close()
    outf3_2.close()

    # generate the detailed potential tRFs for each sample
    tRF_dir = Path(workDir)/'tRFs.samples.tmp'
    os.system('mkdir %s'%(tRF_dir))
    sampletRFDic = {}
    for sample in sampleList:
        sampletRFDic.update({sample:{}})
    for seqKey in trfContentDic.keys():
        for i in range(len(sampleList)):
            if trfContentDic[seqKey]['count'][i] > 0:
                for subkey in trfContentDic[seqKey].keys():
                    if subkey not in ['uid', 'RPM', 'count']:
                        if 'pre' not in subkey:
                            templateSeq = trnaStruDic[subkey]['seq']
                        else:
                            templateSeq = pretrnaNameSeqDic[subkey]
                        filledSeq = trfContentDic[seqKey][subkey]['start']*'-'+seqKey+(len(templateSeq)-trfContentDic[seqKey][subkey]['start']-len(seqKey))*'-'
                        try:
                            sampletRFDic[sampleList[i]][subkey].append((trfContentDic[seqKey]['count'][i], trfContentDic[seqKey][subkey]['start'], seqKey, filledSeq, trfContentDic[seqKey][subkey]['tRFType'], trfContentDic[seqKey]['RPM'][i]))
                        except KeyError:
                            sampletRFDic[sampleList[i]][subkey] = [(trfContentDic[seqKey]['count'][i], trfContentDic[seqKey][subkey]['start'], seqKey, filledSeq, trfContentDic[seqKey][subkey]['tRFType'], trfContentDic[seqKey]['RPM'][i])]
    for sample in sampleList:
        outf_tmp = open((Path(tRF_dir)/(sample+'.aligned_tRFs.summary.report')), 'w')
        outf_tmp.write('amino acid\tCounts\tRP100K\tUnique reads\n')
        aaTypeList = []
        aaTypeDict = {}
        with open((Path(tRF_dir)/(sample+'.aligned_tRFs.report')), 'w') as outf:
            sumtRNAList = []
            for tRNAName in sampletRFDic[sample].keys():
                readSum = sum([trfSet[0] for trfSet in sampletRFDic[sample][tRNAName]])
                RPMSum = sum([trfSet[5] for trfSet in sampletRFDic[sample][tRNAName]])
                sumtRNAList.append((readSum, tRNAName, RPMSum))
            sumtRNAList.sort(reverse=True)
            for (readSum, tRNAName, RPMSum) in sumtRNAList:
                aaType = trnaAAanticodonDic[tRNAName]['aaType']
                if 'pre_' in tRNAName:
                    aaType = 'pre:'+aaType
                aaAnticodonPair = aaType+'-'+trnaAAanticodonDic[tRNAName]['anticodon']
                if 'pre:' in aaType:
                    if aaType+' tRF-1' not in aaTypeList:
                        aaTypeList.append(aaType+' tRF-1')
                        aaTypeDict.update({aaType+' tRF-1':[0,0,0]})
                else:
                    if aaType+" 5'" not in aaTypeList:
                        aaTypeList.append(aaType+" 5'")
                        aaTypeDict.update({aaType+" 5'":[0,0,0]})
                    if aaType+" 3'" not in aaTypeList:
                        aaTypeList.append(aaType+" 3'")
                        aaTypeDict.update({aaType+" 3'":[0,0,0]})
                    if aaType+" other" not in aaTypeList:
                        aaTypeList.append(aaType+" other")
                        aaTypeDict.update({aaType+" other":[0,0,0]})
                sampletRFDic[sample][tRNAName].sort(reverse=True)
                outf.write(tRNAName+'\t'+'read count sum:'+str(readSum)+'\tRP100K sum:'+'%.3f'%(round(RPMSum, 3))+'\n')
                if 'pre' not in tRNAName:
                    templateSeq = trnaStruDic[tRNAName]['seq']
                    templateType = 'mature tRNA'
                else:
                    templateSeq = pretrnaNameSeqDic[tRNAName]
                    templateType = 'primary tRNA trailer'
                for trfSet in sampletRFDic[sample][tRNAName]:
                    if len(trfSet[3]) == len(templateSeq):
                        outf.write(trfSet[3]+'\t'+trfSet[4]+'\t'+str(trfSet[0])+'\t'+'%.3f'%(round(trfSet[5], 3))+'\n')
                        if 'pre:' in aaType:
                            aaTypeDict[aaType+' tRF-1'][0] = aaTypeDict[aaType+' tRF-1'][0] + trfSet[0]
                            aaTypeDict[aaType+' tRF-1'][1] = aaTypeDict[aaType+' tRF-1'][1] + trfSet[5]
                            aaTypeDict[aaType+' tRF-1'][2] = aaTypeDict[aaType+' tRF-1'][2] + 1
                        else:
                            leftCount, rightCount = locationOfStartEnd(trfSet[3])
                            if leftCount <= 2:
                                aaTypeDict[aaType+" 5'"][0] = aaTypeDict[aaType+" 5'"][0] + trfSet[0]
                                aaTypeDict[aaType+" 5'"][1] = aaTypeDict[aaType+" 5'"][1] + trfSet[5]
                                aaTypeDict[aaType+" 5'"][2] = aaTypeDict[aaType+" 5'"][2] + 1
                            else:
                                if rightCount <= 2:
                                    aaTypeDict[aaType+" 3'"][0] = aaTypeDict[aaType+" 3'"][0] + trfSet[0]
                                    aaTypeDict[aaType+" 3'"][1] = aaTypeDict[aaType+" 3'"][1] + trfSet[5]
                                    aaTypeDict[aaType+" 3'"][2] = aaTypeDict[aaType+" 3'"][2] + 1
                                else:
                                    aaTypeDict[aaType+" other"][0] = aaTypeDict[aaType+" other"][0] + trfSet[0]
                                    aaTypeDict[aaType+" other"][1] = aaTypeDict[aaType+" other"][1] + trfSet[5]
                                    aaTypeDict[aaType+" other"][2] = aaTypeDict[aaType+" other"][2] + 1
                outf.write(templateSeq+'\t'+templateType+'\t'+str(readSum)+'\t'+'%.3f'%(round(RPMSum, 3))+'\n')
        for aaType in aaTypeList:
            outf_tmp.write('\t'.join([aaType, str(aaTypeDict[aaType][0]), '%.3f'%(round(aaTypeDict[aaType][1], 3)), str(aaTypeDict[aaType][2])]))
            outf_tmp.write('\n')
        outf_tmp.close()
        # cluter the reads that are alligned to the specific tRNA.
        readInforList = load_data_new((Path(tRF_dir)/(sample+'.aligned_tRFs.report')))
        dcValue, rhomin, deltamin, RP100KCutoff = 3.0, 5.0, 8.0, 10.0
        tRNANameList = []
        tRNANameDic = {}
        # in Unified_tRNANameList, 'pre_Homo_sapiens_tRNA-Cys-GCA-2-4_trailer' will be treated as 'Homo_sapiens_tRNA-Cys-GCA-2-4'
        Unified_tRNANameList = []
        tRNANameList_selected = []
        Unified_tRNANameList_selected = []
        potential_tRFs_rep1 = sample+'.aligned_tRFs.report'
        with open(Path(tRF_dir)/potential_tRFs_rep1, 'r') as inf:
            for line in inf:
                if 'RP100K sum:' in line:
                    tRNAName = line.split('\t')[0]
                    tRNANameList.append(tRNAName)
                    tRNANameDic.update({tRNAName:[]})
                    if 'pre' in tRNAName:
                        tRNAName = '_'.join(tRNAName.split('_')[1:-1])
                    if tRNAName not in Unified_tRNANameList:
                        Unified_tRNANameList.append(tRNAName)
        potential_tRFs_file = sample+'.aligned_tRFs.clusters.detail'
        tRFs_report_file = sample+'.tRFs.report.tsv'
        outf_1 = open(Path(tRF_dir)/potential_tRFs_file, 'w')
        outf_2 = open(Path(tRF_dir)/tRFs_report_file, 'w')
        for indexId, readInforDic in enumerate(readInforList):
            tRNANameTmp = ""
            tRNANameTmp = tRNANameList[indexId]
            distanceDic, max_dis, min_dis, max_id = getDistance(readInforDic)
            rho = local_density(distanceDic, readInforDic, max_id, dc=dcValue)
            delta, nneigh, nneighDistance, sort_rho_idx = min_distance(distanceDic, max_dis, max_id, rho)
            # assign cluster center
            NCLUST = 0
            cl = np.zeros(max_id+1)-1
            ccenter = {}
            for idx, (ldensity, mdistance, nneigh_item) in enumerate(zip(rho, delta, nneigh)):
                if idx == 0:
                    continue
                if ldensity >= rhomin and mdistance >= deltamin:
                    NCLUST = NCLUST + 1
                    cl[idx] = NCLUST
                    ccenter[NCLUST] = idx
            # special situation: all of the reads are in one cluster.
            if NCLUST == 0:
                ldensityList = [-1.0]
                mdistanceList = [-1.0]
                for idx, (ldensity, mdistance, nneigh_item) in enumerate(zip(rho, delta, nneigh)):
                    if idx == 0:
                        continue
                    else:
                        ldensityList.append(ldensity)
                        mdistanceList.append(mdistance)
                if max(mdistanceList) <= deltamin and  max(ldensityList) >= rhomin:
                    NCLUST = 1
                    idx = ldensityList.index(max(ldensityList))
                    cl[idx] = NCLUST
                    ccenter[NCLUST] = idx
            # assignation
            for i in range(max_id):
                    if cl[sort_rho_idx[i]]==-1:
                        cl[sort_rho_idx[i]] = cl[nneigh[sort_rho_idx[i]]]
            # transform cl into int type:
            cl = cl.astype(np.int32)
            # assign cluster core and cluster halo for the points in each cluster
            halo = np.zeros(max_id+1)
            halo[:] = cl
            if NCLUST>1:
                bord_rho = np.zeros(NCLUST+1)
                # calclate the the average of rho in the border region of each cluster
                for i in range(1, max_id):
                    for j in range(i+1, max_id+1):
                        if cl[i]!=cl[j] and distanceDic[(i,j)] <= dcValue:
                            rho_aver = (rho[i]+rho[j])/2
                            if rho_aver > bord_rho[cl[i]]:
                                bord_rho[cl[i]] = rho_aver
                            if rho_aver > bord_rho[cl[j]]:
                                bord_rho[cl[j]] = rho_aver
                # if halo[i] of point i is 0, this point is a cluster halo (outlier)
                for i in range(1, max_id+1):
                    if rho[i]<bord_rho[cl[i]]:
                        halo[i] = 0
                    # if the distance of the point to its cluster center is larger than the deltamin, it's a cluster halo as well.
                    try:
                        if distanceDic[(i,ccenter[cl[i]])] > deltamin:
                            halo[i] = 0
                    except KeyError:
                        continue
            elif NCLUST==1:
                for i in range(1, max_id+1):
                    try:
                        if distanceDic[(i,ccenter[cl[i]])] > deltamin:
                            halo[i] = 0
                    except KeyError:
                        continue
            else:
                halo[:] = np.zeros(max_id+1)
            paraClusterCon = ''
            paraClusterCon += tRNANameTmp+':\n'
            
            sumCount = 0
            sumRP100K = 0.0
            coreCountList  = []
            haloCountList = []
            coreRP100KList = []
            haloRP100KList = []
            ClusterCount = 1
            clusterContentList = []
            if NCLUST >= 1:
                for i in range(1, NCLUST+1):
                    nc = 0
                    nh = 0
                    center_idx = ccenter[i]
                    select_idx = []
                    totalHaloCount = 0
                    totalHaloRP100K = 0.0
                    for j in range(1, max_id+1):
                        if cl[j] == i:
                            nc = nc+1
                        if halo[j] == i:
                            nh = nh+1
                            select_idx.append(j)
                        if cl[j] == i and halo[j] != i:
                            totalHaloCount = totalHaloCount + readInforDic[j]['count']
                            totalHaloRP100K = totalHaloRP100K + readInforDic[j]['RPM']
                    totalCoreCount = 0
                    totalCoreRP100K = 0.0
                    for idx in select_idx:
                        totalCoreCount = totalCoreCount + readInforDic[idx]['count']
                        totalCoreRP100K = totalCoreRP100K + readInforDic[idx]['RPM']
                    
                    if totalCoreCount > 0:
                        abundantSeqCount = 0
                        abundantSeqRPM = 0.0
                        # The abundant seq should be from the center point or the sencond point
                        tmpPair = []
                        tmpPair.append((readInforDic[center_idx]['count'], readInforDic[center_idx]['allignedSeq'], readInforDic[center_idx]['type']))
                        for idx in select_idx:
                            if idx != center_idx:
                                tmpPair.append((readInforDic[idx]['count'], readInforDic[center_idx]['allignedSeq'], readInforDic[center_idx]['type']))
                        tmpPair.sort(reverse=True)
                        abundantSeq = removeDash(tmpPair[0][1])
                        abundantSeqType = tmpPair[0][2]

                        headDashCount1, tailDashCount1 = DashCount(tmpPair[0][1])
                        pos_Infor = ':'.join([str(headDashCount1+1), str(len(tmpPair[0][1])-tailDashCount1)])
                        
                        paraClusterCon += 'Cluster: %d Total Read Count in Core: %d Total Read Count in Halo: %d Total RP100K in Core: %.2f Total RP100K in Halo: %.2f Center Index: %d Elements: %d Core: %d Halo: %d\n'%(ClusterCount, totalCoreCount, totalHaloCount, totalCoreRP100K, totalHaloRP100K, center_idx, nc, nh, nc-nh)
                        paraClusterCon += 'Center:\n'
                        paraClusterCon += '%s\t%s\t%d\t%.2f\n'%(readInforDic[center_idx]['allignedSeq'], readInforDic[center_idx]['type'], readInforDic[center_idx]['count'], readInforDic[center_idx]['RPM'])
                        abundantSeqCount = abundantSeqCount + readInforDic[center_idx]['count']
                        abundantSeqRPM = abundantSeqRPM + readInforDic[center_idx]['RPM']
                        for idx in select_idx:
                            if idx != center_idx:
                                paraClusterCon += '%s\t%s\t%d\t%.2f\n'%(readInforDic[idx]['allignedSeq'], readInforDic[idx]['type'], readInforDic[idx]['count'], readInforDic[idx]['RPM'])
                                abundantSeqCount = abundantSeqCount + readInforDic[idx]['count']
                                abundantSeqRPM = abundantSeqRPM + readInforDic[idx]['RPM']
                        paraClusterCon += '**********************************\n'
                        clusterContentList.append((abundantSeq, abundantSeqType, pos_Infor, abundantSeqCount, abundantSeqRPM))
                        ClusterCount = ClusterCount + 1
                    sumCount = sumCount + totalHaloCount
                    sumCount = sumCount + totalCoreCount
                    sumRP100K = sumRP100K + totalHaloRP100K
                    sumRP100K = sumRP100K + totalCoreRP100K
                    coreCountList.append(totalCoreCount)
                    haloCountList.append(totalHaloCount)
                    coreRP100KList.append(totalCoreRP100K)
                    haloRP100KList.append(totalHaloRP100K)
            else:
                for i in range(1, max_id+1):
                    sumCount = sumCount + readInforDic[i]['count']
                    sumRP100K = sumRP100K + readInforDic[i]['RPM']
            paraClusterCon += 'Summary:\nNumber of Clusters: %d\n'%(ClusterCount-1)
            paraClusterCon += 'total Read Count : %d\n'%(sumCount)
            paraClusterCon += 'total RP100K: %.2f\n'%(sumRP100K)
            paraClusterCon += 'total Cluster Core Read Count: %s=%d\n'%('+'.join([str(item) for item in coreCountList]), sum(coreCountList))
            paraClusterCon += 'total Cluster Core RP100K: %s=%.3f\n'%('+'.join([str(item) for item in coreRP100KList]), sum(coreRP100KList))
            paraClusterCon += 'total Cluster Halo Read Count: %s=%d\n'%('+'.join([str(item) for item in haloCountList]), sum(haloCountList))
            paraClusterCon += 'total Cluster Halo RP100K: %s=%.3f\n'%('+'.join([str(item) for item in haloRP100KList]), sum(haloRP100KList))
            paraClusterCon += '##################################\n'
            outf_1.write(paraClusterCon)
            clusterInfor = (ClusterCount-1, float(sum(coreCountList))/sumCount, sumCount, sumRP100K)
            tRNANameDic[tRNANameTmp]= clusterContentList
            # Only keep the tRNAs if it's RP100K is no less than RP100KCutoff
            if sumRP100K >= RP100KCutoff:
                tRNANameList_selected.append(tRNANameTmp)
                if 'pre' in tRNANameTmp:
                    tRNANameTmp = '_'.join(tRNANameTmp.split('_')[1:-1])
                if tRNANameTmp not in Unified_tRNANameList_selected:
                    Unified_tRNANameList_selected.append(tRNANameTmp)
        outf_1.close()
        outf_2.write('tRNA name\ttRNA sequence\ttRF sequence\ttRF mismatch\ttRF type\ttRF coordinate\tRead count\tRP100K\n')
        for tRNA_selected in Unified_tRNANameList_selected:
            if tRNA_selected in tRNANameList_selected:
                tRNA_seq = trnaStruDic[tRNA_selected]['seq']
                for contentTmp in tRNANameDic[tRNA_selected]:
                    mimatchState, mismatchPosition = detectMismach(contentTmp[0], tRNA_seq, contentTmp[2])
                    outf_2.write(tRNA_selected+'\t'+tRNA_seq+'\t'+contentTmp[0]+'\t'+':'.join([mimatchState, mismatchPosition])+'\t'+contentTmp[1]+'\t'+contentTmp[2]+'\t'+str(contentTmp[3])+'\t'+'%.2f'%(round(contentTmp[4], 2))+'\n')
            if 'pre_'+tRNA_selected+'_trailer' in tRNANameList_selected:
                tRNA_seq = pretrnaNameSeqDic['pre_'+tRNA_selected+'_trailer']
                for contentTmp in tRNANameDic['pre_'+tRNA_selected+'_trailer']:
                    newPos = ':'.join([contentTmp[2].split(':')[0], str(int(contentTmp[2].split(':')[1])-3)])
                    mimatchState, mismatchPosition = detectMismach(contentTmp[0][:-3], tRNA_seq, newPos)
                    outf_2.write(tRNA_selected+'\t'+tRNA_seq+'\t'+contentTmp[0]+'\t'+':'.join([mimatchState, mismatchPosition])+'\t'+contentTmp[1]+'\t'+contentTmp[2]+'\t'+str(contentTmp[3])+'\t'+'%.2f'%(round(contentTmp[4], 2))+'\n')
        outf_2.close()


def refineName(name):
    if '.fastq' in name:
        label = name.index('.fastq')
        newName = ''
        for i,item in enumerate(name):
            if i in range(label, label+len('.fastq')):
                pass
            else:
                newName = newName + item
    else:
        newName = name
    return newName


def checkSeqList(seqList, seqDic, onlyCanmiRNA):
    # Check whether there are at least one sequences belonging to known miRNA. 
    state = False
    for seq in seqList:
        seqNew = removeDash(seq)
        try:
            if seqNew in onlyCanmiRNA:
                state = True
                break
        except KeyError:
            print('seq %s and %s does not exit.'%(seqNew, seq))
    return state



def a2i_editing(args, cannonical, isomirs, base_names, workDir, Filtered_miRNA_Reads, mirMergedNameDic, mirDic, ref_db, seqDic, onlyCanmiRNA):
    a2IPercentageCutoff = 2.0
    startBase = 'A'
    endBase = 'G'
    mirNameSeqDicTmp = {}
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
    for can_list_element in canonical_gff:
        if can_list_element[1] in mirMergedNameDic:
            miRName = mirMergedNameDic[can_list_element[1]]
        else:
            miRName = can_list_element[1]
        try:
            mirNameSeqDicTmp[miRName].append(can_list_element[0])
        except KeyError:
            mirNameSeqDicTmp.update({miRName:[can_list_element[0]]})

    for iso_list_element in isomir_gff:
        if iso_list_element[1] in mirMergedNameDic:
            miRName = mirMergedNameDic[iso_list_element[1]]
        else:
            miRName = iso_list_element[1]
        iso_sample_values = iso_list_element[2:]
        if any([(iso_sample_values[i]*freq_list[i]) >=1 for i in range(len(base_names))]):
            try:
                mirNameSeqDicTmp[miRName].append(iso_list_element[0])
            except KeyError:
                mirNameSeqDicTmp.update({miRName:[iso_list_element[0]]})
        else:
            pass
    samToMapFasta = Path(workDir)/"SeqToMap.fasta"
    with open(samToMapFasta,'w') as outf:
        for mirName in mirNameSeqDicTmp.keys():
            for seqTmp in mirNameSeqDicTmp[mirName]:
                outf.write('>'+seqTmp+'\n'+seqTmp+'\n')

    ###ARUN
    mirNameSeqDic={}
    #human_mirna_SNP_pseudo_miRBase.fa mirNameSeqDic 
    fastaName = args.organism_name+"_mirna_SNP_pseudo_"+ref_db+".fa"
    mirNameSeqDic_fasta = Path(args.libraries_path)/args.organism_name/"fasta.Libs"/fastaName
    
    
    #Getting list of miRNAs which overlap with repeat elements
    repeatElementName = args.organism_name+'_miRNAs_in_repetitive_element_'+ref_db+".csv"
    miRNAs_in_repetitive_element = Path(args.libraries_path)/args.organism_name/'annotation.Libs'/repeatElementName
    removedMiRNAList = []
    try:
        with open(miRNAs_in_repetitive_element, 'r') as inf:
            for line in inf:
                if line.strip().split(',')[0] not in removedMiRNAList:
                    removedMiRNAList.append(line.strip().split(',')[0])
    except FileNotFoundError:
        pass
    

    with open(mirNameSeqDic_fasta,'r') as inf:
        line = inf.readline()
        while line != '':
            mirName = line.strip()[1:]
            line = inf.readline()
            sequence = line.strip()
            mirNameSeqDic.update({mirName:sequence})
            line = inf.readline()
    
    bwtCommand = str(Path(args.bowtie_path)/"bowtie ") if args.bowtie_path else "bowtie "
    if args.bowtieVersion == "False": # That is if version is v1.3.0
        bwtCommand = str(bwtCommand) + " --threads " + str(args.threads) + " -x "
    else:
        bwtCommand = str(bwtCommand) + " --threads " + str(args.threads) + " "

    
    indexName  = str(args.organism_name) + str("_genome") 
    genome_index = Path(args.libraries_path)/args.organism_name/"index.Libs"/indexName
    bwtCommand = str(bwtCommand) + str(genome_index) + ' -n 1 -f -a -3 2 ' + str(samToMapFasta)
    retainedSeqDic = {}
    retainedSeqContentDicTmp = {}
    if args.phred64:
        bwtCommand = str(bwtCommand) + ' --phred64-quals '
    bowtie = subprocess.run(str(bwtCommand), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
    if bowtie.returncode==0:
        bwtOut = bowtie.stdout
        bwtErr = bowtie.stderr
    for srow in bwtOut.split('\n'):
        if not srow.startswith('@'):
            sam_line = srow.split('\t')
            if sam_line != ['']:
                try:
                    retainedSeqContentDicTmp[sam_line[0]].append(sam_line[-1].count(':'))
                except KeyError:
                    retainedSeqContentDicTmp.update({sam_line[0]:[sam_line[-1].count(':')]})
                #['AAAAACTGAGACTACTTTTG', '+', 'chr10', '110988977', 'AAAAACTGAGACTACTTT', 'IIIIIIIIIIIIIIIIII', '0', '']
                #['AAAAACTGAGACTACTTTTG', '-', 'chr2', '236950434', 'AAAGTAGTCTCAGTTTTT', 'IIIIIIIIIIIIIIIIII', '0', '']
                #['AAAAACTGAGACTACTTTTG', '-', 'chr14', '65022420', 'AAAGTAGTCTCAGTTTTT', 'IIIIIIIIIIIIIIIIII', '0', '14:T>G']
                #['AAAAACTGAGACTACTTTTG', '-', 'chr4', '134615933', 'AAAGTAGTCTCAGTTTTT', 'IIIIIIIIIIIIIIIIII', '0', '14:A>G']
    for seq_tmp in retainedSeqContentDicTmp.keys():
        kept = False
        content = retainedSeqContentDicTmp[seq_tmp]
        if len(content) == 1:
            kept = True
        else:
            minimum = min(content)
            if len([subitem for subitem in content if subitem == minimum]) == 1:
                kept = True
        if kept:
            retainedSeqDic.update({seq_tmp:True})
    
    mirNamePositionValueDict = {}
    mirNamePositionList = []
    a2IEditingFileTmp1 = Path(workDir)/'a2IEditing.report.tmp1.csv'
    a2IEditingFileTmp2 = Path(workDir)/'a2IEditing.report.tmp2.csv'
    a2IEditingFileTmp3 = Path(workDir)/'a2IEditing.report.tmp3.csv'
    a2IEditingFile = Path(workDir)/'a2IEditing.report.csv'
    a2IEditingFileTrans = Path(workDir)/'a2IEditing.report.newform.csv'
    a2IEditingDetailFile = Path(workDir)/'a2IEditing.detail.txt'
    mismatchCountFile = Path(workDir)/'mismatchCount.csv'

# (args, cannonical, isomirs, base_names, workDir, Filtered_miRNA_Reads, mirMergedNameDic):
    readCounts = list(Filtered_miRNA_Reads.values()) # YOU MAY NEED TO DELETE THIS VARIABLE AS I FORGOT TO USE IT
    sampleList = base_names
    outfTmp = open(a2IEditingDetailFile, 'w')
    for mirName in mirNameSeqDicTmp.keys():
        targetSeq = mirNameSeqDic[mirName]
        Filtered_miRNA_Reads_list = list(Filtered_miRNA_Reads.values())
        #human_mirna_SNP_pseudo_miRBase.fa mirNameSeqDic 
        #mirDic miRCounts 
        seqList = mirNameSeqDicTmp[mirName]
        for i in range(len(base_names)):
            try:
                if 1000000.0*mirDic[mirName][i]/Filtered_miRNA_Reads_list[i] >= 1:
                    countListTmp = []
                    seqListTmp = []
                    for seqTmp in seqList:
                        if 1000000.0*seqDic[seqTmp][i]/Filtered_miRNA_Reads_list[i] >= 1 or (seqDic[seqTmp][i] > 0):
                            seqListTmp.append(seqTmp)
                            countListTmp.append(seqDic[seqTmp][i])
                    if len(countListTmp) >1:
                        if min(countListTmp) > 0:
                            alignSeqListKept, positionList, a2IPositionCountDic, a2IPositionRatioDic, a2IPositionPvalueDic, countTrue, seqCountTrue, canonicalSeqCount = A2IEditing(targetSeq, seqListTmp, countListTmp, mirName, outfTmp, retainedSeqDic, startBase, endBase)
                            for position in positionList:
                                mirNamePosition = ':'.join([mirName, str(position)])
                                if mirNamePosition not in mirNamePositionValueDict.keys():
                                    mirNamePositionList.append(mirNamePosition)
                                    mirNamePositionValueDict.update({mirNamePosition:[[] for t in range(len(sampleList))]})
                                    mirNamePositionValueDict[mirNamePosition][i] = mirNamePositionValueDict[mirNamePosition][i] + [alignSeqListKept, str(sum(countListTmp)), str(len([item1 for item1 in countListTmp if item1 >0])), str(countTrue), str(seqCountTrue), str(canonicalSeqCount), str(a2IPositionCountDic[position]), a2IPositionRatioDic[position], a2IPositionPvalueDic[position]]
                                else:
                                    mirNamePositionValueDict[mirNamePosition][i] = mirNamePositionValueDict[mirNamePosition][i] + [alignSeqListKept, str(sum(countListTmp)), str(len([item1 for item1 in countListTmp if item1 >0])), str(countTrue), str(seqCountTrue), str(canonicalSeqCount), str(a2IPositionCountDic[position]), a2IPositionRatioDic[position], a2IPositionPvalueDic[position]]
            except ZeroDivisionError:
                pass
    outfTmp.close()

    basePairMismachCountList_total = []
    for i in range(len(base_names)):
        basePairMismachCountListTmp = [[0, 0, 0] for j in range(12)]
        for mirName in mirNameSeqDicTmp.keys():
            targetSeq = mirNameSeqDic[mirName]
            seqList = mirNameSeqDicTmp[mirName]
            try:
                if 1000000.0*mirDic[mirName][i]/Filtered_miRNA_Reads_list[i] >= 1:
                    countListTmp = []
                    seqListTmp = []
                    for seqTmp in seqList:
                        if 1000000.0*seqDic[seqTmp][i]/Filtered_miRNA_Reads_list[i] >= 1 or (seqDic[seqTmp][i] > 0):
                            seqListTmp.append(seqTmp)
                            countListTmp.append(seqDic[seqTmp][i])
                    if len(countListTmp) >1:
                        if min(countListTmp) > 0:
                            basePairMismachCountList = mismatchCountAnalysis(targetSeq, seqListTmp, countListTmp, retainedSeqDic)
                            for indexTmp, s in enumerate(basePairMismachCountList):
                                basePairMismachCountListTmp[indexTmp][0] = basePairMismachCountListTmp[indexTmp][0] + s[0]
                                basePairMismachCountListTmp[indexTmp][1] = basePairMismachCountListTmp[indexTmp][1] + s[1]
                                basePairMismachCountListTmp[indexTmp][2] = basePairMismachCountListTmp[indexTmp][2] + s[2]
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
            except ZeroDivisionError:
                pass
        basePairMismachCountList_total.append(basePairMismachCountListTmp)
    #print(basePairMismachCountList_total)

        # P-value required after multiple testing correction. Here, Benjamini-Hochberg correction is performed.
    for i in range(len(sampleList)):
        p_value_list = []
        for mirNamePosition in mirNamePositionList:
            if len(mirNamePositionValueDict[mirNamePosition][i]) != 0:
                p_value_list.append([mirNamePositionValueDict[mirNamePosition][i][8], mirNamePosition])
        p_value_list.sort()
        for index_tmp, p_value in enumerate(p_value_list):
            adjusted_p_value = p_value[0]*len(p_value_list)/(index_tmp+1)
            mirNamePositionValueDict[p_value[1]][i].append(adjusted_p_value)

    with open(a2IEditingFileTmp1, 'w') as outf:
        outf.write('miRNA,A-to-I position in the miRNA,miRNA sequence')
        for i in range(len(sampleList)):
            name_new = refineName(sampleList[i])
            outf.write(','+','.join([name_new+'.readCount', name_new+'.readCount.canonical', name_new+'.RPM.canonical', name_new+'.readCount.mismatch', name_new+'.RPM.mismatch', name_new+'.AtoI.percentage', name_new+'.AtoI.adjusted.pValue']))
        outf.write('\n')
        for mirNamePosition in mirNamePositionList:
            mirName = mirNamePosition.split(':')[0].strip()
            position = mirNamePosition.split(':')[1].strip()
            targetSeq = mirNameSeqDic[mirName]
            a2Icontent = mirNamePositionValueDict[mirNamePosition]
            outf.write(','.join([mirName, position, targetSeq]))
            for i in range(len(sampleList)):
                outf.write(',')
                if len(a2Icontent[i]) != 0:
                    alignSeqKept = a2Icontent[i][0]
                    try:
                        if checkSeqList(alignSeqKept, seqDic, onlyCanmiRNA):
                            outf.write('%s'%(a2Icontent[i][3]))
                            outf.write(',%s'%(a2Icontent[i][5]))
                            outf.write(',%.2f'%(1000000.0*int(a2Icontent[i][5])/Filtered_miRNA_Reads_list[i]))
                            outf.write(',%s'%(a2Icontent[i][6]))
                            outf.write(',%.2f'%(1000000.0*int(a2Icontent[i][6])/Filtered_miRNA_Reads_list[i]))
                            outf.write(',%.2f%%'%(a2Icontent[i][7]*100))
                            if a2Icontent[i][9] <= 0.05:
                                outf.write(',%.2E'%(a2Icontent[i][9]))
                            else:
                                outf.write(',NS')
                        else:
                            #outf.write(','.join(a2Icontent[i][1:6]))
                            #outf.write(',NE,NE,NE')
                            outf.write('%s'%(a2Icontent[i][3]))
                            outf.write(',%s'%(a2Icontent[i][5]))
                            outf.write(',%.2f'%(1000000.0*int(a2Icontent[i][5])/Filtered_miRNA_Reads_list[i]))
                            outf.write(',%s'%(a2Icontent[i][6]))
                            outf.write(',%.2f'%(1000000.0*int(a2Icontent[i][6])/Filtered_miRNA_Reads_list[i]))
                            outf.write(',NE,NE')
                    except KeyError:
                        print('Error happens at %s: %s'%(mirName, '\t'.join(alignSeqKept)))
                        outf.write('%s'%(a2Icontent[i][3]))
                        outf.write(',%s'%(a2Icontent[i][5]))
                        outf.write(',%.2f'%(1000000.0*int(a2Icontent[i][5])/Filtered_miRNA_Reads_list[i]))
                        outf.write(',%s'%(a2Icontent[i][6]))
                        outf.write(',%.2f'%(1000000.0*int(a2Icontent[i][6])/Filtered_miRNA_Reads_list[i]))
                        outf.write(',NE,NE')
                else:
                    outf.write(','.join(['NE','NE','NE','NE','NE','NE','NE']))
            outf.write('\n')
        
    # Remove some unimportant rows in *.a2IEditing.report.tmp1.csv
    outf = open(a2IEditingFileTmp2, 'w')
    with open(a2IEditingFileTmp1, 'r') as inf:
        line = inf.readline()
        outf.write(line)
        line = inf.readline()
        while line != '':
            content = line.strip().split(',')
            portion = ((len(content)-9)/7)+1
            keptState = False
            for i in range(int(portion)):
                try:
                    tmpVale = float(content[9+i*7])
                    keptState = True
                    break
                except ValueError:
                    pass
            if keptState:
                outf.write(line)
            line = inf.readline()
    outf.close()
    # Sort the report file based on the first column.
    os.system("(head -n 1 %s && tail -n +2 %s | sort -t',' -k1,1 -k2,2n) > %s"%(a2IEditingFileTmp2, a2IEditingFileTmp2, a2IEditingFileTmp3))
    # Remove the miRNAs if:
    # 1) their canonical sequences' RPM < 1 across all of the samples
    # 2) miRNAs are located at repetitive element region.
    # 3) the correponding 1 nucluetide mismathed sequence can be alligned to no less than 1 locations on the genome with trimming the last two nucluetides at 3' in the setting of Bowtie.
    miRNAPositionDic = {}
    seqtojudge=Path(workDir)/'SeqToJudge.fasta'
    with open(seqtojudge, 'w') as outf:
        with open(a2IEditingFileTmp3, 'r') as inf:
            line = inf.readline()
            line = inf.readline()
            while line != '':
                content = line.strip().split(',')
                miRNAPositionName = ':'.join([content[0], content[1]])
                mismathedSeq = ''
                for index, necleotide in enumerate(content[2]):
                    if index != int(content[1])-1:
                        mismathedSeq = mismathedSeq + necleotide
                    else:
                        mismathedSeq = mismathedSeq + endBase
                canonicalRPMList = []
                for i in range(int(((len(content)-9)/7)+1)):
                    try :
                        canonicalRPMList.append(float(content[5+i*7]))
                    except ValueError:
                        canonicalRPMList.append(0)
                canonicalRPMState = not any([item >= 1 for item in canonicalRPMList])
                dicTmp = {}
                dicTmp.update({'mismathedSeq':mismathedSeq})
                dicTmp.update({'RPMRemoveState':canonicalRPMState})
                dicTmp.update({'ReElementRemoveState':content[0] in removedMiRNAList})
                dicTmp.update({'SeqRemoveState':False})
                miRNAPositionDic.update({miRNAPositionName:dicTmp})
                outf.write('>'+mismathedSeq+'\n'+mismathedSeq+'\n')
                line = inf.readline()
    bwtCommand = str(Path(args.bowtie_path)/"bowtie ") if args.bowtie_path else "bowtie "
    if args.bowtieVersion == "False": # That is if version is v1.3.0
        bwtCommand = str(bwtCommand) + " --threads " + str(args.threads) + " -x "
    else:
        bwtCommand = str(bwtCommand) + " --threads " + str(args.threads) + " "

    indexName  = str(args.organism_name) + str("_genome") 
    genome_index = Path(args.libraries_path)/args.organism_name/"index.Libs"/indexName
    bwtCommand = str(bwtCommand) + str(genome_index) + ' -n 0 -f -a -3 2 ' + str(seqtojudge)
    if args.phred64:
        bwtCommand = bwtCommand + ' --phred64-quals '
    bowtie = subprocess.run(str(bwtCommand), shell=True, check=True, stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE, universal_newlines=True)
    if bowtie.returncode==0:
        bwtOut = bowtie.stdout
        bwtErr = bowtie.stderr
    removedSeqList = []
    for srow in bwtOut.split('\n'):
        if not srow.startswith('@'):
            sam_line = srow.split('\t')
            if sam_line != [''] and sam_line[0] not in removedSeqList:
                removedSeqList.append(sam_line[0])
    os.system('rm %s '%(seqtojudge))

    for key in miRNAPositionDic.keys():
        if miRNAPositionDic[key]['mismathedSeq'] in removedSeqList:
            miRNAPositionDic[key]['SeqRemoveState'] = True
            
    with open(a2IEditingFile, 'w') as outf:
        with open(a2IEditingFileTmp3, 'r') as inf:
            line = inf.readline()
            outf.write(line)
            line = inf.readline()
            while line != '':
                miRNAPositionName = ':'.join([line.strip().split(',')[0], line.strip().split(',')[1]])
                if miRNAPositionDic[miRNAPositionName]['RPMRemoveState'] or miRNAPositionDic[miRNAPositionName]['ReElementRemoveState'] or miRNAPositionDic[miRNAPositionName]['SeqRemoveState']:
                    pass
                else:
                    outf.write(line)
                line = inf.readline()
    # Remove the *.a2IEditing.report.tmp1.csv *.a2IEditing.report.tmp2.csv *.a2IEditing.report.tmp3.csv
    os.system('rm %s %s %s'%(a2IEditingFileTmp1, a2IEditingFileTmp2, a2IEditingFileTmp3))
    with open(mismatchCountFile, 'w') as outf:
        f1 = ['>'.join(item)+'_raw' for item in [('A', 'G'),('A', 'C'),('A', 'T'),('T', 'G'),('T', 'A'),('T', 'C'),('C', 'G'),('C', 'A'),('C', 'T'),('G', 'A'),('G', 'C'),('G', 'T')]]
        f2 = ['>'.join(item)+'' for item in [('A', 'G'),('A', 'C'),('A', 'T'),('T', 'G'),('T', 'A'),('T', 'C'),('C', 'G'),('C', 'A'),('C', 'T'),('G', 'A'),('G', 'C'),('G', 'T')]]
        f3 = ['>'.join(item)+'_filtered' for item in [('A', 'G'),('A', 'C'),('A', 'T'),('T', 'G'),('T', 'A'),('T', 'C'),('C', 'G'),('C', 'A'),('C', 'T'),('G', 'A'),('G', 'C'),('G', 'T')]]
        outf.write('sample,'+','.join(f1)+','+','.join(f2)+','+','.join(f3)+'\n')
        for i in range(len(sampleList)):
            outf.write(sampleList[i]+',')
            s1 = [str(item[0]) for item in basePairMismachCountList_total[i]]
            s2 = [str(item[1]) for item in basePairMismachCountList_total[i]]
            s3 = [str(item[2]) for item in basePairMismachCountList_total[i]]
            outf.write(','.join(s1))
            outf.write(',')
            outf.write(','.join(s2))
            outf.write(',')
            outf.write(','.join(s3))
            outf.write('\n')
    # Remove the mismatchCountFile 
    os.system('rm %s'%(mismatchCountFile))
    

    # Transform the format of a2IEditing.report.csv in order to plot heatmap.
    sampleDicTmp = {}
    sampleDicTmp.update({'SRR837842':'Colon 1'})
    sampleDicTmp.update({'SRR837839':'Colon 2'})
    sampleDicTmp.update({'SRR5127219':'Colon cell'})
    sampleDicTmp.update({'SRR1646473':'Colon cancer 1'})
    sampleDicTmp.update({'SRR1646493':'Colon cancer 2'})
    sampleDicTmp.update({'SRR1917324':'DKO1'})
    sampleDicTmp.update({'SRR1917336':'DLD1'})
    sampleDicTmp.update({'SRR1917329':'DKS8'})
    sampleDicTmp.update({'SRR567638':'Placenta 2'})

    with open(a2IEditingFileTrans, 'w') as outf:
        outf.write('miRNA:position,sample,A-to-I percentage,log2RPM\n')
        sampleListTmp = []
        miRNAListTmp = []
        miRNAContentListTmp = []
        miRNAContentDicTmp = {}
        with open(a2IEditingFile, "r") as inf:
            line = inf.readline()
            for item in line.strip().split(','):
                if '.AtoI.percentage' in item and item.split('.')[0] not in sampleListTmp:
                    sampleListTmp.append(item.split('.')[0] )
            line = inf.readline()
            while line != '':
                content = line.strip().split(',')
                miRNANameNew = ":".join(content[:2])
                tmp = content[3:]
                contentTmpList = []
                for item in [tmp[i:i + 7] for i in range(0, len(tmp), 7)]:
                    if item[6] == 'NE' or item[6] == 'NS':
                        contentTmp = ('NA', 'NA')
                    else:
                        contentTmp = (item[5][:-1], str(math.log(float(item[4]),2)))
                    contentTmpList.append(contentTmp)
                # The miRNA will be kept if the percentage of at least one sample is larger than 5.0%
                state = False
                for item in contentTmpList:
                    try:
                        if float(item[0]) >= a2IPercentageCutoff:
                            state = True
                            break
                    except ValueError:
                        pass
                if state:
                    miRNAListTmp.append(miRNANameNew)
                    miRNAContentListTmp.append(contentTmpList)
                    miRNAContentDicTmp.update({miRNANameNew:contentTmpList})
                line = inf.readline()
        miRNAList2 = []
        for index, item in enumerate(miRNAListTmp):
            count = 0
            for item2 in miRNAContentListTmp[index]:
                if item2[0] != 'NA':
                    count = count + 1
            miRNAList2.append((count, item))
        miRNAList2.sort(reverse=True)

        for index1, sample in enumerate(sampleListTmp):
            for tmp in miRNAList2:
                miRNA = tmp[1]
                if sample in sampleDicTmp.keys():
                    outf.write(miRNA+','+sampleDicTmp[sample]+','+','.join(miRNAContentDicTmp[miRNA][index1])+'\n')
                else:
                    outf.write(miRNA+','+sample+','+','.join(miRNAContentDicTmp[miRNA][index1])+'\n')

    if len(miRNAList2) <= 5:
        print('The number of A-to-I editing sites is less than 10 so that no heatmap is drawn.')
    else:
        RscriptDirTmp = Path(__file__).resolve().parents[1]
        RscriptDir = Path(RscriptDirTmp)/('rScripts')/('A-to-I_plot.R')
        outA2Ipdf = Path(workDir)/('a-to-I.heatmap.pdf')
        os.system('Rscript %s %s %s'%(RscriptDir, a2IEditingFileTrans, outA2Ipdf))
    
    os.remove(samToMapFasta)
