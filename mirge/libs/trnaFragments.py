#!/usr/bin/env python
import time
import math
import random
import re
from pathlib import Path
import numpy as np
import os, sys

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
        if target_seq[i] != seqTmp[i]:
            mimatchState = 'Y'
            mismatchPos.append(str(startPos+1+i))
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
    
    trfFile = Path(workDir)/"tRFs.potential.report.tsv"
    with open(trfFile, 'w') as outf:
        outf.write('read sequence\tuid\tread count(%s)\tRP100K (%s)\tamino acid all hits\tamino acid-anticodon all hits\ttRF information all hits\tamino acid all deduplicated hits\tamino acid-anticodon all deduplicated hits\ttRF information all deduplicated hits\tamino acid one hit\tamino acid-anticodon one hit\ttRF information one hit\n'%(','.join(sampleList), ','.join(sampleList)))
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
                outf.write('\t'.join([seqKey, trfContentDic[seqKey]['uid'], ','.join([str(s) for s in trfContentDic[seqKey]['count']]), ','.join(['%.3f'%(round(s, 3)) for s in trfContentDic[seqKey]['RPM']])]))
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
            readcountListTmp = [int(tmp) for tmp in content[2].split(',')]
            RP100KListTmp = [float(tmp) for tmp in content[3].split(',')]
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
        outf_tmp = open((Path(tRF_dir)/(sample+'.potential_tRFs.summary.report')), 'w')
        outf_tmp.write('amino acid\tCounts\tRP100K\tUnique reads\n')
        aaTypeList = []
        aaTypeDict = {}
        with open((Path(tRF_dir)/(sample+'.potential_tRFs.report')), 'w') as outf:
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
        readInforList = load_data_new((Path(tRF_dir)/(sample+'.potential_tRFs.report')))
        dcValue, rhomin, deltamin, RP100KCutoff = 3.0, 5.0, 8.0, 10.0
        tRNANameList = []
        tRNANameDic = {}
        # in Unified_tRNANameList, 'pre_Homo_sapiens_tRNA-Cys-GCA-2-4_trailer' will be treated as 'Homo_sapiens_tRNA-Cys-GCA-2-4'
        Unified_tRNANameList = []
        tRNANameList_selected = []
        Unified_tRNANameList_selected = []
        potential_tRFs_rep1 = sample+'.potential_tRFs.report'
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
        potential_tRFs_file = sample+'.potential_tRFs.clusters.detail'
        tRFs_report_file = sample+'.tRFs.report.tsv'
        outf_1 = open(Path(tRF_dir)/potential_tRFs_file, 'w')
        outf_2 = open(Path(tRF_dir)/tRFs_report_file, 'w')
        for indexId, readInforDic in enumerate(readInforList):
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
