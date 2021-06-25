import os
from pathlib import Path
import subprocess
from reportlab.lib import colors
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from math import sqrt, asin, cos, sin
from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna
import re
from mirge.classes.exportHTML import FormatJS

def Shifting(xcDic, ycDic):
    minx = min([xcDic[key] for key in xcDic.keys()])
    maxx = max([xcDic[key] for key in xcDic.keys()])
    miny = min([ycDic[key] for key in ycDic.keys()])
    maxy = max([ycDic[key] for key in ycDic.keys()])

    ## now place sec-structure in upper right corner
    shiftx = abs(minx)+10;
    shifty = abs(miny)+10;

    #shift everything to printable area;
    if minx < 0:
        for i in range(len(xcDic)):
            xcDic[i] = xcDic[i] + shiftx 
    if miny < 0:
        for i in range(len(ycDic)):
            ycDic[i] = ycDic[i] + shifty
    if maxx > 595.27:
        for i in range(len(xcDic)):
            xcDic[i] = xcDic[i] - (minx +10)
    if maxy > 841.89:
        for i in range(len(ycDic)):
            ycDic[i] = ycDic[i] - (miny + 10)

def centroidOfRectangle(minx, maxx, miny, maxy):
    return ((minx+maxx)/2, (miny+maxy)/2)

def lenAndWidOfRectangle(minx, maxx, miny, maxy):
    return (maxx-minx, maxy-miny)

def chunkListNew(number):
    outList = []
    if number <= 35:
        outList.append(range(number))
    else:
        outList.append(range(35))
        newList = range(35, number)
        newListTmp = [newList[i:i+67] for i in range(0, len(newList), 67)]
        for itmTmp in newListTmp:
            outList.append(itmTmp)
    return outList


def creatPDF(sampleName, novelmiRNANameNew, probability, chr, startPos, endPos, strand, armType, readCountSumMatureMiRNA, totalReadCountSum, matureMiRNAPrecusorSeq, matureMiRNAPrecusorStr, matureMiRNASeq, passengerMiRNASeq, psFile, outFile, matureReadContentlist, passengerReadContentlist):
    matureMiRNALocationList = range(matureMiRNAPrecusorSeq.index(matureMiRNASeq), matureMiRNAPrecusorSeq.index(matureMiRNASeq)+len(matureMiRNASeq))
    if passengerMiRNASeq != 'None':
        try:
            starMiRNALocationList = range(matureMiRNAPrecusorSeq.index(passengerMiRNASeq), matureMiRNAPrecusorSeq.index(passengerMiRNASeq)+len(passengerMiRNASeq))
        except ValueError:
            starMiRNALocationList = []
    else:
        starMiRNALocationList = []

    c = canvas.Canvas(outFile, pagesize=A4)
    xposshift = 15
    downy = 20
    # Write the titile
    c.setFont('Courier', 15)
    c.drawString(xposshift+180, downy+750, 'miRge 3.0 Novel miRNA')

    # Input the sequence information
    c.setFont('Courier', 7)
    c.drawString(xposshift+15, downy+700, 'Name: ')
    c.drawString(xposshift+148, downy+700, novelmiRNANameNew)
    c.drawString(xposshift+15, downy+690, 'Sample name: ')
    c.drawString(xposshift+148, downy+690, sampleName)
    c.drawString(xposshift+15, downy+680, 'Probability: ')
    c.drawString(xposshift+148, downy+680, '%.5f'%(float(probability)))
    c.drawString(xposshift+15, downy+670, 'Chr(GRCh38): ')
    c.drawString(xposshift+148, downy+670, chr)
    c.drawString(xposshift+15, downy+660, 'Start position of mature miRNA: ')
    c.drawString(xposshift+148, downy+660, str(startPos))
    c.drawString(xposshift+15, downy+650, 'End position of mature miRNA: ')
    c.drawString(xposshift+148, downy+650, str(endPos))
    c.drawString(xposshift+15, downy+640, 'Strand: ')
    c.drawString(xposshift+148, downy+640, strand)
    c.drawString(xposshift+15, downy+630, 'Total read count: ')
    c.drawString(xposshift+148, downy+630, str(totalReadCountSum))
    c.drawString(xposshift+15, downy+620, 'Mature miRNA read count: ')
    c.drawString(xposshift+148, downy+620, str(readCountSumMatureMiRNA))
    c.drawString(xposshift+15, downy+610, 'Passenger miRNA read count: ')
    c.drawString(xposshift+148, downy+610, str(int(totalReadCountSum)-int(readCountSumMatureMiRNA)))
    
    # Draw secondary structures.
    xcDic = {}
    ycDic = {}
    bpDic = {}
    twisted = 0 ## if mature sequence comes after loop ('arm3') twisted is 1
    if armType == 'arm3':
        twisted = 1

    mb = matureMiRNAPrecusorSeq.index(matureMiRNASeq)
    me = mb + len(matureMiRNASeq) - 1     # mature end coordinate
    centering_x = 0
    with open(psFile, 'r') as infTmp:
        line = infTmp.readline()
        flag = 0
        counter = 0
        while line != '':
            if line.strip().startswith('/sequence'):
                line = infTmp.readline()
                sequence = line.strip()[:-1]
            if line.strip() == '/coor [':
                flag = 1
                line = infTmp.readline()
            if line.strip() == '] def':
                flag = 0
            if line.strip() == '/pairs [':
                flag = 2
                line = infTmp.readline()
            if flag == 1:
                xc = float(line.strip()[1:-1].split(' ')[0])
                yc = float(line.strip()[1:-1].split(' ')[1])
                xcDic.update({counter:xc})
                ycDic.update({counter:yc})
                counter = counter + 1
            if flag == 2:
                if twisted == 1:
                    bpo2r = int(line.strip()[1:-1].split(' ')[0])
                    bpo1r = int(line.strip()[1:-1].split(' ')[1])
                else:
                    bpo2r = int(line.strip()[1:-1].split(' ')[1])
                    bpo1r = int(line.strip()[1:-1].split(' ')[0])
                bpDic.update({bpo1r-1:bpo2r-1})
            line = infTmp.readline()
    centering_x = matureMiRNALocationList[len(matureMiRNALocationList)//2]
    
    #print centering_x
    Shifting(xcDic, ycDic)
    
    #Scale the secondary structure.
    #scaling now structure
    # scaling center coordinates on page
    ry = 300
    scx = 0 ## scaling center coordinates on page
    scy = 0
    
    minx = min([xcDic[key] for key in xcDic.keys()])
    maxx = max([xcDic[key] for key in xcDic.keys()])
    miny = min([ycDic[key] for key in ycDic.keys()])
    maxy = max([ycDic[key] for key in ycDic.keys()])
    
    
    scfactor = ry/(maxy-miny)
    if scfactor < 1:
        for i in range(len(xcDic)):
            tx = xcDic[i] - scx
            ty = ycDic[i] - scy
            xcDic[i] = tx*scfactor + scx
            ycDic[i] = ty*scfactor + scy
    
    # mirror sequence so that when it rotate, the 5' will be on the top.
    for i in range(len(xcDic)):
        xcDic[i] = -xcDic[i]
    Shifting(xcDic, ycDic)
    
    # Determine the rotated angle phi
    # middle part of the mature RNA is the rotation center
    ax = xcDic[centering_x]
    ay = ycDic[centering_x]
     
    # point relativ to center
    if twisted == 1:
        bx = xcDic[centering_x-1]
        by = ycDic[centering_x-1]
    else:
        bx = xcDic[centering_x+1]
        by = ycDic[centering_x+1]
    gk = by-ay
    ak = bx-ax
    r = sqrt((ak**2)+(gk**2))
    phi = asin(gk/r)
    if (bx < ax and by > ay):
        phi = 3.141593-phi
    if (bx <= ax and by <= ay):
        phi = phi*(-1)
        phi == phi + 3.141593
    #This phi is abandoned, so only set phi to be 3.141593*1.5
    phi = 3.141593*1.5
    # Now rotate every point in a designated angle of phi
    for i in range(len(xcDic)):
        if i == centering_x:
            pass
        else:
            bx = xcDic[i]
            by = ycDic[i]
            gk = by-ay
            ak = bx-ax
            r = sqrt((ak**2)+(gk**2))
            alpha = asin(gk/r)
            if bx < ax and by > ay:
                alpha = 3.141593-alpha
            if bx <= ax and by <= ay:
                alpha = alpha*(-1)
                alpha = alpha + 3.141593
            alpha = alpha - phi
            xcDic[i] = ax + r*cos(alpha)
            ycDic[i] = ay + r*sin(alpha)
    Shifting(xcDic, ycDic)
    
    # Translocate the structure to the proper postion
    minx = min([xcDic[key] for key in xcDic.keys()])
    maxx = max([xcDic[key] for key in xcDic.keys()])
    miny = min([ycDic[key] for key in ycDic.keys()])
    maxy = max([ycDic[key] for key in ycDic.keys()])
    cnetroidX, centroidY = centroidOfRectangle(minx, maxx, miny, maxy)
    #print (cnetroidX, centroidY)
    scx = 95
    scy = 380
    tx = scx - cnetroidX
    ty = scy - centroidY
    for i in range(len(xcDic)):
        xcDic[i] = xcDic[i] + tx
        ycDic[i] = ycDic[i] + ty

    # Draw other parts of the pdf
    y = 300 #300;
    x = 300
    c.setLineWidth(0.5)
    c.setStrokeColor(colors.grey)

    # Draw lines
    for i in range(len(xcDic)):
        if i > 0:
            c.line(xcDic[i-1]+2.5+x,ycDic[i-1]+2+y,xcDic[i]+2.5+x,ycDic[i]+2+y)

    # Draw bases
    c.setFont('Courier', 10*scfactor)
    for i in range(len(sequence)):
        if i in matureMiRNALocationList:
            c.setFillColor(colors.red)
        elif i in starMiRNALocationList:
            c.setFillColor(colors.lightskyblue)
        else:
            c.setFillColor(colors.black)
        c.drawString(xcDic[i]+x,ycDic[i]+y,sequence[i])

    # Draw base pairing
    scfactorl= 0.4
    c.setLineWidth(0.6)
    c.setStrokeColor(colors.black)
    for key in bpDic.keys():
        dx = abs(xcDic[key] - xcDic[bpDic[key]])
        dy = abs(ycDic[key] - ycDic[bpDic[key]])

        dx1 = (dx-scfactorl*dx)/2
        dy1 = (dy-scfactorl*dy)/2

        if xcDic[key] > xcDic[bpDic[key]]:
            fx = xcDic[key] - dx1
            tox = xcDic[bpDic[key]] + dx1
        else:
            fx = xcDic[key] + dx1
            tox = xcDic[bpDic[key]] - dx1

        if ycDic[key] > ycDic[bpDic[key]]:
            fy = ycDic[key] - dy1
            toy = ycDic[bpDic[key]] + dy1
        else:
            fy = ycDic[key] + dy1
            toy = ycDic[bpDic[key]] - dy1
        #c.line(fx+2.5+x,fy+2+y,tox+2.5+x,toy+2+y)
        c.line(fx+2.5*scfactor+x,fy+2*scfactor+y,tox+2.5*scfactor+x,toy+2*scfactor+y)
    # Draw 5' and 3'
    c.setFillColor(colors.black)
    c.setFont('Courier', 8)
    c.drawString(xcDic[0]+2.5+x-12,ycDic[0]+2+y-5, "5'")
    c.drawString(xcDic[len(xcDic)-1]+2.5+x-12,ycDic[len(xcDic)-1]+2+y-5, "3'")
    
    # Draw histogram.
    # draw y axis
    multiplier = 3.6 # minimal distance between two letters
    lstruct_multi = (len(matureMiRNAPrecusorSeq)+2)*multiplier
    c.setLineWidth(2.5)
    c.setStrokeColor(colors.black)
    c.line(xposshift+15, downy+550+1, xposshift+15, downy+450-1)
    # draw x axis
    c.setStrokeColor(colors.grey)
    c.line(xposshift+15, downy+450, xposshift+15+lstruct_multi, downy+450)
    splitXScale = len(matureMiRNAPrecusorSeq)//10
    c.setFont('Courier', 5)
    c.setFillColor(colors.black)
    c.setStrokeColor(colors.grey)
    for i in range(1, splitXScale+1):
        c.line(xposshift+15+i*multiplier*10, downy+450, xposshift+15+i*multiplier*10, downy+453)
        c.drawString(xposshift+12+i*multiplier*10, downy+440, str(i)+'0')
    # draw x axis arrow
    c.setLineWidth(2.5)
    c.setStrokeColor(colors.black)
    c.line(xposshift+12+lstruct_multi, downy+453, xposshift+15+lstruct_multi, downy+450)
    c.line(xposshift+15+lstruct_multi, downy+450, xposshift+12+lstruct_multi, downy+447)
    # draw y axis arrow
    #c.line(xposshift+12, downy+557, xposshift+15, downy+560)
    #c.line(xposshift+15, downy+560, xposshift+18, downy+557)
    # add label of x axis and y axis
    c.setFont('Courier', 8)
    c.drawString(xposshift+15+lstruct_multi*0.5-10, downy+430, "Length(bp)")
    c.setFont('Courier', 6)
    c.drawString(xposshift+0, downy+560, "Percentage")
    # draw scale in y axis
    c.setFont('Courier', 6)
    # 1.0
    c.line(xposshift+15-1.0, downy+550, xposshift+18, downy+550)
    c.drawString(xposshift-5, downy+548, "100%")
    # 0.75
    c.line(xposshift+15, downy+525, xposshift+18, downy+525)
    c.drawString(xposshift-3, downy+523, "75%")
    # 0.5
    c.line(xposshift+15, downy+500, xposshift+18, downy+500)
    c.drawString(xposshift-3, downy+498, "50%")
    # 0.25
    c.line(xposshift+15, downy+475, xposshift+18, downy+475)
    c.drawString(xposshift-3, downy+473, "25%")
    # 0
    c.drawString(xposshift+7, downy+448, "0")
    # Draw the detailed reads information in the cluster.
    # Calculate the frequency of the nucleotide on genome for each postion
    totalReadCount = 0
    for item in matureReadContentlist:
        totalReadCount = totalReadCount + int(item[1])
    if passengerReadContentlist != 'None':
        for item in passengerReadContentlist:
            totalReadCount = totalReadCount + int(item[1])
    frequencyList = []
    # Calculate predefined pdf loci for alignment characters
    position_Dic = {}
    counter = 0
    for i in range(200):
        position_Dic.update({i:xposshift+multiplier+15+i*multiplier})
    
    for index, nucleotide in enumerate(matureMiRNAPrecusorSeq):
        countTmp = 0
        for matureRead in matureReadContentlist:
            if matureRead[0][index] == nucleotide:
                countTmp = countTmp + int( matureRead[1])
        if passengerReadContentlist != 'None':
            for passengerRead in passengerReadContentlist:
                if passengerRead[0][index] == nucleotide:
                    countTmp = countTmp + int( passengerRead[1])
        frequencyList.append(float(countTmp)/totalReadCount)
    # Draw histogram with different clolor
    #print matureMiRNALocationList
    for i in range(len(xcDic)):
        # Assign colors
        if i in matureMiRNALocationList:
            c.setStrokeColor(colors.red)
        elif i in starMiRNALocationList:
            c.setStrokeColor(colors.lightskyblue)
        else:
            c.setStrokeColor(colors.black)
        if i == 0:
            c.line(xposshift+15, downy+450, position_Dic[i], downy+450+frequencyList[i]*100)
        else:
            c.line(position_Dic[i-1], downy+450+frequencyList[i-1]*100, position_Dic[i], downy+450+frequencyList[i]*100)
    # Draw alignment between precursor sequence and read sequences
    for i in range(len(matureMiRNAPrecusorSeq)):
        c.setFont('Courier', 5)
        if i in matureMiRNALocationList:
            c.setFillColor(colors.red)
        elif i in starMiRNALocationList:
            c.setFillColor(colors.lightskyblue)
        else:
            c.setFillColor(colors.black)
        c.drawString(position_Dic[i],downy+400,matureMiRNAPrecusorSeq[i])
    c.setFont('Courier', 5)
    c.setFillColor(colors.black)
    c.drawString(xposshift+20+lstruct_multi, downy+400, "-3'")
    c.drawString(xposshift+5, downy+400, "5'-")
    for i, subitem in enumerate(matureMiRNAPrecusorStr):
        c.drawString(position_Dic[i], downy+390, subitem)
    c.setFont('Courier-Bold', 5)
    c.setFillColor(colors.black)
    c.drawString(xposshift+20+ lstruct_multi, downy+390,'reads')
    #c.drawString(xposshift+60+ lstruct_multi, downy+390,'sample')
    c.setFont('Courier', 5)
    c.setFillColor(colors.black)
    minusValue = 10
    lowLimit = downy+75
    if passengerReadContentlist != 'None':
        totalSeqCount = len(matureReadContentlist) + len(passengerReadContentlist)
        totalReadContentlist = matureReadContentlist + passengerReadContentlist
        chunkLabelList = chunkListNew(totalSeqCount)
    else:
        totalSeqCount = len(matureReadContentlist)
        totalReadContentlist = matureReadContentlist
        chunkLabelList = chunkListNew(totalSeqCount)
    for i, chunk in enumerate(chunkLabelList):
        minusValue = 10
        if i == 0:
            for index in chunk:
                for j, subitem in enumerate(totalReadContentlist[index][0]):
                    c.drawString(position_Dic[j], downy+390-minusValue, subitem)
                c.drawString(xposshift+20+ lstruct_multi, downy+390-minusValue,str(totalReadContentlist[index][1]))
                #c.drawString(xposshift+60+ lstruct_multi, downy+390-minusValue,sampleName)
                minusValue = minusValue + 10
        else:
            c.showPage()
            c.setFont('Courier', 5)
            c.setFillColor(colors.black)
            for index in chunk:
                for j, subitem in enumerate(totalReadContentlist[index][0]):
                    c.drawString(position_Dic[j], downy+710-minusValue, subitem)
                c.drawString(xposshift+20+ lstruct_multi, downy+710-minusValue,str(totalReadContentlist[index][1]))
                #c.drawString(xposshift+60+ lstruct_multi, downy+710-minusValue,sampleName)
                minusValue = minusValue + 10
    c.showPage()
    c.save()

def getMiRNAPosition(miRNASeq, RNASequence):
    startPos = RNASequence.find(miRNASeq)
    endPos = startPos + len(miRNASeq)
    return [startPos, endPos]

def getMiRNAStructure(miRNASeq, RNASequence, RNAStructure):
    [startPos, endPos] = getMiRNAPosition(miRNASeq, RNASequence)
    mirnaStruct = RNAStructure[startPos:endPos]
    return mirnaStruct

def removeDash(seq):
    newSeq = ''
    for unit in seq:
        if unit != '-':
            newSeq = newSeq + unit
    return newSeq


def headDashCount2(clusterSeq):
    # clusterSeq is: CTAATACTGCCTGGTAATGATGACGG-
    headDashCount = 0
    for symbol in clusterSeq:
        if symbol == '-':
            headDashCount = headDashCount + 1
        else:
            break
    return headDashCount

def tailDashCount2(clusterSeq):
    # clusterSeq is: CTAATACTGCCTGGTAATGATGACGG-
    tailDashCount = 0
    for symbol in clusterSeq[::-1]:
        if symbol == '-':
            tailDashCount = tailDashCount + 1
        else:
            break
    return tailDashCount

def chunkInto2(arr):
    # length of list must be even and larger than 2
    return [arr[i:i+2] for i in range(0, len(arr), 2)]

def padClusteredList(matureReadContentlistRaw, matureMiRNASeq, matureMiRNAPrecusorSeq, startPos, endPos, strand):
    startTmp = matureMiRNAPrecusorSeq.index(matureMiRNASeq)
    endTmp = matureMiRNAPrecusorSeq.index(matureMiRNASeq) + len(matureMiRNASeq)
    headLen = startTmp
    tailLen = len(matureMiRNAPrecusorSeq) - endTmp
    if strand == '+':
        precursorHeadLocation = int(startPos) - headLen
        precursorTailLocation = int(endPos) + tailLen
    else:
        precursorHeadLocation = int(startPos) - tailLen
        precursorTailLocation = int(endPos) + headLen
    newReadContentList = []
    for readContent in matureReadContentlistRaw:
        readSeq = readContent[0]
        readSeqStart = int(readContent[1])
        readSeqEnd = int(readContent[2])
        readCount = readContent[3]
        if strand == '+':
            readSeqPadded = '-'*(readSeqStart-precursorHeadLocation)+readSeq+'-'*(precursorTailLocation-readSeqEnd)
        else:
            readSeqPadded = '-'*(precursorTailLocation-readSeqEnd)+readSeq+'-'*(readSeqStart-precursorHeadLocation)
        newReadContentList.append((readSeqPadded, readCount))
    return newReadContentList


def write_novel_report(novelmiRNALListFile, featureFile, clusterFile, rnafoldCmdTmp, outputdir2, files):
    dir_tmp = Path((Path(outputdir2).resolve().parents[0]))/(files+'_novel_miRNAs')
    os.mkdir(dir_tmp)
    #samppleName_tmp = '_'.join(novelmiRNALListFile.split('_')[:-3])
    novelmiRNALOriginalList = []
    clusterNameProbabilityDic = {}
    with open(novelmiRNALListFile, 'r') as inf:
        line = inf.readline()
        probabilityIndex = line.strip().split(',').index('probability')
        clusterNameIndex = line.strip().split(',').index('clusterName')
        line = inf.readline()
        while line != '':
            content = line.strip().split(',')
            probability = content[probabilityIndex]
            clusterName = content[clusterNameIndex]
            if clusterName not in novelmiRNALOriginalList:
                novelmiRNALOriginalList.append(clusterName)
                clusterNameProbabilityDic.update({clusterName:probability})
            line = inf.readline()
    # The majority of the paired miRNAs' pruned precursor RNA sequences should be identical. However, the pruned precursor RNA sequences might be slightly different
    # at tail and head part. This situation happens at a quite low probility due to pruning error.
    # Therefore, the paired miRNA's precursor miRNA of the identified novel miRNA should be corrected in order to plot.
    with open(featureFile, 'r') as inf:
        totalContent = inf.readlines()

    with open(featureFile, 'r') as inf:
        line = inf.readline()
        content = line.strip().split('\t')
        clusterNameIndex = content.index('clusterName')
        pruned_precusor_seqIndex = content.index('pruned_precusor_seq')
        pruned_precusor_strIndex = content.index('pruned_precusor_str')
        upstreamDistanceIndex = content.index('upstreamDistance')
        downstreamDistanceIndex = content.index('downstreamDistance')
        armTypeIndex = content.index('armType')
        pair_stateIndex = content.index('pair_state')
        i = 1
        line = inf.readline()
        while line != '':
            content = line.strip().split('\t')
            if content[clusterNameIndex] in novelmiRNALOriginalList:
                armTypeTmp = content[armTypeIndex]
                pair_state = content[pair_stateIndex]
                pruned_precusor_seq = content[pruned_precusor_seqIndex]
                pruned_precusor_str = content[pruned_precusor_strIndex]
                upstreamDistance = content[upstreamDistanceIndex]
                downstreamDistance = content[downstreamDistanceIndex]
                if pair_state == 'Yes':
                    if upstreamDistance != 'None':
                        if int(upstreamDistance) <= 44:
                            correctLineContent = totalContent[i-1].strip().split('\t')
                            if (correctLineContent[pair_stateIndex]) == 'Yes' and (correctLineContent[armTypeIndex] in ['arm5', 'arm3']) and (correctLineContent[armTypeIndex] != armTypeTmp):
                                correctLineContent[pruned_precusor_seqIndex] = pruned_precusor_seq
                                correctLineContent[pruned_precusor_strIndex] = pruned_precusor_str
                            totalContent[i-1] = '\t'.join(correctLineContent)+'\n'
                        elif int(downstreamDistance) <= 44:
                            correctLineContent = totalContent[i+1].strip().split('\t')
                            if (correctLineContent[pair_stateIndex]) == 'Yes' and (correctLineContent[armTypeIndex] in ['arm5', 'arm3']) and (correctLineContent[armTypeIndex] != armTypeTmp):
                                correctLineContent[pruned_precusor_seqIndex] = pruned_precusor_seq
                                correctLineContent[pruned_precusor_strIndex] = pruned_precusor_str
                            totalContent[i+1] = '\t'.join(correctLineContent)+'\n'
                        else:
                            pass
            line = inf.readline()
            i = i + 1

    corrFile = str(Path(outputdir2)/(files + '_corrected.tsv'))
    with open(corrFile, 'w') as outf:
        for t, totalContentTmp in enumerate(totalContent):
            # Correct the armType if it is wrongly alocated because the armType is based on the stableClusterSeq. This sequnce may be too short.
            # So use the clusterSeq to recorrect the armType.
            if t == 0:
                content = totalContentTmp.strip().split('\t')
                clusterNameIndex = content.index('clusterName')
                clusterSeqIndex = content.index('clusterSeq')
                stableClusterSeqIndex = content.index('stableClusterSeq')
                pruned_precusor_seqIndex = content.index('pruned_precusor_seq')
                pruned_precusor_strIndex = content.index('pruned_precusor_str')
                armTypeIndex = content.index('armType')
                outf.write(totalContentTmp)
            else:
                content = totalContentTmp.strip().split('\t')
                clusterSeq = content[clusterSeqIndex]
                clusterSeqNew = ''
                for nucl in clusterSeq:
                    if nucl == 'T':
                        clusterSeqNew = clusterSeqNew + 'U'
                    else:
                        clusterSeqNew = clusterSeqNew + nucl
                stableClusterSeq = content[stableClusterSeqIndex]
                pruned_precusor_seq = content[pruned_precusor_seqIndex]
                pruned_precusor_str = content[pruned_precusor_strIndex]
                clusterSeqStr = getMiRNAStructure(clusterSeqNew, pruned_precusor_seq, pruned_precusor_str)
                pattern = re.compile('\(.*\)')
                if pattern.search(clusterSeqStr) is not None:
                    if content[armTypeIndex] == 'loop':
                        outf.write(totalContentTmp)
                    else:
                        content[armTypeIndex] = 'loop'
                        outf.write('\t'.join(content)+'\n')
                else:
                    outf.write(totalContentTmp)

    clusterNameFeatureDic = {}
    precursorSeqclusterNameDic = {}
    with open(corrFile, 'r') as inf:
        line = inf.readline()
        content = line.strip().split('\t')
        clusterNameIndex = content.index('clusterName')
        seqCountIndex = content.index('seqCount')
        readCountSumIndex = content.index('readCountSum')
        stableClusterSeqIndex = content.index('stableClusterSeq')
        alignedClusterSeqLabel = content.index('alignedClusterSeq')
        headUnstableLengthLabel = content.index('headUnstableLength')
        tailUnstableLengthLabel = content.index('tailUnstableLength')
        precusorSeqIndex = content.index('pruned_precusor_seq')
        precusorStrIndex = content.index('pruned_precusor_str')
        armTypeIndex = content.index('armType')
        line = inf.readline()
        while line != '':
            content = line.strip().split('\t')
            clusterName = content[clusterNameIndex]
            seqCount = content[seqCountIndex]
            readCountSum = content[readCountSumIndex]
            stableClusterSeq = str(Seq(content[stableClusterSeqIndex]).transcribe())
            alignedClusterSeq = content[alignedClusterSeqLabel]
            headUnstableLength = int(content[headUnstableLengthLabel])
            tailUnstableLength = int(content[tailUnstableLengthLabel])
            precusorSeq = content[precusorSeqIndex]
            precusorStr = content[precusorStrIndex]
            armType = content[armTypeIndex]
            strand = clusterName[-1]
            chr = clusterName.split(':')[2]
            startPos = int(clusterName.split(':')[3][:-1].split('_')[0].strip())
            endPos = int(clusterName.split(':')[3][:-1].split('_')[1].strip())
            headDashCountTmp = headDashCount2(alignedClusterSeq)
            tailDashCountTmp = tailDashCount2(alignedClusterSeq)
            if strand == '+':
                startPos = startPos - headDashCountTmp + headUnstableLength
                endPos = endPos + tailDashCountTmp - tailUnstableLength
            else:
                startPos = startPos - tailDashCountTmp + tailUnstableLength
                endPos = endPos + headDashCountTmp - headUnstableLength
            if precusorSeq != 'None':
                if clusterName not in clusterNameFeatureDic.keys():
                    clusterNameFeatureDic.update({clusterName:[chr, startPos, endPos, strand, stableClusterSeq, seqCount, readCountSum, armType, precusorSeq, precusorStr]})
                if (precusorSeq, chr, strand) not in precursorSeqclusterNameDic.keys():
                    precursorSeqclusterNameDic.update({(precusorSeq, chr, strand):[clusterName]})
                else:
                    precursorSeqclusterNameDic[(precusorSeq, chr, strand)].append(clusterName)
            line = inf.readline()

    clusterNameClusterSeqDic = {}
    # Parse the clusterFile to get the detailed reads information for each cluster.
    with open(clusterFile, 'r') as inf:
        contentTmp = inf.readlines()
    labelList = []
    for index, item in enumerate(contentTmp):
        if 'Cluster Name:' in item:
            labelList.append(index)

    for k in range(len(labelList)):
        if k != len(labelList)- 1:
            subContentTmp = contentTmp[labelList[k]:labelList[k+1]]
        else:
            subContentTmp = contentTmp[labelList[k]:]
        clusterNameTmp = subContentTmp[0].strip().split(' ')[2]
        if clusterNameTmp in clusterNameFeatureDic.keys():
            #print clusterNameTmp
            startTmp = int(clusterNameTmp[:-1].split(':')[-1].split('_')[0])
            endTmp = int(clusterNameTmp[:-1].split(':')[-1].split('_')[1])
            if clusterNameTmp[-1] == '+':
                dashedClusterSeqStart = startTmp - headDashCount2(subContentTmp[4].strip())
                dashedClusterSeqEnd = endTmp + tailDashCount2(subContentTmp[4].strip())
            else:
                dashedClusterSeqStart = startTmp - tailDashCount2(subContentTmp[4].strip())
                dashedClusterSeqEnd = endTmp + headDashCount2(subContentTmp[4].strip())
            readContentlist = []
            for subitem in subContentTmp[5:]:
                if clusterNameTmp[-1] == '+':
                    ReadSeqStart = dashedClusterSeqStart + headDashCount2(subitem.split('\t')[0])
                    ReadSeqEnd = dashedClusterSeqEnd - tailDashCount2(subitem.split('\t')[0])
                else:
                    #ReadSeqStart = dashedClusterSeqStart + headDashCount2(subitem.split('\t')[0])-2
                    #ReadSeqEnd = dashedClusterSeqEnd - tailDashCount2(subitem.split('\t')[0])-2
                    
                    ReadSeqStart = dashedClusterSeqStart + tailDashCount2(subitem.split('\t')[0])
                    ReadSeqEnd = dashedClusterSeqEnd - headDashCount2(subitem.split('\t')[0])
                readCount = int(subitem.strip().split('\t')[1])
                #readSeq = removeDash(subitem.split('\t')[0])
                readSeq = str(Seq(removeDash(subitem.split('\t')[0])).transcribe())
                readContentlist.append((readSeq, ReadSeqStart, ReadSeqEnd, readCount))
            clusterNameClusterSeqDic.update({clusterNameTmp:readContentlist})
    
    #print(clusterNameClusterSeqDic)
    #print()
    
    # Output the novel miRNA report csv file
    #print precursorSeqclusterNameDic[('CUGACUGCCGAGGGGGCCCUGGCCUGGAUCCAUGCUGGGCAGAAGCAGCUGGACACUGACCAGGACCCCCCAGGGCCGGAGGAACC', 'chr9', '+')]
    finalFile = str(Path(dir_tmp)/(files+'_novel_miRNAs_report.csv'))
    #finalFile = str(Path(outputdir2)/(files+'_novel_miRNAs_report.csv'))
    html_data = FormatJS(str(Path(outputdir2).resolve().parents[0]))
    #print(str(Path(outputdir2).resolve().parents[0])) files+'_novel_miRNAs'
    #"14","<a href=\"SRR8557389_novel_miRNAs/SRR8557389_novel_miRNA_14.pdf\" target=\"_blank\">SRR8557389_novel_miRNA_14</a>","0.829119033","chr1","26554575","26554595","CUCCUGCCCUCCUUGCUGUAG","29"
    outf1 = open(finalFile, 'w')
    outf1.write('Novel miRNA name,Probability,Chr,Start Pos,End Pos,Strand,Mature miRNA sequence,Arm type,Passenger miRNA sequence,Mature miRNA read Count,Passenger miRNA read Count,Precursor miRNA sequence,Precursor miRNA structure\n')
    i = 1
    while len(novelmiRNALOriginalList) >= 1:
        js_nmirVar= ""
        novelmiRNA = novelmiRNALOriginalList[0]
        sampleName = str(files)
        #sampleName = '_'.join(novelmiRNA.split(':')[0].split('_')[2:])
        precursorSeq = clusterNameFeatureDic[novelmiRNA][-2]
        clusterNameList = precursorSeqclusterNameDic[(precursorSeq, clusterNameFeatureDic[novelmiRNA][0], clusterNameFeatureDic[novelmiRNA][3])]
        for clusterNameListTmp in chunkInto2(clusterNameList):
            flag_padClust = 0
            if len(clusterNameListTmp) == 1:
                matureMiRNAName = clusterNameListTmp[0]
                passengerMiRNAName = 'None'
            elif len(clusterNameListTmp) == 2:
                if clusterNameListTmp[0] in clusterNameProbabilityDic.keys() and clusterNameListTmp[1] in clusterNameProbabilityDic.keys():
                    if int(clusterNameFeatureDic[clusterNameListTmp[0]][6]) >= int(clusterNameFeatureDic[clusterNameListTmp[1]][6]):
                        matureMiRNAName = clusterNameListTmp[0]
                        passengerMiRNAName = clusterNameListTmp[1]
                    else:
                        matureMiRNAName = clusterNameListTmp[1]
                        passengerMiRNAName = clusterNameListTmp[0]
                elif clusterNameListTmp[0] in clusterNameProbabilityDic.keys() and clusterNameListTmp[1] not in clusterNameProbabilityDic.keys():
                    matureMiRNAName = clusterNameListTmp[0]
                    passengerMiRNAName = clusterNameListTmp[1]
                else:
                    matureMiRNAName = clusterNameListTmp[1]
                    passengerMiRNAName = clusterNameListTmp[0]
                if clusterNameFeatureDic[passengerMiRNAName][7] == 'loop':
                    passengerMiRNAName = 'None'
            chr = clusterNameFeatureDic[matureMiRNAName][0]
            startPos = clusterNameFeatureDic[matureMiRNAName][1]
            endPos = clusterNameFeatureDic[matureMiRNAName][2]
            strand = clusterNameFeatureDic[matureMiRNAName][3]
            matureMiRNASeq = clusterNameFeatureDic[matureMiRNAName][4]
            readCountSumMatureMiRNA = clusterNameFeatureDic[matureMiRNAName][6]
            armType = clusterNameFeatureDic[matureMiRNAName][7]
            matureMiRNAPrecusorSeq = clusterNameFeatureDic[matureMiRNAName][8]
            matureMiRNAPrecusorStr = clusterNameFeatureDic[matureMiRNAName][9]
            matureReadContentlistRaw = clusterNameClusterSeqDic[matureMiRNAName]
            if passengerMiRNAName == 'None':
                passengerMiRNASeq = 'None'
                passengerReadContentlistRaw = 'None'
                passengerReadContentlist = 'None'
                readCountSumPassengerMiRNA = 'None'
                passengerStrand = 'None'
            else:
                passengerMiRNASeq = clusterNameFeatureDic[passengerMiRNAName][4]
                readCountSumPassengerMiRNA = clusterNameFeatureDic[passengerMiRNAName][6]
                passengerReadContentlistRaw = clusterNameClusterSeqDic[passengerMiRNAName]
                passengerStartPos = clusterNameFeatureDic[passengerMiRNAName][1]
                passengerEndPos = clusterNameFeatureDic[passengerMiRNAName][2]
                passengerStrand = clusterNameFeatureDic[passengerMiRNAName][3]
            # Pad the head and tail part with '.' according to the allignment to the precursor miRNA.
            try:
                startTmp = matureMiRNAPrecusorSeq.index(matureMiRNASeq)
                flag_padClust = 0
            except ValueError:
                flag_padClust = 1
            
            if flag_padClust == 1:
                continue
            else:
                matureReadContentlist = padClusteredList(matureReadContentlistRaw, matureMiRNASeq, matureMiRNAPrecusorSeq, startPos, endPos, strand)

            if passengerReadContentlistRaw != 'None':
                try:
                    #startTmp = matureMiRNAPrecusorSeq.index(passengerMiRNASeq)
                    passengerReadContentlist = padClusteredList(passengerReadContentlistRaw, passengerMiRNASeq, matureMiRNAPrecusorSeq, passengerStartPos, passengerEndPos, passengerStrand)
                    totalReadCountSum = int(readCountSumMatureMiRNA) + int(readCountSumPassengerMiRNA)
                except ValueError:
                    passengerReadContentlist = 'None'
                    totalReadCountSum = int(readCountSumMatureMiRNA)
            else:
                passengerReadContentlist = 'None'
                totalReadCountSum = int(readCountSumMatureMiRNA)
            if armType != 'loop':
                novelmiRNANameNew = str(files)+'_novel_miRNA_'+str(i)
                try:
                    outf1.write(','.join([novelmiRNANameNew, clusterNameProbabilityDic[matureMiRNAName], chr, str(startPos), str(endPos), strand, matureMiRNASeq, armType, passengerMiRNASeq, str(readCountSumMatureMiRNA), str(readCountSumPassengerMiRNA), matureMiRNAPrecusorSeq, matureMiRNAPrecusorStr]))
                    js_nmirVar = '"<a href=\\"'+str(files)+'_novel_miRNAs/'+str(files)+'_novel_miRNA_'+str(i)+'.pdf\\" target=\\"_blank\\">'+novelmiRNANameNew+'</a>","'+str(clusterNameProbabilityDic[matureMiRNAName])+'","'+chr+'","'+str(startPos)+'","'+str(endPos)+'","'+matureMiRNASeq+'","'+str(readCountSumMatureMiRNA)+'"'
                    html_data.addSerialNum(js_nmirVar)
                except KeyError:
                    outf1.write(','.join([novelmiRNANameNew, "1", chr, str(startPos), str(endPos), strand, matureMiRNASeq, armType, passengerMiRNASeq, str(readCountSumMatureMiRNA), str(readCountSumPassengerMiRNA), matureMiRNAPrecusorSeq, matureMiRNAPrecusorStr]))
                    js_nmirVar = '"<a href=\\"'+str(files)+'_novel_miRNAs/'+str(files)+'_novel_miRNA_'+str(i)+'.pdf\\" target=\\"_blank\\">'+novelmiRNANameNew+'</a>","1","'+chr+'","'+str(startPos)+'","'+str(endPos)+'","'+matureMiRNASeq+'","'+str(readCountSumMatureMiRNA)+'"'
                    html_data.addSerialNum(js_nmirVar)
                outf1.write('\n')
                # Prepare to plot the precurosr sturcuture, cluster seuqences into a pdf file.
                with open(str(Path(outputdir2)/(files+'_precusorTmp.fa')), 'w') as outf:
                    fa_tmp = '\n'.join(['>'+novelmiRNANameNew, matureMiRNAPrecusorSeq, matureMiRNAPrecusorStr])
                    outf.write(fa_tmp+'\n')
                f1 = str(files+'_precusorTmp.fa')
                f2 = str(files+'_precusorTmp.str')
                os.system('cd %s && %s -d 0 < %s > %s'%(Path(outputdir2), rnafoldCmdTmp, f1, f2))
                f3 = str(Path(outputdir2)/(files+'_novel_miRNA_'+str(i)+'_ss.ps'))
                f4 = str(Path(dir_tmp)/(files+'_novel_miRNA_'+str(i)+'.pdf'))
                try:
                    creatPDF(sampleName, novelmiRNANameNew, clusterNameProbabilityDic[matureMiRNAName], chr, startPos, endPos, strand, armType, readCountSumMatureMiRNA, totalReadCountSum, matureMiRNAPrecusorSeq, matureMiRNAPrecusorStr, matureMiRNASeq, passengerMiRNASeq, f3, f4, matureReadContentlist, passengerReadContentlist)
                except KeyError:
                    creatPDF(sampleName, novelmiRNANameNew, "1", chr, startPos, endPos, strand, armType, readCountSumMatureMiRNA, totalReadCountSum, matureMiRNAPrecusorSeq, matureMiRNAPrecusorStr, matureMiRNASeq, passengerMiRNASeq, f3, f4, matureReadContentlist, passengerReadContentlist)
                
                i = i + 1
            # Delete the clusterNames in clusterNameListTmp
            for clusterName in clusterNameListTmp:
                #print '%s is remvoed'%(clusterName)
                if clusterName in novelmiRNALOriginalList:
                    novelmiRNALOriginalList.remove(clusterName)     
    outf1.close()
