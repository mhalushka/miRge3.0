from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
import sys

def calculate_identity(sequenceA, sequenceB):
	#sequenceA is the clustered seq and sequenceB is the read
	#Returns the percentage of identical characters between two sequences.
	sa, sb, sl = sequenceA, sequenceB, len(sequenceA)
	matches = 0
	for i in range(sl):
		if sa[i] == sb[i] and sa[i] != '-' and  sb[i] != '-':
			matches = matches + 1
	return matches

def chunks(seq, n):
	return (seq[i:i+n] for i in range(0, len(seq), n))

class ReadCluster(object):
	threshold = 0.8
	headShiftStep = 3
	tailShiftStep = 6
	def __init__(self, chrSeqDic, chr, strand, startPos, endPos, clusterName, clusterSeq, readNameList, readSeqList, readCountList):
		self.chrSeqDic = chrSeqDic
		self.chr = chr
		self.strand = strand
		self.startPos = startPos
		self.endPos = endPos
		self.clusterName = clusterName
		self.clusterSeq = clusterSeq
		self.readNameList = readNameList
		self.readSeqList = readSeqList
		self.readCountList = readCountList

	def totalreadCount(self):
		return sum(self.readCountList)

	def totalseqCount(self):
		return len(self.readSeqList)

	def align2Standard(self):
		#align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
		alignSeqList = []
		clusterSeqTmp = self.clusterSeq
		exactMatchCount = 0
		#for seq in self.readSeqList:
		for i in range(len(self.readSeqList)):
			# Perform pairwise local alignment. And increase the penalty of gap opening and extending into -10 and -10,
			# so that there will be no gap in the alignned result.
			seq = self.readSeqList[i]
			alignTmp = pairwise2.align.localms(clusterSeqTmp, seq, 2, -1, -20, -20)
			clusterSeqTmpNew = alignTmp[0][0]
			seqTmp = alignTmp[0][1]
			identity = calculate_identity(clusterSeqTmpNew, seqTmp)
			if identity == len(seq):
				exactMatchCount = exactMatchCount + self.readCountList[i]
			else:
				pass
			if len(alignSeqList) == 0:
				alignSeqList.append(clusterSeqTmpNew)
				alignSeqList.append(seqTmp)
			else:
				if clusterSeqTmpNew == alignSeqList[0]:
					alignSeqList.append(seqTmp)
				else:
					headAdd1 = clusterSeqTmpNew.index(clusterSeqTmp)
					tailAdd1 = len(clusterSeqTmpNew)-headAdd1-len(clusterSeqTmp)
					headAdd2 = alignSeqList[0].index(clusterSeqTmp)
					tailAdd2 = len(alignSeqList[0])-headAdd2-len(clusterSeqTmp)
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
		return (alignSeqList, exactMatchCount)
	
	def locateStartPosition(self):
		alignSeqList1, exactMatchCount1 = self.align2Standard()
		clusterSeq = alignSeqList1[0]
		transformedListTmp = []
		for i in range(len(clusterSeq)):
			for j in range(1, len(alignSeqList1)):
				transformedListTmp.append(alignSeqList1[j][i])
		transformedList = list(chunks(transformedListTmp, len(alignSeqList1)-1))
		readCountSum = sum(self.readCountList)
		clusterSeqRatioList = []
		majorSeqRationList = []
		majorSeq = ""
		for i in range(len(transformedList)):
			item = transformedList[i]
			count_A = 0
			count_T = 0
			count_C = 0
			count_G = 0
			for j in range(len(item)):
				if item[j] == 'A':
					count_A = count_A + self.readCountList[j]
				elif item[j] == 'T':
					count_T = count_T + self.readCountList[j]
				elif item[j] == 'C':
					count_C = count_C + self.readCountList[j]
				elif item[j] == 'G':
					count_G = count_G + self.readCountList[j]
				else:
					pass
			if clusterSeq[i] == 'A':
				clusterSeqPositionRatio = float(count_A)/readCountSum
			elif clusterSeq[i] == 'T':
				clusterSeqPositionRatio = float(count_T)/readCountSum
			elif clusterSeq[i] == 'C':
				clusterSeqPositionRatio = float(count_C)/readCountSum
			elif clusterSeq[i] == 'G':
				clusterSeqPositionRatio = float(count_G)/readCountSum
			else:
				clusterSeqPositionRatio = 0
			clusterSeqRatioList.append(clusterSeqPositionRatio)
			listTmp = [[count_A, 'A'], [count_T, 'T'], [count_C, 'C'], [count_G, 'G']]
			listTmp.sort(reverse=True)
			majorSeqRationList.append(float(listTmp[0][0])/readCountSum)
			majorSeq = majorSeq + listTmp[0][1]
		headStartPosition = None
		tailStartPosition = None
		for i in range(len(majorSeqRationList)):
			if float(majorSeqRationList[i]) >= ReadCluster.threshold:
				headStartPosition = i
				break
		for i in range(-1, -len(majorSeqRationList)-1, -1):
			if float(majorSeqRationList[i]) >= ReadCluster.threshold:
				tailStartPosition = i
				break
		return (alignSeqList1, exactMatchCount1, clusterSeq, clusterSeqRatioList, majorSeq, majorSeqRationList, headStartPosition, tailStartPosition)

	def calculateFeature(self):
		alignSeqList2, exactMatchCount2, clusterSeq, clusterSeqRatioList, majorSeq, majorSeqRationList, headStartPosition, tailStartPosition = self.locateStartPosition()
		#alignSeqList2, exactMatchCount2 = self.align2Standard()
		#for subseq in alignSeqList2:
		#	print subseq
		#print ReadCluster.headShiftStep
		#print ReadCluster.tailShiftStep
		neclueotidePositionList = []
		for i in range(ReadCluster.headShiftStep+ReadCluster.tailShiftStep):
			if i <= ReadCluster.headShiftStep-1:
				postionLable = 'head_minus%d_'%(ReadCluster.headShiftStep-i)
			else:
				postionLable = 'tail_plus%d_'%(i+1-ReadCluster.headShiftStep)
			for neclueotide in ['templateNucleotide', 'TemplateNucleotide_percentage', 'nonTemplateNucleotide_percentage', 'A_percentage', 'T_percentage', 'C_percentage', 'G_percentage']:
				neclueotidePositionList.append(postionLable+neclueotide)

		if (headStartPosition is not None )and (tailStartPosition is not None):
			readCountSum = sum(self.readCountList)
			seqCount = len(self.readSeqList)
			exactMatchRatio = float(exactMatchCount2)/readCountSum
			headAdd = headStartPosition
			tailAdd = -tailStartPosition-1
			headUnstableLength = headAdd
			tailUnstableLength = tailAdd
			#print tailAdd
			if headAdd < ReadCluster.headShiftStep:
				for i in range(len(alignSeqList2)):
					alignSeqList2[i] = (ReadCluster.headShiftStep-headAdd)*'-'+alignSeqList2[i]
				headStartPosition = ReadCluster.headShiftStep
			else:
				pass
			if tailAdd < ReadCluster.tailShiftStep:
				for i in range(len(alignSeqList2)):
					alignSeqList2[i] = alignSeqList2[i]+(ReadCluster.tailShiftStep-tailAdd)*'-'
				tailStartPosition = 0 - ReadCluster.tailShiftStep
			else:
				pass
			adjustedClusterSeq = alignSeqList2[0]
			headDashCount = 0
			for i in alignSeqList2[0]:
				if i == '-':
					headDashCount += 1
				else:
					break
			tailDashCount = 0
			for i in reversed(alignSeqList2[0]):
				if i == '-':
					tailDashCount += 1
				else:
					break

			if self.strand == '+':
				genomeSeqStart = self.startPos - headDashCount
				genomeSeqEnd = self.endPos + tailDashCount
				tmp_new = self.chrSeqDic[self.chr][genomeSeqStart-1:genomeSeqEnd]
				templateSeq = tmp_new.upper()
			else:
				genomeSeqStart = self.startPos - tailDashCount
				genomeSeqEnd = self.endPos + headDashCount
				tmp_new = self.chrSeqDic[self.chr][genomeSeqStart-1:genomeSeqEnd]
				templateSeq = str(Seq(tmp_new.upper()).reverse_complement())
			alignSeqList2.insert(0, templateSeq)
			clusterSecondSeq = alignSeqList2[2]

			colList = list(range(headStartPosition-ReadCluster.headShiftStep, headStartPosition)) + list(range(tailStartPosition, tailStartPosition+ReadCluster.tailShiftStep))
			transformedListTmp = []
			for i in range(len(colList)):
				for j in range(2, len(alignSeqList2)):
					transformedListTmp.append(alignSeqList2[j][colList[i]])
			transformedList = list(chunks(transformedListTmp, len(alignSeqList2)-2))
			neclueotideCountList = []
			#print len(self.readCountList)
			#print len(transformedList)
			for i in range(len(transformedList)):
				try:
					templateNucleotide = templateSeq[colList[i]]
				except IndexError:
					print("Error happends at: %s"%(self.clusterName))
					print("Please check it.")
					print("templateSeq is: %s"%(templateSeq))
					print("index is:")
					print(colList[i])
					sys.exit(1)
				item = transformedList[i]
				#print len(item)
				count_A = 0
				count_T = 0
				count_C = 0
				count_G = 0
				count_Dash = 0
				for j in range(len(item)):
					if item[j] == 'A':
						count_A = count_A + self.readCountList[j]
					elif item[j] == 'T':
						count_T = count_T + self.readCountList[j]
					elif item[j] == 'C':
						count_C = count_C + self.readCountList[j]
					elif item[j] == 'G':
						count_G = count_G + self.readCountList[j]
					else:
						count_Dash = count_Dash + self.readCountList[j]

				if templateNucleotide == 'A':
					count_templateNucleotide = count_A
					count_nonTemplateNucleotide = count_T+count_C+count_G
				elif templateNucleotide == 'T':
					count_templateNucleotide = count_T
					count_nonTemplateNucleotide = count_A+count_C+count_G
				elif templateNucleotide == 'C':
					count_templateNucleotide = count_C
					count_nonTemplateNucleotide = count_A+count_T+count_G
				elif templateNucleotide == 'G':
					count_templateNucleotide = count_G
					count_nonTemplateNucleotide = count_A+count_T+count_C
				else:
					count_templateNucleotide = 0
					count_nonTemplateNucleotide = 0

				if count_A+count_T+count_C+count_G != 0:
					count_TemplateNucleotide_percentage = float(count_templateNucleotide)/(count_A+count_T+count_C+count_G)
					count_nonTemplateNucleotide_percentage = float(count_nonTemplateNucleotide)/(count_A+count_T+count_C+count_G)
					if templateNucleotide == 'A':
						if count_nonTemplateNucleotide != 0:
							count_A_percentage = 0
							count_T_percentage = float(count_T)/count_nonTemplateNucleotide
							count_C_percentage = float(count_C)/count_nonTemplateNucleotide
							count_G_percentage = float(count_G)/count_nonTemplateNucleotide
						else:
							count_A_percentage = 0
							count_T_percentage = 0
							count_C_percentage = 0
							count_G_percentage = 0
					elif templateNucleotide == 'T':
						if count_nonTemplateNucleotide != 0:
							count_A_percentage = float(count_A)/count_nonTemplateNucleotide
							count_T_percentage = 0
							count_C_percentage = float(count_C)/count_nonTemplateNucleotide
							count_G_percentage = float(count_G)/count_nonTemplateNucleotide
						else:
							count_A_percentage = 0
							count_T_percentage = 0
							count_C_percentage = 0
							count_G_percentage = 0
					elif templateNucleotide == 'C':
						if count_nonTemplateNucleotide != 0:
							count_A_percentage = float(count_A)/count_nonTemplateNucleotide
							count_T_percentage = float(count_T)/count_nonTemplateNucleotide
							count_C_percentage = 0
							count_G_percentage = float(count_G)/count_nonTemplateNucleotide
						else:
							count_A_percentage = 0
							count_T_percentage = 0
							count_C_percentage = 0
							count_G_percentage = 0
					elif templateNucleotide == 'G':
						if count_nonTemplateNucleotide != 0:
							count_A_percentage = float(count_A)/count_nonTemplateNucleotide
							count_T_percentage = float(count_T)/count_nonTemplateNucleotide
							count_C_percentage = float(count_C)/count_nonTemplateNucleotide
							count_G_percentage = 0
						else:
							count_A_percentage = 0
							count_T_percentage = 0
							count_C_percentage = 0
							count_G_percentage = 0
					else:
						pass
					
				else:
					count_TemplateNucleotide_percentage = 0
					count_nonTemplateNucleotide_percentage = 0
					count_A_percentage = 0
					count_T_percentage = 0
					count_C_percentage = 0
					count_G_percentage = 0

				neclueotideCountList = neclueotideCountList + [templateNucleotide, count_TemplateNucleotide_percentage, count_nonTemplateNucleotide_percentage, count_A_percentage, count_T_percentage, count_C_percentage, count_G_percentage]
				
			outList = [self.clusterName, self.clusterSeq, adjustedClusterSeq, clusterSecondSeq, templateSeq, seqCount, readCountSum, exactMatchRatio, headUnstableLength, tailUnstableLength, neclueotidePositionList, neclueotideCountList, True]
		else:
			outList = [self.clusterName, self.clusterSeq, None, None, None, None, None, None, None, None, neclueotidePositionList, None, False]
		return outList

	def calculateFeature_old_version(self):
		alignSeqList2, exactMatchCount2, clusterSeq, clusterSeqRatioList, majorSeq, majorSeqRationList, headStartPosition, tailStartPosition = self.locateStartPosition()
		#alignSeqList2, exactMatchCount2 = self.align2Standard()
		#for subseq in alignSeqList2:
		#	print subseq
		#print ReadCluster.headShiftStep
		#print ReadCluster.tailShiftStep
		readCountSum = sum(self.readCountList)
		seqCount = len(self.readSeqList)
		exactMatchRatio = float(exactMatchCount2)/readCountSum
		headAdd = alignSeqList2[0].index(self.clusterSeq)
		tailAdd = len(alignSeqList2[0])-headAdd-len(self.clusterSeq)
		#print tailAdd
		if headAdd < ReadCluster.headShiftStep:
			for i in range(len(alignSeqList2)):
				alignSeqList2[i] = (ReadCluster.headShiftStep-headAdd)*'-'+alignSeqList2[i]
		else:
			pass
		if tailAdd < ReadCluster.tailShiftStep:
			for i in range(len(alignSeqList2)):
				alignSeqList2[i] = alignSeqList2[i]+(ReadCluster.tailShiftStep-tailAdd)*'-'
		else:
			pass
		start = alignSeqList2[0].index(self.clusterSeq)
		end = alignSeqList2[0].index(self.clusterSeq)+len(self.clusterSeq)-1
		colList = list(range(start-ReadCluster.headShiftStep, start+ReadCluster.headShiftStep+1)) + list(range(end-ReadCluster.tailShiftStep, end+ReadCluster.tailShiftStep+1))
		transformedListTmp = []
		for i in range(len(colList)):
			for j in range(1, len(alignSeqList2)):
				transformedListTmp.append(alignSeqList2[j][colList[i]])
		transformedList = list(chunks(transformedListTmp, len(alignSeqList2)-1))
		neclueotideCountList = []
		neclueotidePositionList = []
		for i in range(len(transformedList)):
			item = transformedList[i]
			count_A = 0
			count_T = 0
			count_C = 0
			count_G = 0
			for j in range(len(item)):
				if item[j] == 'A':
					count_A = count_A + self.readCountList[j]
				elif item[j] == 'T':
					count_T = count_T + self.readCountList[j]
				elif item[j] == 'C':
					count_C = count_C + self.readCountList[j]
				elif item[j] == 'G':
					count_G = count_G + self.readCountList[j]
				else:
					pass
			neclueotideCountList = neclueotideCountList + [count_A,count_T,count_C,count_G]
			#neclueotideCountList.append(count_A,count_T,count_C,count_G)
			if i <= ReadCluster.headShiftStep-1:
				postionLable = 'head_minus%d_'%(ReadCluster.headShiftStep-i)
			elif i == ReadCluster.headShiftStep:
				postionLable = 'head_0_'
			elif i > ReadCluster.headShiftStep and i <= ReadCluster.headShiftStep*2:
				postionLable = 'head_plus%d_'%(i-ReadCluster.headShiftStep)
			elif i > ReadCluster.headShiftStep*2 and i < ReadCluster.headShiftStep*2+ReadCluster.tailShiftStep+1:
				postionLable = 'tail_minus%d_'%(ReadCluster.headShiftStep*2+ReadCluster.tailShiftStep+1-i)
			elif i == ReadCluster.headShiftStep*2+ReadCluster.tailShiftStep+1:
				postionLable = 'tail_0_'
			elif i > ReadCluster.headShiftStep*2+ReadCluster.tailShiftStep+1:
				postionLable = 'tail_plus%d_'%(i-ReadCluster.headShiftStep*2-ReadCluster.tailShiftStep-1)
			else:
				pass
			for neclueotide in ['A', 'T', 'C', 'G']:
				neclueotidePositionList.append(postionLable+neclueotide)
		return (self.clusterName, self.clusterSeq, seqCount, readCountSum, exactMatchRatio, neclueotidePositionList, neclueotideCountList)
