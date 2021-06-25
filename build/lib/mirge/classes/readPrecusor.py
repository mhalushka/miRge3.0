import mirge.forgi.graph.bulge_graph as fgb
import os
import re

def isisParenthese(i):
	flag = False
	if i == '(' or i == ')':
		flag = True
	return flag

#get miRNA position on RNASequence, 0-based
def getMiRNAPosition(miRNASeq, RNASequence):
	startPos = RNASequence.find(miRNASeq)
	endPos = startPos + len(miRNASeq)
	return [startPos, endPos]

#get miRNA structure
def getMiRNAStructure(miRNASeq, RNASequence, RNAStructure):
	[startPos, endPos] = getMiRNAPosition(miRNASeq, RNASequence)
	mirnaStruct = RNAStructure[startPos:endPos]
	return mirnaStruct

#get GC content of a sequence
def getGCContent(sequence):
	gc = sequence.count('G')
	gc = gc + sequence.count('g')
	gc = gc + sequence.count('C')
	gc = gc + sequence.count('c')
	GCContent = float(gc)/len(sequence)
	return GCContent

#check whether there is a loop in mature miRNA
def getMirnaIncludedInLoop( RNASequence, RNAStructure, miRNASeq):
	flag = False
	mirnaStruct = getMiRNAStructure(miRNASeq, RNASequence, RNAStructure)
	#bg2 = fgb.BulgeGraph()
	#fa2 = '\n'.join(['>miRNAStr', miRNASeq, mirnaStruct])
	#bg2.from_fasta(fa2)
	#for i in range(1, len(miRNASeq)+1):
	#	if bg2.pairing_partner(i)
	pattern = re.compile("\(.*\)")
	pattern.search(mirnaStruct)
	if(('(' in mirnaStruct) and (')' in mirnaStruct)):
		if pattern.search(mirnaStruct) is None:
			pass
		else:
			flag = True
	return flag

#check miRNA located in which arm
def checkArm(RNASequence, RNAStructure, miRNASeq):
	armDetailedType = ''
	if(getMirnaIncludedInLoop(RNASequence, RNAStructure, miRNASeq)):
		armDetailedType = 'loop'
	#check arms
	else:
		mirnaStruct = getMiRNAStructure(miRNASeq, RNASequence, RNAStructure)
		if '(' in mirnaStruct :
			if mirnaStruct.count('(') >= mirnaStruct.count(')'):
				armDetailedType = 'arm5'
			else:
				armDetailedType = 'arm3'
		elif ')' in mirnaStruct :
			if mirnaStruct.count(')') >= mirnaStruct.count('('):
				armDetailedType = 'arm3'
			else:
				armDetailedType = 'arm5'
		else:
			armDetailedType = 'unmatchedRegion'
	return armDetailedType
#end fun

def decorate_str(RNAStructue):
	open_brackets1 = '('
	close_brackets1 = ')'
	brackets_map1 = {')': '('}
	stack1 = []
	errorFrontPosition1 = []
	for i, char in enumerate(RNAStructue):
		if char in open_brackets1:
			stack1.append(char)
		elif char in close_brackets1:
			if len(stack1) < 1:
				stack1.append(char)
				errorFrontPosition1.append(i)
			elif brackets_map1[char] == stack1[-1]:
				stack1.pop()
			else:
				stack1.append(char)
				errorFrontPosition1.append(i)
		else:
			continue
	open_brackets2 = ')'
	close_brackets2 = '('
	brackets_map2 = {'(': ')'}
	stack2 = []
	errorFrontPosition2 = []
	for i in range(len(RNAStructue)-1, -1, -1):
		char = RNAStructue[i]
		if char in open_brackets2:
			stack2.append(char)
		elif char in close_brackets2:
			if len(stack2) < 1:
				stack2.append(char)
				errorFrontPosition2.append(i)
			elif brackets_map2[char] == stack2[-1]:
				stack2.pop()
			else:
				stack2.append(char)
				errorFrontPosition2.append(i)
		else:
			continue
	errorFrontPosition = errorFrontPosition1+errorFrontPosition2
	decorated_RNAStr = ''
	for i, char in enumerate(RNAStructue):
		if i in errorFrontPosition:
			decorated_RNAStr = decorated_RNAStr + '.'
		else:
			decorated_RNAStr = decorated_RNAStr + char
	return decorated_RNAStr

# Statisitcally analyze the list such as: ['s0', 's0', 'i0', 'i0','s1', 's1'...], location is 0-based
def statisList(elementIndexList, RNASeq):
# Here , the length of elementIndexList and RNASeq must be equal.
	elementSeq = []
	elementDic = {}
	i = 0
	for item in elementIndexList:
		if item not in elementSeq:
			elementSeq.append(item)
			elementDic.update({item:[]})
		i = i + 1

# Rewrote function of 'to_element_string', because the lower version forgi 0.2 lack the option of with_numbers=True
def to_element_string(bulgeGraphInstance, with_numbers=False):
	output_str = [' '] * (bulgeGraphInstance.seq_length + 1)
	output_nr = [' '] * (bulgeGraphInstance.seq_length + 1)
	for d in bulgeGraphInstance.defines.keys():
		for resi in bulgeGraphInstance.define_residue_num_iterator(d, adjacent=False):
			output_str[resi] = d[0]
			output_nr[resi] = d[-1]
	if with_numbers:
		return "".join(output_str).strip()+"\n"+"".join(output_nr).strip()
	else:
		return "".join(output_str).strip()

def chunks(arr, n):
	return [arr[i:i+n] for i in range(0, len(arr), n)]

class ReadPrecusor(object):
	def __init__(self, rnafoldCmdTmp, premiRNASeqName, premiRNASeq, premiRNAStructure, mfe, tmpdir, prunedLength, *miRNASeq):
		# The first element in the list of *miRNASeq must be the miRNAs in the pre-MiRNA.
		# The features about the miRNA should be about the first element.
		# If *miRNASeq has two elements, the second one will be the other one which is only used to generate the corePrecusor.
		self.rnafoldCmdTmp = rnafoldCmdTmp
		self.premiRNASeqName = premiRNASeqName
		self.premiRNASeq = premiRNASeq
		#self.direction = premiRNASeq
		self.premiRNAStructure = premiRNAStructure
		self.mfe = mfe
		self.tmpdir = tmpdir
		self.miRNASeq = miRNASeq
		self.prunedLength = int(prunedLength)
		#Here, self.miRNASeq is a tuple which stores one or two mirNA sequences.

	def getPrunedPrecusor(self):
		if len(self.miRNASeq) == 1:
			miRNASeq1 = self.miRNASeq[0]
			selectedMiRNASeq = miRNASeq1
		else:
			miRNASeq1 = self.miRNASeq[0]
			miRNASeq2 = self.miRNASeq[1]
			startPos1, endPos1 = getMiRNAPosition(miRNASeq1, self.premiRNASeq)
			if startPos1 == 20:
				selectedMiRNASeq = miRNASeq1
			else:
				selectedMiRNASeq = miRNASeq2
		armType = checkArm(self.premiRNASeq, self.premiRNAStructure, selectedMiRNASeq)
		bg1 = fgb.BulgeGraph()
		fa = '\n'.join(['>PrecusorStr', self.premiRNASeq, self.premiRNAStructure])
		bg1.from_fasta(fa)
		if armType == 'arm5':
			startPos1, endPos1 = getMiRNAPosition(selectedMiRNASeq, self.premiRNASeq)
			otherCount = 0 
			for i in range(startPos1, endPos1):
				if self.premiRNAStructure[i] == '(':
					break
				otherCount = otherCount + 1
			# Locate the partner of a base pair
			try:
				pairedPos = bg1.pairing_partner(startPos1+otherCount+1) - 1
				if otherCount + pairedPos >= len(self.premiRNASeq):
					newEndPos = len(self.premiRNASeq) - 1
				else:
					newEndPos = otherCount + pairedPos
				prunedPrecursorSeq = self.premiRNASeq[startPos1-self.prunedLength : newEndPos+1+self.prunedLength]
				prunedPrecusorStr = self.premiRNAStructure[startPos1-self.prunedLength : newEndPos+1+self.prunedLength]
			except TypeError:
				prunedPrecursorSeq = None
				prunedPrecusorStr = None
		else:
			startPos1, endPos1 = getMiRNAPosition(selectedMiRNASeq, self.premiRNASeq)
			otherCount = 0 
			for i in range(endPos1-1, -1, -1):
				if self.premiRNAStructure[i] == ')':
					break
				otherCount = otherCount + 1
			# Locate the partner of a base pair
			try:
				pairedPos = bg1.pairing_partner(endPos1-otherCount) - 1
				if pairedPos - otherCount < 0:
					newStartPos = 0
				else:
					newStartPos = pairedPos - otherCount
				prunedPrecursorSeq = self.premiRNASeq[newStartPos-self.prunedLength : endPos1+self.prunedLength]
				prunedPrecusorStr = self.premiRNAStructure[newStartPos-self.prunedLength : endPos1+self.prunedLength]
			except TypeError:
				prunedPrecursorSeq = None
				prunedPrecusorStr = None
		# Decorate the CorePrecusor structure.
		if not (prunedPrecusorStr is None):
			prunedPrecusorStr = decorate_str(prunedPrecusorStr)
			
		# Predict the new Structure for the prunedPrecursorSeq
		if not (prunedPrecursorSeq is None):
			with open(os.path.join(self.tmpdir, 'prunedPrecusor.fa'), 'w') as outf:
				fa_tmp = '\n'.join(['>'+self.premiRNASeqName, prunedPrecursorSeq])
				outf.write(fa_tmp+'\n')
			f1 = os.path.join(self.tmpdir, 'prunedPrecusor.fa')
			f2 = os.path.join(self.tmpdir, 'prunedPrecusor.str')
			os.system('%s < %s --noPS --noLP > %s'%(self.rnafoldCmdTmp, f1, f2))
			try:	
				with open(f2, 'r') as inf:
					try:
						tmp = inf.readlines()
						prunedPrecusorStrNew = tmp[2].strip().split()[0]
						if len(tmp[2].strip().split()[1]) >= 3:
							mFE = float(tmp[2].strip().split()[1][1:-1].strip())
						else:
							mFE = float(tmp[2].strip().split()[2][:-1].strip())
					except:
						mFE = 0
						prunedPrecusorStrNew = prunedPrecusorStr
			except FileNotFoundError:
				prunedPrecusorStrNew = None
			os.remove(f1)
			os.remove(f2)
		else:
			prunedPrecusorStrNew = None
			mFE = None
		return (prunedPrecursorSeq, prunedPrecusorStrNew, mFE)

	def featuresInPrecusor(self):
		if not (self.premiRNAStructure  is None) and '(' in self.premiRNAStructure  and ')' in self.premiRNAStructure:
			bg = fgb.BulgeGraph()
			fa = '\n'.join(['>prunedPrecursorStr', self.premiRNASeq, self.premiRNAStructure])
			bg.from_fasta(fa)
			bindingCount1 = 0
			bindingCount2 = 0
			for item in self.premiRNAStructure:
				if item == '(' :
					bindingCount1 = bindingCount1 + 1
				if item == ')' :
					bindingCount2 = bindingCount2 + 1
			# Calculate bindingCount in the CorePrecusor
			bindingCount = max([bindingCount1, bindingCount2])
			# Calculate the hairpinCount and interiorLoopCount
			hairpinCount = 0
			interiorLoopCount = 0
			for string in bg.to_bg_string().strip().split('\n'):
				if 'define' in string:
					contentTmp = string.split(' ')
					elementType = contentTmp[1]
					if elementType[0] == 'h':
						hairpinCount = hairpinCount + 1
					if elementType[0] == 'i':
						interiorLoopCount = interiorLoopCount + 1
			# Calculate apicalLoopSize of the CorePrecusor
			# To get compatible with the older version of forgi, the below statement is abandoned.
			#elementTypeSeq, indeSeq = bg.to_element_string(with_numbers=True).split('\n')
			elementTypeSeq, indeSeq = to_element_string(bg, with_numbers=True).split('\n')
			
			elementIndexList = [elementTypeSeq[i]+indeSeq[i] for i in range(len(elementTypeSeq))]
			elementSeqList = []
			elementlocationList = []
			# the position is 0-based
			for index, item in enumerate(elementIndexList):
				if len(elementlocationList) == 0:
					elementSeqList.append(item)
					elementlocationList.append([index, index])
				elif item != elementSeqList[-1]:
					elementSeqList.append(item)
					elementlocationList.append([index, index])
				else:
					elementlocationList[-1][1] = index
			elementDic = {}
			for i in range(len(elementSeqList)):
				if elementSeqList[i] not in elementDic.keys():
					elementDic.update({elementSeqList[i] : [elementlocationList[i]]})
				else:
					elementDic[elementSeqList[i]].append(elementlocationList[i])
			hairpinRagion = []
			for d in elementDic.keys():
				if d[0] == 'h':
					hairpinRagion = hairpinRagion + elementDic[d]
			if len(hairpinRagion) >= 1:
				apicalLoopSize = max([len(range(item[0], item[1]+1)) for item in hairpinRagion])
			else:
				apicalLoopSize = 0
			armTypeOverlenDic = {}
			# Calculate the armType of the miRNA
			armType = checkArm(self.premiRNASeq, self.premiRNAStructure, self.miRNASeq[0])
			# Calculate the overlap of miRNA with the hairpin
			# startPos and endPos are the locations of the miRNA in the corePrecusor. 0-based
			startPos, endPos = getMiRNAPosition(self.miRNASeq[0], self.premiRNASeq)
			endPos = endPos - 1
			if len(hairpinRagion) >= 1:
				overlapLen = max([len(set(range(startPos, endPos+1)) & set(range(item[0], item[1]+1))) for item in hairpinRagion])
			else:
				overlapLen = 0
			if overlapLen == 0:
				if armType == 'arm5':
					distanceToLoop = min([item[0]-1-endPos for item in hairpinRagion])
				else:
					distanceToLoop = min([startPos-item[1]-1 for item in hairpinRagion])
			else:
				for item in hairpinRagion:
					if len(set(range(startPos, endPos+1)) & set(range(item[0], item[1]+1))) == overlapLen:
						if armType == 'arm5':
							distanceToLoop = item[0]-1-endPos
						else:
							distanceToLoop = startPos - item[1] - 1
						break
			armTypeOverlenDic.update({armType: [overlapLen, distanceToLoop]})
			if len(self.miRNASeq) == 2:
				armType2 = checkArm(self.premiRNASeq, self.premiRNAStructure, self.miRNASeq[1])
				startPos2, endPos2 = getMiRNAPosition(self.miRNASeq[1], self.premiRNASeq)
				endPos2 = endPos2 - 1
				if len(hairpinRagion) >= 1:
					overlapLen2 = max([len(set(range(startPos2, endPos2+1)) & set(range(item[0], item[1]+1))) for item in hairpinRagion])
				else:
					overlapLen2 = 0
				if overlapLen2 == 0:
					if armType2 == 'arm5':
						distanceToLoop2 = min([item[0]-1-endPos2 for item in hairpinRagion])
					else:
						distanceToLoop2 = min([startPos2-item[1]-1 for item in hairpinRagion])
				else:
					for item in hairpinRagion:
						if len(set(range(startPos2, endPos2+1)) & set(range(item[0], item[1]+1))) == overlapLen2:
							if armType2 == 'arm5':
								distanceToLoop2 = item[0]-1-endPos2
							else:
								distanceToLoop2 = startPos2 - item[1]-1
							break
				armTypeOverlenDic.update({armType2: [overlapLen2, distanceToLoop2]})
			# Calculate the stem length of precusor.
			stemLen = 0
			for d in elementDic.keys():
				if d[0] == 's':
					stemLen = stemLen + max([item[1]+1-item[0] for item in elementDic[d]])
			# Check whether apical UGU/UGUG motif exist in hairpin region
			flag = 'No'
			for d in elementDic.keys():
				if d[0] == 'h':
					for item in elementDic[d]:
						seqTmp = self.premiRNASeq[item[0]:item[1]]
						if 'UGU' in seqTmp or 'UGUG' in seqTmp:
							flag = 'Yes'
							break
			# pair state
			if len(self.miRNASeq) == 2:
				pairState = 'Yes'
			else:
				pairState = 'No'
		else:
			hairpinCount = None
			bindingCount = None
			interiorLoopCount = None
			armType = None
			apicalLoopSize = None
			armTypeOverlenDic = None
			stemLen= None
			flag= None
			pairState = None
		return (hairpinCount, bindingCount, interiorLoopCount, armType, apicalLoopSize, armTypeOverlenDic, stemLen, flag, self.mfe, pairState)

	def getPrunedPrecusor_old(self):
		if len(self.miRNASeq) == 1:
			miRNASeq1 = self.miRNASeq[0]
			selectedMiRNASeq = miRNASeq1
		else:
			miRNASeq1 = self.miRNASeq[0]
			miRNASeq2 = self.miRNASeq[1]
			startPos1, endPos1 = getMiRNAPosition(miRNASeq1, self.premiRNASeq)
			if startPos1 == 20:
				selectedMiRNASeq = miRNASeq1
			else:
				selectedMiRNASeq = miRNASeq2
		armType = checkArm(self.premiRNASeq, self.premiRNAStructure, selectedMiRNASeq)
		bg1 = fgb.BulgeGraph()
		fa = '\n'.join(['>PrecusorStr', self.premiRNASeq, self.premiRNAStructure])
		bg1.from_fasta(fa)
		if armType == 'arm5':
			startPos1, endPos1 = getMiRNAPosition(selectedMiRNASeq, self.premiRNASeq)
			startPos1 = 0
			otherCount = 0 
			for i in range(startPos1, endPos1):
				if self.premiRNAStructure[i] == '(':
					break
				otherCount = otherCount + 1
			# Locate the partner of a base pair
			try:
				pairedPos = bg1.pairing_partner(startPos1+otherCount+1) - 1
				if otherCount + pairedPos >= len(self.premiRNASeq):
					newEndPos = len(self.premiRNASeq) - 1
				else:
					newEndPos = otherCount + pairedPos
				prunedPrecusorSeq = self.premiRNASeq[startPos1 : newEndPos+1]
				prunedPrecusorStr = self.premiRNAStructure[startPos1 : newEndPos+1]
			except TypeError:
				prunedPrecusorSeq = None
				prunedPrecusorStr = None
		else:
			startPos1, endPos1 = getMiRNAPosition(selectedMiRNASeq, self.premiRNASeq)
			endPos1 = len(self.premiRNASeq)
			print(startPos1)
			print(endPos1)
			otherCount = 0 
			for i in range(endPos1-1, -1, -1):
				if self.premiRNAStructure[i] == ')':
					break
				otherCount = otherCount + 1
			# Locate the partner of a base pair
			print(otherCount)
			try:
				pairedPos = bg1.pairing_partner(endPos1-otherCount) - 1
				if pairedPos - otherCount < 0:
					newStartPos = 0
				else:
					newStartPos = pairedPos - otherCount
				prunedPrecusorSeq = self.premiRNASeq[newStartPos : endPos1]
				prunedPrecusorStr = self.premiRNAStructure[newStartPos : endPos1]
			except TypeError:
				prunedPrecusorSeq = None
				prunedPrecusorStr = None
		# Decorate the CorePrecusor structure.
		if not (prunedPrecusorStr is None):
			prunedPrecusorStr = decorate_str(prunedPrecusorStr)
		return (prunedPrecusorSeq, prunedPrecusorStr)

	def featuresInPrunedPrecusor(self):
		prunedPrecusorSeq, prunedPrecusorStr = self.getPrunedPrecusor()
		if not (prunedPrecusorStr is None) and '(' in prunedPrecusorStr and ')' in prunedPrecusorStr:
			bg = fgb.BulgeGraph()
			fa = '\n'.join(['>prunedPrecusorStr', prunedPrecusorStr, prunedPrecusorStr])
			bg.from_fasta(fa)
			bindingCount1 = 0
			bindingCount2 = 0
			for item in prunedPrecusorStr:
				if item == '(' :
					bindingCount1 = bindingCount1 + 1
				if item == ')' :
					bindingCount2 = bindingCount2 + 1
			# Calculate bindingCount in the CorePrecusor
			bindingCount = max([bindingCount1, bindingCount2])
			# Calculate the hairpinCount and interiorLoopCount
			hairpinCount = 0
			interiorLoopCount = 0
			for string in bg.to_bg_string().strip().split('\n'):
				if 'define' in string:
					contentTmp = string.split(' ')
					elementType = contentTmp[1]
					if elementType[0] == 'h':
						hairpinCount = hairpinCount + 1
					if elementType[0] == 'i':
						interiorLoopCount = interiorLoopCount + 1
			# Calculate apicalLoopSize of the CorePrecusor
			# To get compatible with the older version of forgi, the below statement is abandoned.
			#elementTypeSeq, indeSeq = bg.to_element_string(with_numbers=True).split('\n')
			elementTypeSeq, indeSeq = to_element_string(bg, with_numbers=True).split('\n')
			
			elementIndexList = [elementTypeSeq[i]+indeSeq[i] for i in range(len(elementTypeSeq))]
			elementSeqList = []
			elementlocationList = []
			# the position is 0-based
			for index, item in enumerate(elementIndexList):
				if len(elementlocationList) == 0:
					elementSeqList.append(item)
					elementlocationList.append([index, index])
				elif item != elementSeqList[-1]:
					elementSeqList.append(item)
					elementlocationList.append([index, index])
				else:
					elementlocationList[-1][1] = index
			elementDic = {}
			for i in range(len(elementSeqList)):
				if elementSeqList[i] not in elementDic.keys():
					elementDic.update({elementSeqList[i] : [elementlocationList[i]]})
				else:
					elementDic[elementSeqList[i]].append(elementlocationList[i])
			hairpinRagion = []
			for d in elementDic.keys():
				if d[0] == 'h':
					hairpinRagion = hairpinRagion + elementDic[d]
			if len(hairpinRagion) >= 1:
				apicalLoopSize = max([len(range(item[0], item[1]+1)) for item in hairpinRagion])
			else:
				apicalLoopSize = 0
			armTypeOverlenDic = {}
			# Calculate the armType of the miRNA
			armType = checkArm(self.premiRNASeq, self.premiRNAStructure, self.miRNASeq[0])
			# Calculate the overlap of miRNA with the hairpin
			# startPos and endPos are the locations of the miRNA in the corePrecusor. 0-based
			startPos, endPos = getMiRNAPosition(self.miRNASeq[0], prunedPrecusorSeq)
			endPos = endPos - 1
			if len(hairpinRagion) >= 1:
				overlapLen = max([len(set(range(startPos, endPos+1)) & set(range(item[0], item[1]+1))) for item in hairpinRagion])
			else:
				overlapLen = 0
			if overlapLen == 0:
				if armType == 'arm5':
					distanceToLoop = min([item[0]-1-endPos for item in hairpinRagion])
				else:
					distanceToLoop = min([startPos-item[1]-1 for item in hairpinRagion])
			else:
				for item in hairpinRagion:
					if len(set(range(startPos, endPos+1)) & set(range(item[0], item[1]+1))) == overlapLen:
						if armType == 'arm5':
							distanceToLoop = item[0]-1-endPos
						else:
							distanceToLoop = startPos - item[1] - 1
						break
			armTypeOverlenDic.update({armType: [overlapLen, distanceToLoop]})
			if len(self.miRNASeq) == 2:
				armType2 = checkArm(self.premiRNASeq, self.premiRNAStructure, self.miRNASeq[1])
				startPos2, endPos2 = getMiRNAPosition(self.miRNASeq[1], prunedPrecusorSeq)
				endPos2 = endPos2 - 1
				if len(hairpinRagion) >= 1:
					overlapLen2 = max([len(set(range(startPos2, endPos2+1)) & set(range(item[0], item[1]+1))) for item in hairpinRagion])
				else:
					overlapLen2 = 0
				if overlapLen2 == 0:
					if armType2 == 'arm5':
						distanceToLoop2 = min([item[0]-1-endPos2 for item in hairpinRagion])
					else:
						distanceToLoop2 = min([startPos2-item[1]-1 for item in hairpinRagion])
				else:
					for item in hairpinRagion:
						if len(set(range(startPos2, endPos2+1)) & set(range(item[0], item[1]+1))) == overlapLen2:
							if armType2 == 'arm5':
								distanceToLoop2 = item[0]-1-endPos2
							else:
								distanceToLoop2 = startPos2 - item[1]-1
							break
				armTypeOverlenDic.update({armType2: [overlapLen2, distanceToLoop2]})
			# Calculate the stem length of precusor.
			stemLen = 0
			for d in elementDic.keys():
				if d[0] == 's':
					stemLen = stemLen + max([item[1]+1-item[0] for item in elementDic[d]])
			# Check whether apical UGU/UGUG motif exist in hairpin region
			flag = 'No'
			for d in elementDic.keys():
				if d[0] == 'h':
					for item in elementDic[d]:
						seqTmp = prunedPrecusorSeq[item[0]:item[1]]
						if 'UGU' in seqTmp or 'UGUG' in seqTmp:
							flag = 'Yes'
							break
			# minimun free energy
			with open(os.path.join(self.tmpdir, 'prunedPrecusor.fa'), 'w') as outf:
				fa_tmp = '\n'.join(['>'+self.premiRNASeqName, prunedPrecusorSeq, prunedPrecusorStr])
				outf.write(fa_tmp+'\n')
			f1 = os.path.join(self.tmpdir, 'prunedPrecusor.fa')
			f2 = os.path.join(self.tmpdir, 'prunedPrecusor.str')
			os.system('%s < %s --noPS --noLP > %s'%(self.rnafoldCmdTmp, f1, f2))
			with open(f2, 'r') as inf:
				try:
					tmp = inf.readlines()
					if len(tmp[-1].strip().split(' ')[1]) >= 3:
						mFE = float(tmp[-1].strip().split(' ')[1][1:-1].strip())
					else:
						mFE = float(tmp[-1].strip().split(' ')[2][:-1].strip())
				except:
					mFE = 0
			os.remove(f1)
			os.remove(f2)
			# pair state
			if len(self.miRNASeq) == 2:
				pairState = 'Yes'
			else:
				pairState = 'No'
		else:
			hairpinCount = None
			bindingCount = None
			interiorLoopCount = None
			armType = None
			apicalLoopSize = None
			armTypeOverlenDic = None
			stemLen= None
			flag= None
			mFE = None
			pairState = None
		return (hairpinCount, bindingCount, interiorLoopCount, armType, apicalLoopSize, armTypeOverlenDic, stemLen, flag, mFE, pairState)

	def percentageOfPairedInMiRNA(self):
		miRNAStr = getMiRNAStructure(self.miRNASeq[0], self.premiRNASeq, self.premiRNAStructure)
		try:
			percentage = (miRNAStr.count('(')+miRNAStr.count(')'))/float(len(miRNAStr))
		except:
			percentage = 0
		return percentage

	def bindingsInMiRNA(self):
		miRNAStr = getMiRNAStructure(self.miRNASeq[0], self.premiRNASeq, self.premiRNAStructure)
		countBindingsInMiRNA = miRNAStr.count('(')+miRNAStr.count(')')
		return countBindingsInMiRNA

	def getCorePrecusor(self):
		if len(self.miRNASeq) == 1:
			miRNASeq1 = self.miRNASeq[0]
			armType = checkArm(self.premiRNASeq, self.premiRNAStructure, miRNASeq1)
			bg1 = fgb.BulgeGraph()
			fa = '\n'.join(['>PrecusorStr', self.premiRNASeq, self.premiRNAStructure])
			bg1.from_fasta(fa)
			if armType == 'arm5':
				startPos1, endPos1 = getMiRNAPosition(miRNASeq1, self.premiRNASeq)
				otherCount = 0 
				for i in range(startPos1, endPos1):
					if self.premiRNAStructure[i] == '(':
						break
					otherCount = otherCount + 1
				# Locate the partner of a base pair
				try:
					pairedPos = bg1.pairing_partner(startPos1+otherCount+1) - 1
					if otherCount + pairedPos >= len(self.premiRNASeq):
						newEndPos = len(self.premiRNASeq) - 1
					else:
						newEndPos = otherCount + pairedPos
					corePrecusorSeq = self.premiRNASeq[startPos1 : newEndPos+1]
					corePrecusorStr = self.premiRNAStructure[startPos1 : newEndPos+1]
				except TypeError:
					corePrecusorSeq = None
					corePrecusorStr = None
			else:
				startPos1, endPos1 = getMiRNAPosition(miRNASeq1, self.premiRNASeq)
				otherCount = 0 
				for i in range(endPos1-1, -1, -1):
					if self.premiRNAStructure[i] == ')':
						break
					otherCount = otherCount + 1
				# Locate the partner of a base pair
				try:
					pairedPos = bg1.pairing_partner(endPos1-otherCount) - 1
					if pairedPos - otherCount < 0:
						newStartPos = 0
					else:
						newStartPos = pairedPos - otherCount
					corePrecusorSeq = self.premiRNASeq[newStartPos : endPos1]
					corePrecusorStr = self.premiRNAStructure[newStartPos : endPos1]
				except TypeError:
					corePrecusorSeq = None
					corePrecusorStr = None
		else:
			miRNASeq1 = self.miRNASeq[0]
			miRNASeq2 = self.miRNASeq[1]
			startPos1, endPos1 = getMiRNAPosition(miRNASeq1, self.premiRNASeq)
			startPos2, endPos2 = getMiRNAPosition(miRNASeq2, self.premiRNASeq)
			newStartPos = min([startPos1, endPos1, startPos2, endPos2])
			newEndPos = max([startPos1, endPos1, startPos2, endPos2])
			corePrecusorSeq = self.premiRNASeq[newStartPos : newEndPos]
			corePrecusorStr = self.premiRNAStructure[newStartPos : newEndPos]
		# Decorate the CorePrecusor structure.
		if not (corePrecusorStr is None):
			corePrecusorStr = decorate_str(corePrecusorStr)
		return (corePrecusorSeq, corePrecusorStr)

	def featuresInCorePrecusor(self):
		corePrecusorSeq, corePrecusorStr = self.getCorePrecusor()
		if not (corePrecusorStr is None) and '(' in corePrecusorStr and ')' in corePrecusorStr:
			bg = fgb.BulgeGraph()
			fa = '\n'.join(['>corePrecusorStr', corePrecusorSeq, corePrecusorStr])
			bg.from_fasta(fa)
			bindingCount1 = 0
			bindingCount2 = 0
			for item in corePrecusorStr:
				if item == '(' :
					bindingCount1 = bindingCount1 + 1
				if item == ')' :
					bindingCount2 = bindingCount2 + 1
			# Calculate bindingCount in the CorePrecusor
			bindingCount = max([bindingCount1, bindingCount2])
			# Calculate the hairpinCount and interiorLoopCount
			hairpinCount = 0
			interiorLoopCount = 0
			for string in bg.to_bg_string().strip().split('\n'):
				if 'define' in string:
					contentTmp = string.split(' ')
					elementType = contentTmp[1]
					if elementType[0] == 'h':
						hairpinCount = hairpinCount + 1
					if elementType[0] == 'i':
						interiorLoopCount = interiorLoopCount + 1
			# Calculate apicalLoopSize of the CorePrecusor
			# To get compatible with the older version of forgi, the below statement is abandoned.
			#elementTypeSeq, indeSeq = bg.to_element_string(with_numbers=True).split('\n')
			elementTypeSeq, indeSeq = to_element_string(bg, with_numbers=True).split('\n')
			
			elementIndexList = [elementTypeSeq[i]+indeSeq[i] for i in range(len(elementTypeSeq))]
			elementSeqList = []
			elementlocationList = []
			# the position is 0-based
			for index, item in enumerate(elementIndexList):
				if len(elementlocationList) == 0:
					elementSeqList.append(item)
					elementlocationList.append([index, index])
				elif item != elementSeqList[-1]:
					elementSeqList.append(item)
					elementlocationList.append([index, index])
				else:
					elementlocationList[-1][1] = index
			elementDic = {}
			for i in range(len(elementSeqList)):
				if elementSeqList[i] not in elementDic.keys():
					elementDic.update({elementSeqList[i] : [elementlocationList[i]]})
				else:
					elementDic[elementSeqList[i]].append(elementlocationList[i])
			hairpinRagion = []
			for d in elementDic.keys():
				if d[0] == 'h':
					hairpinRagion = hairpinRagion + elementDic[d]
			if len(hairpinRagion) >= 1:
				apicalLoopSize = max([len(range(item[0], item[1]+1)) for item in hairpinRagion])
			else:
				apicalLoopSize = 0
			armTypeOverlenDic = {}
			# Calculate the armType of the miRNA
			armType = checkArm(self.premiRNASeq, self.premiRNAStructure, self.miRNASeq[0])
			# Calculate the overlap of miRNA with the hairpin
			# startPos and endPos are the locations of the miRNA in the corePrecusor. 0-based
			startPos, endPos = getMiRNAPosition(self.miRNASeq[0], corePrecusorSeq)
			endPos = endPos - 1
			if len(hairpinRagion) >= 1:
				overlapLen = max([len(set(range(startPos, endPos+1)) & set(range(item[0], item[1]+1))) for item in hairpinRagion])
			else:
				overlapLen = 0
			if overlapLen == 0:
				if armType == 'arm5':
					distanceToLoop = min([item[0]-1-endPos for item in hairpinRagion])
				else:
					distanceToLoop = min([startPos-item[1]-1 for item in hairpinRagion])
			else:
				for item in hairpinRagion:
					if len(set(range(startPos, endPos+1)) & set(range(item[0], item[1]+1))) == overlapLen:
						if armType == 'arm5':
							distanceToLoop = item[0]-1-endPos
						else:
							distanceToLoop = startPos - item[1] - 1
						break
			armTypeOverlenDic.update({armType: [overlapLen, distanceToLoop]})
			if len(self.miRNASeq) == 2:
				armType2 = checkArm(self.premiRNASeq, self.premiRNAStructure, self.miRNASeq[1])
				startPos2, endPos2 = getMiRNAPosition(self.miRNASeq[1], corePrecusorSeq)
				endPos2 = endPos2 - 1
				if len(hairpinRagion) >= 1:
					overlapLen2 = max([len(set(range(startPos2, endPos2+1)) & set(range(item[0], item[1]+1))) for item in hairpinRagion])
				else:
					overlapLen2 = 0
				if overlapLen2 == 0:
					if armType2 == 'arm5':
						distanceToLoop2 = min([item[0]-1-endPos2 for item in hairpinRagion])
					else:
						distanceToLoop2 = min([startPos2-item[1]-1 for item in hairpinRagion])
				else:
					for item in hairpinRagion:
						if len(set(range(startPos2, endPos2+1)) & set(range(item[0], item[1]+1))) == overlapLen2:
							if armType2 == 'arm5':
								distanceToLoop2 = item[0]-1-endPos2
							else:
								distanceToLoop2 = startPos2 - item[1]-1
							break
				armTypeOverlenDic.update({armType2: [overlapLen2, distanceToLoop2]})
			# Calculate the stem length of precusor.
			stemLen = 0
			for d in elementDic.keys():
				if d[0] == 's':
					stemLen = stemLen + max([item[1]+1-item[0] for item in elementDic[d]])
			# Check whether apical UGU/UGUG motif exist in hairpin region
			flag = 'No'
			for d in elementDic.keys():
				if d[0] == 'h':
					for item in elementDic[d]:
						seqTmp = corePrecusorSeq[item[0]:item[1]]
						if 'UGU' in seqTmp or 'UGUG' in seqTmp:
							flag = 'Yes'
							break
			# minimun free energy
			with open(os.path.join(self.tmpdir, 'corePrecusor.fa'), 'w') as outf:
				fa_tmp = '\n'.join(['>'+self.premiRNASeqName, corePrecusorSeq, corePrecusorStr])
				outf.write(fa_tmp+'\n')
			f1 = os.path.join(self.tmpdir, 'corePrecusor.fa')
			f2 = os.path.join(self.tmpdir, 'corePrecusor.str')
			os.system('%s < %s --noPS > %s'%(self.rnafoldCmdTmp, f1, f2))
			with open(f2, 'r') as inf:
				try:
					tmp = inf.readlines()
					if len(tmp[-1].strip().split(' ')[1]) >= 3:
						mFE = float(tmp[-1].strip().split(' ')[1][1:-1].strip())
					else:
						mFE = float(tmp[-1].strip().split(' ')[2][:-1].strip())
				except:
					mFE = 0
			os.remove(f1)
			os.remove(f2)
			# pair state
			if len(self.miRNASeq) == 2:
				pairState = 'Yes'
			else:
				pairState = 'No'
		else:
			hairpinCount = None
			bindingCount = None
			interiorLoopCount = None
			armType = None
			apicalLoopSize = None
			armTypeOverlenDic = None
			stemLen= None
			flag= None
			mFE = None
			pairState = None
		return (hairpinCount, bindingCount, interiorLoopCount, armType, apicalLoopSize, armTypeOverlenDic, stemLen, flag, mFE, pairState)
