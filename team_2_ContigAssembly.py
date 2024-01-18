'''
A program that askes for a fasta file, reads in the fragments
condenses down to just the useful frags and then proceeds with
the haplotype assembly problem and prints each haplotype and
its complement.

Let's put the puzzle back together!!

Eric Allen, Ben Julius, Jennifer Flynn, Kyle Smith
'''

from readfasta import readfasta
from Final_Semiglobal import align
from build_Consensus_Seq import buildConsensusSeq, merge
from haplotypeAssembly import haplotypeAssembly
from math import ceil

'''
Generates the pairwise matrix used in determining how to align the seqs
'''
def generateMatrix(frags):
    pairwiseMatrix = [ [0 for i in range(len(frags))] for j in range(len(frags))]
    for row in range(len(pairwiseMatrix)):
        for col in range(len(pairwiseMatrix)):
            if col > row:
                alignment = align(frags[row],frags[col])
                pairwiseMatrix[row][col] = alignment[2]/len(alignment[0])
    return pairwiseMatrix

'''
Takes a list of fragments and checks to see if any of them exist completely
in other stings. If so then condense down until only longer, more useful
strings remain.
'''
def consolidate(frags):

    comboList = [ [] for val in range(len(frags))]
    for i in range(len(frags)):
        frag1 = frags[i]
        if frag1:
            for j in range(len(frags)):
                if i != j:
                    if frag1 in frags[j]:
                        comboList[j].append(i)

    for row in range(len(comboList)):
        for val in range(len(comboList[row])):
            comboList[comboList[row][val]] = [-1]

    allFrags = []
    finalFrags = []
    for row in range(len(comboList)):
        if -1 not in comboList[row] and len(comboList[row]) > 0:
            finalFrags.append(frags[row])


    return finalFrags


'''
Returns the largest value in the matrix and is x,y position
'''
def findLargestVal(matrix):
    largestElement = [matrix[0][1],0,1]

    for row in range(len(matrix)):
        for col in range(len(matrix)):
            if col > row:
                if matrix[row][col] > largestElement[0]:
                    largestElement[0] = matrix[row][col]
                    largestElement[1] = row
                    largestElement[2] = col
    return largestElement

'''
Using the smallest value and its position we figure the new values for all of
the nonzero elements in the appropriate row and column. It then deletes the
values that are no longer needed
'''
def condense(distanceMatrix, i, j):
    
    distanceIJ = distanceMatrix[i][j]
    distanceIK = 0
    distanceJK = 0
    for k in range(len(distanceMatrix)):
        if (k != i) and (k != j):
            if j < k:
                distanceJK = distanceMatrix[k][j]
            else:
                distanceJK = distanceMatrix[j][k]
            if i < k:
                distanceIK = distanceMatrix[k][i]
                distanceMatrix[k][i] = ceil((distanceIK + distanceJK - distanceIJ)/2)
            else:
                distanceIK = distanceMatrix[i][k]
                distanceMatrix[i][k] = ceil((distanceIK + distanceJK - distanceIJ)/2)

    for row in range(len(distanceMatrix)):
        del distanceMatrix[row][j]
    del distanceMatrix[j]

    return distanceMatrix


'''
Obtains the opitimally aligned seqs
'''
def alignSeqs(origFrags):
    frags = []
    tempFrags = []
    guideTree = []
    consensusSeq = []
    finalSeqs = []

    #### Make a new clean list of the frags so that
    #### we have one to modify and one to keep constant
    for row in range(len(origFrags)):
        frags.append(origFrags[row])
        tempFrags.append(origFrags[row])

    ## Builds the guide tree for how to align frags
    pairwiseMatrix = generateMatrix(frags)
    while len(pairwiseMatrix) > 2:
        val = findLargestVal(pairwiseMatrix)
        guideTree.append(val)
        alignedSeqs = align(frags[val[1]], frags[val[2]])
        joined = merge(alignedSeqs[0], alignedSeqs[1])
        frags[val[1]] = joined
        del frags[val[2]]
        pairwiseMatrix = condense(pairwiseMatrix, val[1], val[2])

    ## Just a formality to grab the last to to be clustered and put onto the tree
    lastTwo = findLargestVal(pairwiseMatrix)
    alignedSeqs = align(frags[lastTwo[1]], frags[lastTwo[2]])
    finalJoined = merge(alignedSeqs[0], alignedSeqs[1])
    guideTree.append(lastTwo)

    #### Builds the consensus sequence following the guide tree and pulls out the
    #### consensus seq and the percent list
    consensusSeqData = buildConsensusSeq(guideTree, tempFrags)
    consensusSeq = consensusSeqData[0]
    percentList = consensusSeqData[1]

    
    ## Aligns all of the seqs according to the consensus sequence
    for i in range(len(origFrags)):
        alignedSeq = align(consensusSeq, origFrags[i])
        finalSeqs.append(alignedSeq[1])
    
    return [finalSeqs,percentList,consensusSeq]

'''
Assuming four fragments are in a single file then pull them apart
into two different groups. Both of which will undergo haplotype assembly
'''
def pullApart(frags):

    samePairs = []
    #### Much like building a distance matrix align all fragments
    #### to one another and find the ones that have the greatest
    #### overlap
    for row in range(len(frags)):
        for col in range(len(frags)):
            if col > row:
                frag1 = ''
                frag2 = ''
                alignment = align(frags[row],frags[col])
                overlapCounter = 0
                for letter in range(len(alignment[0])):
                    if alignment[0][letter] != '-' and alignment[1][letter] != '-':
                        overlapCounter += 1
                if overlapCounter > 13:
                    samePairs.append([row, col])
    #### Pull out the first two pairs that go togther and put them
    #### onto one "side" from there check to see if the next element
    #### hase a pair already in the side list. Afterwards the remaining
    #### positions not is side1 are in side 2
    side1 = []
    side2 = []
    side1.append(samePairs[0][0])
    side1.append(samePairs[0][1])
    tempPairs = []
    for i in range(len(samePairs)):
        val1 = samePairs[i][0]
        val2 = samePairs[i][1]
        for j in range(len(side1)):
            if val1 == side1[j]:
                if val2 not in side1:
                    side1.append(val2)
            if val2 == side1[j]:
                if val1 not in side1:
                    side1.append(val1)
    fullList = []
    for i in range(len(frags)):
        fullList.append(i)
        
    for j in range(len(fullList)):
        inSide1 = False
        for k in range(len(side1)):
            if fullList[j] == side1[k]:
                inSide1 = True
        if not inSide1:
            side2.append(fullList[j])

    return [side1, side2]

def complement(seq1, seq2, seq3, seq4):
    compDict = {'A' : 'T', 'G': 'C', 'T' : 'A', 'C' : 'G'}

    #### First two string are reversed 3' ---- 5' need to
    #### complement seqs 3 and 4 and compare
    compSeq1 = ''
    for i in range(len(seq3)):
        compSeq1 += compDict[seq3[i]]

    compSeq2 = ''
    for j in range(len(seq4)):
        compSeq2 += compDict[seq4[i]]

    #### Compare the now complemented strings to
    #### one of the reverse complement strings.
    #### if the errors are low then they match
    #### otherwise it is paired with the other
    numMismatch = 0
    for letter in range(len(compSeq1)):
        if compSeq1[letter] != seq1[letter]:
            numMismatch += 1

    if numMismatch < 3:
        print("Haplotype 1 and complement")
        print(seq3, "Seq 3")
        print(seq1, "Seq 1")
        print("")
        print("Haplotype 2 and complement")
        print(seq4, "Seq 4")
        print(seq2, "Seq 2")
    else:
        print("Haplotype 1 and complement")
        print(seq3, "Seq 3")
        print(seq2, "Seq 2")
        print("")
        print("Haplotype 2 and complement")
        print(seq4, "Seq 4")
        print(seq1, "Seq 1")
        
'''
Main body of the code that runs the rest of it
'''
def main():
    file = input("Please input the file name you wish to work with: ")
    allFrags = readfasta(file)
    origFrags = []

    #### Pull all of the frags out and put them into a new list
    stringFrags = []
    for row in range(len(allFrags)):
        stringFrags.append(allFrags[row][2])

    #### Conslidate the frags list to reduce the number of fragments
    consolidateFrags = consolidate(stringFrags)


    #### Build lists containing only the consolidated frags list
    tempFrags = []
    for row in range(len(consolidateFrags)):
        tempFrag = []
        if consolidateFrags[row]:
            for letter in range(len(consolidateFrags[row])):
                tempFrag.append(consolidateFrags[row][letter])
            origFrags.append(tempFrag)

    twoHalves = pullApart(origFrags)
    twoHalves[0].sort()
    half1Pos = twoHalves[0]
    half1 = []

    ######### Using one half of the sets find the two haplotypes
    for row in range(len(half1Pos)):
        half1.append(origFrags[half1Pos[row]])
    finalFragData = alignSeqs(half1)
    haplotypes = haplotypeAssembly(finalFragData)

    print("First Half")
    ######### Hap 1 ########
    hap1 = []
    for pos in range(len(haplotypes[1])):
        hap1.append(half1[haplotypes[1][pos]])
    finalFragData1 = alignSeqs(hap1)
    half1Hap1 = ''
    for letter in range(len(finalFragData1[2])):
        half1Hap1 += finalFragData1[2][letter]
    print(half1Hap1, "Seq 1")

    ######### Hap 2 #########
    hap2 = []
    for pos in range(len(haplotypes[2])):
        hap2.append(half1[haplotypes[2][pos]])
    finalFragData2 = alignSeqs(hap2)
    half1Hap2 = ''
    for letter in range(len(finalFragData2[2])):
        half1Hap2 += finalFragData2[2][letter]
    print(half1Hap2, "Seq 2")

    ######## Now for the other half
    twoHalves[1].sort()
    half2Pos = twoHalves[1]
    half2 = []
    
    for row in range(len(half2Pos)):
        half2.append(origFrags[half2Pos[row]])
    finalFragData = alignSeqs(half2)
    otherHaplotypes = haplotypeAssembly(finalFragData)

    print("")
    print("Second Half")
    ######### Hap 1 ########
    otherHap1 = []
    for pos in range(len(otherHaplotypes[1])):
        otherHap1.append(half2[otherHaplotypes[1][pos]])
    finalFragData3 = alignSeqs(otherHap1)
    half2Hap1 = ''
    for letter in range(len(finalFragData3[2])):
        half2Hap1 += finalFragData3[2][letter]
    print(half2Hap1, "Seq 3")

    ######### Hap 2 #########
    otherHap2 = []
    for pos in range(len(otherHaplotypes[2])):
        otherHap2.append(half2[otherHaplotypes[2][pos]])
    finalFragData4 = alignSeqs(otherHap2)
    half2Hap2 = ''
    for letter in range(len(finalFragData4[2])):
        half2Hap2 += finalFragData4[2][letter]
    print(half2Hap2, "Seq 4")
        
    reverseHalf1Hap1 = half1Hap1[::-1]
    reverseHalf1Hap2 = half1Hap2[::-1]

    print("")
    complement(reverseHalf1Hap1, reverseHalf1Hap2, half2Hap1, half2Hap2)

    
    
main()
