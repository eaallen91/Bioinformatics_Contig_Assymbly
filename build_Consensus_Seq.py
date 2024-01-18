'''''
Test program to build consensus seqs

Eric Allen
'''''
from Final_Semiglobal import align

'''
Take two aligned seqs and merge them into a single
consensus seq
'''
def merge(frag1,frag2):
    consensusSeq = []
    
    for pos in range(len(frag1)):
        if frag1[pos] == '-':
            consensusSeq.append(frag2[pos])
        elif frag2[pos] == '-':
            consensusSeq.append(frag1[pos])
        elif frag1[pos] == frag2[pos]:
            consensusSeq.append(frag1[pos])
        else:
            consensusSeq.append(frag1[pos]+frag2[pos])

    return consensusSeq

'''
Returns the string containing the letters that appear
the most often at any given spot
'''
def findMax(A,C,T,G):

    maxLetter = ''

    maxVal = max(A,C,T,G)
    
    if maxVal == A:
        maxLetter += 'A'
    elif maxVal == C:
        maxLetter += 'C'
    elif maxVal == T:
        maxLetter += 'T'
    elif maxVal == G:
        maxLetter += 'G'

    return maxLetter

'''
Goes through the consensus seq and find positions that
are longer than 1 character and condenses it by making
that position equal to the letter that appears most often
'''
def condense(startConsensus):
    
    consensusLetterPercent = [ [0 for i in range(0,3)] for j in range(len(startConsensus)) ]
    for pos in range(len(startConsensus)):
        freqA = 0
        freqC = 0
        freqT = 0
        freqG = 0
        for letter in startConsensus[pos]:
            if letter == 'A':
                freqA += 1
            elif letter == 'C':
                freqC += 1
            elif letter == 'T':
                freqT += 1
            elif letter == 'G':
                freqG += 1
        total = len(startConsensus[pos])
        consensusLetterPercent[pos] = [freqA/total, freqC/total, freqT/total, freqG/total]
        mostFreq = findMax(freqA,freqC,freqT,freqG)
        startConsensus[pos] = mostFreq
        
    results = [startConsensus, consensusLetterPercent]
        
    return results
                    
            
'''
Build the consensusSeq based upon the guide tree
'''
def buildConsensusSeq(guideTree, frags):

    for row in range(len(guideTree)):
        frag1 = frags[guideTree[row][1]]
        frag2 = frags[guideTree[row][2]]

        alignedFrags = align(frag1, frag2)
        consensus = merge(alignedFrags[0],alignedFrags[1])
        
        frags[guideTree[row][1]] = consensus
        del frags[guideTree[row][2]]

    finalConsensusData = condense(frags[0])
    return finalConsensusData

        

        
