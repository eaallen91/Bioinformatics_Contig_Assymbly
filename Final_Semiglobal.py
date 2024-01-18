'''
Eric Allen

Optimal alignment of two sequences
'''
from readfasta import readfasta

def align(dnaSeq1, dnaSeq2):
    gap = -100000
    terminalGap = 0
    mismatch = -11
    match = 5
    
    dnaMatrix1 = []
    dnaMatrix2 = []

    for i in dnaSeq1:
        dnaMatrix1.append(i)

    for j in dnaSeq2:
        dnaMatrix2.append(j)

    ## Add the empty spot onto the beginning of the strings. So that they may
    ## match up with the matrix
    dnaMatrix1.insert(0, "-")
    dnaMatrix2.insert(0, "-")

    scoringMatrix = [ [0 for i in range(len(dnaMatrix1)) ] for j in range(len(dnaMatrix2))]

    scoringMatrix[0][0] = 0

    ## Fill the empty string col
    for col in range(1, len(dnaMatrix1)):
        scoringMatrix[0][col] = 0

    ## Fill the empty string row
    for row in range(1, len(dnaMatrix2)):
        scoringMatrix[row][0] = 0

    ## Fill the middle of the matrix
    for col in range(1, len(dnaMatrix1)):
        for row in range(1, len(dnaMatrix2)):
            noGap = 0
            ### Checks to see if there are multiple values it can take
            if len(dnaMatrix1[col]) != 1:
                noGap = scoringMatrix[row-1][col-1]
                for letter in dnaMatrix2[row]:
                    if letter in dnaMatrix1[col]:
                        noGap += match
                    else:
                        noGap += mismatch
                    
            elif len(dnaMatrix2[row]) != 1:
                noGap = scoringMatrix[row-1][col-1]
                for letter in dnaMatrix1[col]:
                    if letter in dnaMatrix2[row]:
                        noGap += match
                    else:
                        noGap += mismatch
            ### Common mismatch
            if dnaMatrix1[col] != dnaMatrix2[row]:
                noGap = scoringMatrix[row-1][col-1] + mismatch
            ### Common Match
            elif dnaMatrix1[col] == dnaMatrix2[row]:
                noGap = scoringMatrix[row-1][col-1] + match
            if col == len(scoringMatrix[0]):
                scoringMatrix[row][col] = max(noGap,
                                              scoringMatrix[row][col-1] + terminalGap,
                                              scoringMatrix[row-1][col] + gap)
            elif row == len(scoringMatrix):
                scoringMatrix[row][col] = max(noGap,
                                              scoringMatrix[row][col-1] + gap,
                                              scoringMatrix[row-1][col] + terminalGap)
            else:
                scoringMatrix[row][col] = max(noGap,
                                              scoringMatrix[row][col-1] + gap,
                                              scoringMatrix[row-1][col] + gap)

    ## Find  the largest element in the last row or col. That will be the
    ## starting point
    largestVal = 0
    startingCol = 0
    startingRow = 0

    ## Checks the last col for the largest element
    for row in range(len(dnaMatrix2)):
        if scoringMatrix[row][len(dnaMatrix1)-1] > largestVal:
            largestVal = scoringMatrix[row][len(dnaMatrix1)-1]
            startingRow = row
            startingCol = len(dnaMatrix1)-1
            
    ## Checks the last row for the largest element
    for col in range(len(dnaMatrix1)):
        if scoringMatrix[len(dnaMatrix2)-1][col] > largestVal:
            largestVal = scoringMatrix[len(dnaMatrix2)-1][col]
            startingRow = len(dnaMatrix2)-1
            startingCol = col

            
    col = startingCol
    row = startingRow
    finalVal = scoringMatrix[row][col]

    lastCol = len(dnaMatrix1)-1
    lastRow = len(dnaMatrix2)-1

    ## Do backtracing starting at the bottom left corner.
    resultSeq1 = []
    resultSeq2 = []
    
    ## Inserts gaps at the end of sequence 1 or 2 depending on where the starting
    ## point is
    while((lastRow != row) or (lastCol != col)):
        if lastCol > col:
            resultSeq1.insert(0, dnaMatrix1[lastCol])
            resultSeq2.insert(0, "-")
            lastCol -= 1
        elif lastRow > row:
            resultSeq1.insert(0, "-")
            resultSeq2.insert(0, dnaMatrix2[lastRow])
            lastRow -= 1

    done = False
    while(done == False):
        if (row == 0) and (col == 0):
            break
        ## Gap at the beginning of seq1
        if col == 0:
            resultSeq1.insert(0, "-")
            resultSeq2.insert(0, dnaMatrix2[row])
            row -= 1
        ## Gap at the beginning of seq2
        elif row == 0:
            resultSeq1.insert(0, dnaMatrix1[col])
            resultSeq2.insert(0, "-")
            col -= 1       
        ## Makes it so a gap can only be introduced at the end of seq1
        elif col == lastCol:
            if scoringMatrix[row][col] == (scoringMatrix[row-1][col] + gap):
                resultSeq1.insert(0, "-")
                resultSeq2.insert(0, dnaMatrix2[row])
                row -= 1
            else:
                resultSeq1.insert(0, dnaMatrix1[col])
                resultSeq2.insert(0, dnaMatrix2[row])
                col -= 1
                row -= 1
        ## Makes it so a gap can only be introduces at the end of seq2
        elif row == lastRow:
            if scoringMatrix[row][col] == (scoringMatrix[row][col-1] + gap):
                resultSeq1.insert(0, dnaMatrix1[col])
                resultSeq2.insert(0, "-")
                col -= 1
            else:
                resultSeq1.insert(0, dnaMatrix1[col])
                resultSeq2.insert(0, dnaMatrix2[row])
                col -= 1
                row -= 1
        ## Otherwise it is a diagnoal movement 
        else:
            resultSeq1.insert(0, dnaMatrix1[col])
            resultSeq2.insert(0, dnaMatrix2[row])
            col -= 1
            row -= 1 
    finalSeqs = [resultSeq1, resultSeq2, finalVal]
    return finalSeqs
