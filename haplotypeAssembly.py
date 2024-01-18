def buildConflictMatrix(finalFrags, percentList):
    ##################### Build the conflict matrix with SNPs ##############
    conflictMatrix = [ ['0' for i in range(len(finalFrags[0]))] for j in range(len(finalFrags)) ]
    freqLettersCol = []
    snpLetters = []

    #### Builds the fragment X snp conflict matrix
    snpCounter = 0
    snpCol = []
    for col in range(len(finalFrags[0])):
        SNP1 = ''
        SNP2 = ''
        for row in range(len(finalFrags)):
            #### If there is a col with conflicting letters then it is a snp
            #### the first letter met is A and the second is B
            if finalFrags[row][col] != '-':
                if not SNP1:
                    SNP1 = finalFrags[row][col]
                elif not SNP2 and finalFrags[row][col] != SNP1:
                    SNP2 = finalFrags[row][col]

        
        if SNP1 and SNP2:
            snpCol.append(snpCounter)
            snpCounter += 1
            #### With A and B assigned go back through the col and check to see the
            #### places that equal A, equal B, neither equals, or a third value (X)
            for row in range(len(finalFrags)):
                if finalFrags[row][col] == SNP1:
                    conflictMatrix[row][col] = 'A'
                elif finalFrags[row][col] == SNP2:
                    conflictMatrix[row][col] = 'B'
                elif finalFrags[row][col] == '-':
                    conflictMatrix[row][col] = '0'
                else:
                    conflictMatrix[row][col] = 'X'
        else:
            snpCol.append(-1)

    return conflictMatrix
    ############################## Conflict Matrix built #######################

def buildFragConflict(finalFrags, conflictMatrix):
    ############### Build fragment conflict graph ####################
    #### Builds the fragment conflict graph by finding a snp in one row
    #### and by finding another row with a conflicting letter in the same
    #### col.
    conflictGraph = [ [] for i in range(len(finalFrags)) ]
    for row in range(len(conflictMatrix)):
        for col in range(len(conflictMatrix[row])):
            if conflictMatrix[row][col] != '0':
                SNP = conflictMatrix[row][col]
                for innerRow in range(len(conflictMatrix)):
                    if conflictMatrix[innerRow][col] != '0':
                        if conflictMatrix[innerRow][col] != SNP:
                            if innerRow not in conflictGraph[row] and innerRow != row:
                                conflictGraph[row].append(innerRow)
    return conflictGraph

    ############################# Fragment conflict graph built ################

def buildSNPConflict(conflictMatrix):
    ######################## Build SNP conflict graph ##########################
    #### Builds the SNP conflict graph.
    snpGraph = [ [] for i in range(0,snpCounter) ]
    for col in range(len(conflictMatrix[0])):
        SNP1 = ''
        frag1 = -1
        SNP2 = ''
        frag2 = -1
        #### Finds a col where a SNP conflicts
        for row in range(len(conflictMatrix)):
            if conflictMatrix[row][col] != '0':
                if not SNP1:
                    SNP1 = conflictMatrix[row][col]
                    frag1 = row
                elif not SNP2:
                    SNP2 = conflictMatrix[row][col]
                    frag2 = row
            #### Once the conflicting snps are obtained loop through
            #### both frags and  try to find another conflicting SNP location
            if frag1 != -1 and frag2 != -1: 
                for innerCol in range(len(conflictMatrix[frag1])):
                    if innerCol != col:
                        snpCorners = 0
                        if SNP1 == SNP2:
                            snpCorners = 2
                        else:
                            snpCorners = 1
                        #### Once it finds another SNP with conflicting vals check to see if 3 corners
                        #### are the same value
                        if conflictMatrix[frag1][innerCol] != '0' and conflictMatrix[frag2][innerCol] != '0':
                            SNP3 = conflictMatrix[frag1][innerCol]
                            SNP4 = conflictMatrix[frag2][innerCol]
                            if SNP1 == SNP3:
                                snpCorners += 1
                            if SNP1 == SNP4:
                                snpCorners += 1

                            if snpCorners == 1 or snpCorners == 3:
                                if snpCol[innerCol] not in snpGraph[snpCol[col]]:
                                    snpGraph[snpCol[col]].append(snpCol[innerCol])
                ########## End of inner for loop                    
                SNP1 = ''
                SNP2 = ''
                frag1 = -1
                frag2 = -1              


    return snpGraph
    ################# End of buildSNPGraph

def depthFirstSearch(snpGraph):
    fullGraph = []
    for row in range(len(snpGraph)):
        fullGraph.append(row)
    ########################### Depth-first search for independent set #########
    independentSet = []
    while fullGraph:
        fewestNeighbors = len(snpGraph[fullGraph[0]])
        fewestPos = fullGraph[0]
        
        #### Find which snp has the fewest neighbors
        for i in range(len(fullGraph)):
            if len(snpGraph[fullGraph[i]]) < fewestNeighbors:
                fewestPos = fullGraph[i]
                fewestNeighbors = len(snpGraph[fullGraph[i]])

        #### Put that snp into the independent set list
        setToRemove = []
        independentSet.append(fewestPos)
        setToRemove.append(fewestPos)
        for j in range(len(snpGraph[fewestPos])):
            setToRemove.append(snpGraph[fewestPos][j])
        
        #### Remove nodes from the graph
        setToRemove.sort()
        setToRemove.reverse()
        for i in range(len(setToRemove)):
            for k in range(len(fullGraph)-1,-1,-1):
                if setToRemove[i] == fullGraph[k]:
                    del fullGraph[k]

    #### Sort and reverse the independent set
    independentSet.sort()
    independentSet.reverse()

    #### Build a list with all the snps in it
    nodesToRemove = []
    for row in range(len(snpGraph)):
        nodesToRemove.append(row)

    #### Remove the independent set from all snps
    for snp in range(len(independentSet)):
        del nodesToRemove[independentSet[snp]]

    #### Obtain the col number for each snp
    nodesToRemove.sort()
    nodePos = []
    for col in range(len(nodesToRemove)):
        for item in range(len(snpCol)):
            if snpCol[item] == nodesToRemove[col]:
                nodePos.append(item)

    #### Sort and reverse the pos and delete from the conflict graph
    nodePos.sort()
    nodePos.reverse()
    for i in range(len(nodePos)):
        for row in range(len(conflictMatrix)):
            del conflictMatrix[row][nodePos[i]]

        
    #### Reconstruct the conflict matrix with only the independent set snps
    finalConflict = [ [] for i in range(len(finalFrags)) ]
    for row in range(len(conflictMatrix)):
        for col in range(len(conflictMatrix[row])):
            if conflictMatrix[row][col] != '0':
                SNP = conflictMatrix[row][col]
                for innerRow in range(len(conflictMatrix)):
                    if conflictMatrix[innerRow][col] != '0':
                        if conflictMatrix[innerRow][col] != SNP:
                            if innerRow not in conflictGraph[row] and innerRow != row:
                                finalConflict[row].append(innerRow)

    return finalConflict

def breadthFirstSearch(conflictGraph):
    ########################## Breadth First search to determine bipartitness ##
    #### Building a breadth-first search tree to determine bipartitness
    fullGraph = []
    for row in range(len(conflictGraph)):
        fullGraph.append(row)
        
    tree = []
    visited = []
    level0 = [0]
    tree.append(level0)
    currentLevel = 0
    #### Starting at node 0 put all neighbors in the next level. Check to
    #### see if it has been visited. And pull of its neighbors for the next level
    while tree[currentLevel]:
        nextLevel = []
        
        for vertex in range(len(tree[currentLevel])):
            if tree[currentLevel][vertex] not in visited:
                visited.append(tree[currentLevel][vertex])
            for neighbor in range(len(conflictGraph[tree[currentLevel][vertex]])):
                nextLevel.append(conflictGraph[tree[currentLevel][vertex]][neighbor])
        tree.append(nextLevel)

        if len(visited) == len(fullGraph):
            break
                
        currentLevel += 1

    #### All nodes in the even rows are in the blue list
    #### All nodes in the odd rows are in the red list
    #### check to see if there are any conflictions, print the final strings
    blueList = []
    redList = []
    bipartite = True
    for row in range(len(tree)):
        for pos in range(len(tree[row])):
            if row%2 == 0:
                if tree[row][pos] not in blueList:
                    blueList.append(tree[row][pos])
            else:
                if tree[row][pos] not in redList:
                    redList.append(tree[row][pos])
                    
    #### Checks for any conflicts in the two lists
    for val in range(len(blueList)):
        if blueList[val] in redList:
            print("CONFLICT NOT BIPARTITE:",blueList[val])
            bipartite = False
            
    return [bipartite, blueList, redList]
    ############################# End of breadth-first search

def haplotypeAssembly(finalFragData):
    finalFrags = finalFragData[0]
    percentList = finalFragData[1]

    conflictMatrix = buildConflictMatrix(finalFrags, percentList)
    conflictGraph = buildFragConflict(finalFrags, conflictMatrix)
    bipartite = breadthFirstSearch(conflictGraph)
    if not bipartite[0]:
        snpGraph = buildSNPConflict(conflictMatrix)
        finalConflict = depthFirstSearch(snpGraph)
        bipartite = breadthFirstSearch(finalConflict)

    return bipartite
    
        
        
