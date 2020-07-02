import sys
import numpy as np

def readFile(fileName, linesToSkipInitial,linesToSkipLater, linesToRead,timeStepsToSkip):
    # function takes the file name/location if not in same workspace
    # the number of lines to skip of useless garbage
    # the number of lines we want to get data from
    # a list of col indices that you want to read
    # so if you want to read the 3rd and 4th column you
    # would pass a list [2, 3]
    # 
    
    # open file and read lines
    fileObj = open(fileName, 'r')
    allLinesList = fileObj.readlines()

    # create list for output data to be stored in
    outputRdfList = [0]*linesToRead
    rList = [0]*linesToRead

    for removeIndx in range(linesToSkipInitial+timeStepsToSkip*(linesToSkipLater+linesToRead)):
        del allLinesList[0]
    
    currTimeStep = 0
    for indx in range(len(allLinesList)):
        outIndx = indx%(linesToRead + linesToSkipLater) - 1
        if(outIndx != -1):
            currLineList = allLinesList[indx].split()
            rList[outIndx] = float(currLineList[1])
            tempGr = (float(currLineList[2]) + outputRdfList[outIndx]*currTimeStep)/(currTimeStep + 1)
            outputRdfList[outIndx] = tempGr
            
        else:
            currTimeStep += 1

    
    return rList, outputRdfList


def readMyCOMFile(fileName, linesToRead):


    linesToSkipInitial = 3
    linesToSkipLater = 1
    # open file and read lines
    fileObj = open(fileName, 'r')
    allLinesList = fileObj.readlines()
    
    # create list for output data to be stored in
    
    outputList = []
    tsList = []

    for removeIndx in range(linesToSkipInitial):
        del allLinesList[0]


    for indx in range(len(allLinesList)):
        outIndx = indx%(linesToRead + linesToSkipLater) 
        if(outIndx != 0):
            currLineList = allLinesList[indx].split()
            tsList.append([float(currLineList[1]),float(currLineList[2]),float(currLineList[3])])
        elif(indx>0):
            copyList = tsList.copy()
            outputList.append(copyList)
            tsList = []
    copyList = tsList.copy()
    outputList.append(copyList)        
    
    return outputList

def readMyRdfFileOneAtaTime(fileName, linesToSkipInitial,linesToSkipLater, linesToRead):
    # function takes the file name/location if not in same workspace
    # the number of lines to skip of useless garbage
    # the number of lines we want to get data from
    # a list of col indices that you want to read
    # so if you want to read the 3rd and 4th column you
    # would pass a list [2, 3]
    
    # open file and read lines
    fileObj = open(fileName, 'r')
    allLinesList = fileObj.readlines()
    fileObj.close()

    # create list for output data to be stored in
    outputRdfList = []
    tsList = []
    rList = []

    for removeIndx in range(linesToSkipInitial):
        del allLinesList[0]
    
    currTimeStep = 0
    for lineNum in range(len(allLinesList)):
        outIndx = lineNum%(linesToRead + linesToSkipLater) - 1
        if(outIndx != -1 and currTimeStep == 0):
            currLineList = allLinesList[lineNum].split()
            rList.append(float(currLineList[1]))
            tempGr = float(currLineList[2])
            tsList.append(tempGr)
        elif(outIndx != -1):
            currLineList = allLinesList[lineNum].split() 
            tempGr = float(currLineList[2])
            tsList.append(tempGr)
        elif(lineNum>0):
            addGr = tsList.copy()
            outputRdfList.append(addGr)
            tsList = []
            currTimeStep += 1

    
    return rList, outputRdfList

def readMyPeOUTFile(fileName):
    
    data = np.loadtxt(fileName)
    
    TS = data[:,0]
    PE = data[:,1]
    
    return TS, PE

def readMyDumpFileForXYZ(fileName, numAtoms):
    # function takes the dump file name/location
    # and iterates through to get the coordinates of each atom and puts it in
    # a list with the following format for n total atoms and ts timesteps
    # [[[x,y,z],[x,y,z],...n],[[x,y,z],[x,y,z],...n],...ts]
    # this function assumes the dump command in lammps was x y z ix iy iz

    linesToAlwaysSkip = 9
    # open file and read lines
    fileObj = open(fileName, 'r')
    sys.stdout.write('READING DUMPFILE\n')
    sys.stdout.flush()
    allLinesList = fileObj.readlines()
    sys.stdout.write('SORTING DUMPFILE\n')
    sys.stdout.flush()
    fileObj.close()
        
    # create list for output data to be stored in
    outputCoordsList = [] #[[[x,y,z],[x,y,z]...n],[[x,y,z],[x,y,z]...n],...ts]
    timeStepCoords = [] # [[x,y,z]*nAtoms]
    for indx in range(len(allLinesList)):
        sys.stdout.write('CURRENT LINE: ' + str(indx) + '\n')
        sys.stdout.flush()
        outIndx = indx%(numAtoms+linesToAlwaysSkip)
        if(outIndx not in range(linesToAlwaysSkip)):
            currLineList = allLinesList[indx].split()
            currCoords = [float(currLineList[-6]),float(currLineList[-5]),float(currLineList[-4])]
            timeStepCoords.append(currCoords)
        if(outIndx == numAtoms+linesToAlwaysSkip-1):
            tmpList = timeStepCoords.copy()
            outputCoordsList.append(tmpList)
            timeStepCoords = []
    sys.stdout.write('FINISHED READING DUMPFILE\n')
    sys.stdout.flush()
    return outputCoordsList

def formatDumpFileForXYZ(inFile,outFile, numAtoms):
    # function takes the dump file name/location
    # and iterates through to get the coordinates of each atom and puts it in
    # a list with the following format for n total atoms and ts timesteps
    # [[[x,y,z],[x,y,z],...n],[[x,y,z],[x,y,z],...n],...ts]
    # this function assumes the dump command in lammps was x y z ix iy iz

    linesToAlwaysSkip = 9
    # open file and read lines
    inFileObj = open(inFile, 'r')
    sys.stdout.write('READING DUMPFILE\n')
    sys.stdout.flush()
    allLinesList = inFileObj.readlines()
    sys.stdout.write('SORTING DUMPFILE\n')
    sys.stdout.flush()
    inFileObj.close()
    sys.stdout.write('WRITING TO OUTPUT FILE\n')
    sys.stdout.flush()
    outFileObj = open(outFile,'w')
    
    # create list for output data to be stored in
    outputCoordsList = [] #[[[x,y,z],[x,y,z]...n],[[x,y,z],[x,y,z]...n],...ts]
    timeStepCoords = [] # [[x,y,z]*nAtoms]
    for indx in range(len(allLinesList)):
        sys.stdout.write('CURRENT LINE: ' + str(indx) + '\n')
        sys.stdout.flush()
        outIndx = indx%(numAtoms+linesToAlwaysSkip)
        if(outIndx not in range(linesToAlwaysSkip)):
            currLineList = allLinesList[indx].split()
            outFileObj.write(currLineList[0] + '\t' + "{:.6f}".format(float(currLineList[-6])) + '\t' + "{:.6f}".format(float(currLineList[-5])) + '\t' + "{:.6f}".format(float(currLineList[-4])) + '\n')
        elif(outIndx == 1):
            outFileObj.write('# ' + allLinesList[indx])
        
    sys.stdout.write('FINISHED READING DUMPFILE\n')
    sys.stdout.flush()
    outFileObj.close()

def readMyDumpFileForAvgVal(fileName,linesToAlwaysSkip, numVals, indxWanted):
    # function takes the dump file name/location
    # and iterates through to get a value held at indxWanted of each line
    # for numVals in each time step. linesToAlwaySkip = 9 
    
    # open file and read lines
    fileObj = open(fileName, 'r')
    allLinesList = fileObj.readlines()
    fileObj.close()
    # create list for output data to be stored in
    outputValsList = [] #[[val,val,val...n],[val,val,val...n]...ts]
    timeStepVals = [] #[val,val,val...n]
    for indx in range(len(allLinesList)):
        outIndx = indx%(numVals+linesToAlwaysSkip)
        if(outIndx not in range(linesToAlwaysSkip)):
            currLineList = allLinesList[indx].split()
            currVal = currLineList[indxWanted]
            timeStepVals.append(float(currVal))
        if(outIndx == numVals+linesToAlwaysSkip-1):
            tmplList = timeStepVals.copy()
            outputValsList.append(tmpList)     
            timeStepVals = []
    
    return outputValsList

def readMyVMDradGyrData(fileName):
    fileObj = open(fileName,'r')
    allLinesList = fileObj.readlines()
    radGyr = [float(line) for line in allLinesList]
    fileObj.close()
    return radGyr

def readFileForMultipleLines(fileName,cols):
    # takes a string filename for the file of interest
    # and a list of columns that you want (will return one long list)
    # should be indexed how you want if you want first column, pass 0
    # loads the data from these columns and returns them in a list of lists
    # all # lines will be ignored 

    data = np.loadtxt(fileName)
    outputList = []
    for colNum in cols:
        currData = data[:,colNum]
        copyList = currData.copy()
        outputList.append(copyList)

    return outputList

def readFormattedDump(fileName,atomsPerMol,mols):
    # takes a string filename for the file of interest
    # and a list of columns that you want (will return one long list)
    # should be indexed how you want if you want first column, pass 0
    # loads the data from these columns and returns them in a list of lists
    # all # lines will be ignored 
    
    data = np.loadtxt(fileName)
    timestepList = []
    outputList = []
    totalAtoms = atomsPerMol*mols
    first = 0
    last = atomsPerMol
    currTS = 1
    while(last<= int(len(data[:,0]))):
        tmpList = data[first:last,1:4]
        copyList = tmpList.copy()
        timestepList.append(copyList)
        if(last%totalAtoms == 0):
            sys.stdout.write('CURRENT TS: ' + str(currTS) + '\n')
            sys.stdout.flush()
            currTS += 1
            polsAtTS = timestepList.copy()
            outputList.append(polsAtTS)
            timestepList = []
        first = last
        last += atomsPerMol

    return outputList

def readIndividualDump(fileName,atomsPerMol,totMols):
# function takes the dump file name/location
    # and iterates through to get the coordinates of each atom and puts it in
    # a list with the following format for n total atoms and ts timesteps
    # [[[x,y,z],[x,y,z],...n],[[x,y,z],[x,y,z],...n],...ts]
    # this function assumes the dump command in lammps was x y z ix iy iz

    linesToSkip = 9
    # open file and read lines
    fileObj = open(fileName, 'r')
    allLinesList = fileObj.readlines()
    fileObj.close()
    
    del allLinesList[0:linesToSkip]
    
    # create list for output data to be stored in
    outputCoordsList = [] #[[[x,y,z],[x,y,z]...n],[[x,y,z],[x,y,z]...n]]
    # list for atoms in a polymer
    polymerList = []
    currAtom = 0
    currMol = 1
    while(currMol<=totMols):
        while(currAtom<atomsPerMol):
            currLineList = allLinesList[currAtom+(currMol-1)*atomsPerMol].split()
            currCoords = [float(currLineList[-6]),float(currLineList[-5]),float(currLineList[-4])]
            polymerList.append(currCoords)
            currAtom += 1
        copyList = polymerList.copy()
        outputCoordsList.append(copyList)
        polymerList = []
        currAtom = 0
        currMol += 1
    return outputCoordsList
    