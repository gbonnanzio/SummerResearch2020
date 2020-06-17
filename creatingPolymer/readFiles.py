import sys

def readMyRdfFile(fileName, linesToSkipInitial,linesToSkipLater, linesToRead, numTimeSteps):
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

    for removeIndx in range(linesToSkipInitial):
        del allLinesList[0]
    #sys.stdout.write(allLinesList[0])
    #sys.stdout.flush()
    currTimeStep = 0
    for indx in range(len(allLinesList)):
        outIndx = indx%(linesToRead + linesToSkipLater) - 1
        if(outIndx != -1):
            #sys.stdout.write(str(indx))
            #sys.stdout.flush()
            currLineList = allLinesList[indx].split()
            #sys.stdout.write(currLineList[0]+'\n')
            #sys.stdout.flush()
            rList[outIndx] = float(currLineList[1])
            tempGr = (float(currLineList[2]) + outputRdfList[outIndx]*currTimeStep)/(currTimeStep + 1)
            outputRdfList[outIndx] = tempGr
            
        else:
            currTimeStep += 1
            #tmpLine = allLinesList[indx]
            #sys.stdout.write(tmpLine[0])

    
    return rList, outputRdfList

def readMyPeOUTFile(fileName):
    
    data = np.loadtxt(fileName)
    
    TS = data[2:,0]
    PE = data[2:,1]
    
    return TS, PE

def readMyDumpFile(fileName,linesToAlwaysSkip, numAtoms):
    # function takes the dump file name/location
    # and iterates through to get the coordinates of each atom
    # in a list [x,y,z] held in another list containing all the coordinates
    # of every atom in the sim. held in another list of all the timesteps
    # so it is a triple nested list. linesToAlwaySkip = 9 
    
    # open file and read lines
    fileObj = open(fileName, 'r')
    allLinesList = fileObj.readlines()

    # create list for output data to be stored in
    outputCoordsList = []
    timeStepCoords = []
    for indx in range(len(allLinesList)):
        outIndx = indx%(numAtoms+linesToAlwaysSkip)
        if(outIndx not in range(linesToAlwaysSkip)):
            currLineList = allLinesList[indx].split()
            currCoords = [currLineList[2],currLineList[3],currLineList[4]]
            timeStepCoords.append(currCoords)
        elif(outIndx == numAtoms+linesToAlwaysSkip-1):
            outputCoordsList.append(timeStepCoords)
            timeStepCoords = timeStepCoords.copy()
            timeStepCoords = []
    
    return outputCoordsList

def readMyDumpFileForAvgVal(fileName,linesToAlwaysSkip, numVals):
    # function takes the dump file name/location
    # and iterates through to get the coordinates of each atom
    # in a list [x,y,z] held in another list containing all the coordinates
    # of every atom in the sim. held in another list of all the timesteps
    # so it is a triple nested list. linesToAlwaySkip = 9 
    
    # open file and read lines
    fileObj = open(fileName, 'r')
    allLinesList = fileObj.readlines()
    fileObj.close()
    # create list for output data to be stored in
    outputCoordsList = [] #[[val,val,val...n],[val,val,val...n]...ts]
    timeStepCoords = [] #[val,val,val...n]
    for indx in range(len(allLinesList)):
        outIndx = indx%(numVals+linesToAlwaysSkip)
        if(outIndx not in range(linesToAlwaysSkip)):
            currLineList = allLinesList[indx].split()
            currVal = currLineList[0]
            timeStepCoords.append(float(currVal))
    #    elif(outIndx == numVals+linesToAlwaysSkip-1):
     #       tmplList = timeStepCoords.copy()
      ##     timeStepCoords = []
    
    return timeStepCoords
