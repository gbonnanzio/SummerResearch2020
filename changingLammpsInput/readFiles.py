import sys
import numpy as np
def readMyRdfFile(fileName, linesToSkipInitial,linesToSkipLater, linesToRead,timeStepsToSkip):
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

def readMyDumpFileForXYZ(fileName,linesToAlwaysSkip, numAtoms):
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
    outputCoordsList = [] #[[[x,y,z],[x,y,z]...n],[[x,y,z],[x,y,z]...n],...ts]
    timeStepCoords = [] # [[x,y,z]*nAtoms]
    for indx in range(len(allLinesList)):
        outIndx = indx%(numAtoms+linesToAlwaysSkip)
        if(outIndx not in range(linesToAlwaysSkip)):
            currLineList = allLinesList[indx].split()
            currCoords = [float(currLineList[2]),float(currLineList[3]),float(currLineList[4])]
            timeStepCoords.append(currCoords)
        if(outIndx == numAtoms+linesToAlwaysSkip-1):
            tmpList = timeStepCoords.copy()
            outputCoordsList.append(tmpList)
            timeStepCoords = []
    
    return outputCoordsList

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
