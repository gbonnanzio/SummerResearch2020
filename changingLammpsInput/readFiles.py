import sys
import matplotlib.pyplot as plt

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


