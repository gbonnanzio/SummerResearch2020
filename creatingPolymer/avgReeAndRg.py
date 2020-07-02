#! Program Files (x86)/PythonDownload
import matplotlib.pyplot as plt
import statistics as stat 
import math
import numpy as np
from readFiles import (readFileForMultipleLines, readMyPeOUTFile)


def graphReeVsTime(fileName):
    colsOfInterest = [0,1]
    reeData = readFileForMultipleLines(fileName,colsOfInterest)
    avgRee = reeData[0]
    sdRee = reeData[1]
    timestep = np.linspace(0,100*len(avgRee),len(avgRee))
    plt.errorbar(timestep,avgRee,yerr = sdRee)
    plt.show()

def graphRgVsTime(fileName):
    colsOfInterest = [0,1]
    rgData = readFileForMultipleLines(fileName,colsOfInterest)
    avgRg = rgData[0]
    sdRg = rgData[1]
    timestep = np.linspace(0,100*len(avgRg),len(avgRg))
    plt.errorbar(timestep,avgRg,yerr = sdRg)
    plt.show()

    
def checkInRange(min1,max1,min2,max2):
    check1 = max2 > min1 and max2 < max1
    check2 = min2 < max1 and min2 > min1
    check3 = max1 < max2 and min1 > min2
    finalCheck = check1 or check2 or check3
    return finalCheck 

def checkPrevInList(yList,sdList,myMin,myMax,numToCheck):
    # function takes a list of y values and their corresponding errors
    # and checks whether the extremes of the error of an avg (myMin and myMax)
    # overlap ALL of the previous "numToCheck" blocks in the current segment
    # if it overlaps with all the errors of these blocks it returns true
    # else it breaks and returns false
    lenList = len(yList)
    if(lenList<numToCheck):
        numToCheck = lenList
    for i in range(numToCheck):
        backIndx = -1-i
        tmpMin = yList[backIndx] - sdList[backIndx]
        tmpMax = yList[backIndx] + sdList[backIndx]
        testFlag = checkInRange(tmpMin,tmpMax, myMin,myMax)
        if(testFlag == False):
            return False
    return True

def main():

    # variables to change
    reeFile = 'Calc_Ree.txt'
    rgFile = 'Calc_Rg.txt'
    fileList = [reeFile,rgFile]
    xAxisLabels = ['Timestep','Timestep']
    yAxisLabels = ['$R_{ee}$ [d]','$R_g^2$ [$d^2$]']
    legends = ['$R_{ee}$ [d]','$R_g^2$ [$d^2$]']
    # how often coords were taken from lammps
    coordsEvery_TS = 100
    # how many timesteps you want averaged
    numBlocks = 200
    # how many previous averages you want to check to see if they fall in range
    checkPrev = 3
    
    for fileIndx in range(len(fileList)):
        # list of the average x to plot
        avgXList = []
        # holds the segment currently being built
        tmpXList = []
        # list of the average y to plot
        avgYList = []
        # holds the segment currently being built
        tmpYList = []
        # list of the standard deviations to plot as error bars
        stdevYList = []
        # holds the SD of the segment currently being built
        tmpSDList = []
        # get data in the form [[y1,y2...yn],[SD1,SD2...SDn]]
        fullList = readFileForMultipleLines(fileList[fileIndx],[0,1])
        allY = fullList[0]
        allYSD = fullList[1]
        xVals = np.linspace(0,coordsEvery_TS*len(allY),len(allY))
        indx = 0
        endIndx = math.ceil(len(xVals)/numBlocks)-1
        for tsIndx in range(endIndx):
            # determine where to sample from 
            first = tsIndx*numBlocks
            last = (tsIndx+1)*numBlocks
            # find avg of current block
            tmpX = stat.mean([xVals[first],xVals[last]])
            tmpY = stat.mean(allY[first:last])
            tmpSD = allYSD[first:last]
            # account for propagation of error
            sqrSD = [i**2 for i in tmpSD]
            tmpSD = (math.sqrt(sum(sqrSD)))/numBlocks
            # if first block in segment
            if(indx == 0):
                tmpXList.append(tmpX)
                tmpYList.append(tmpY)
                tmpSDList.append(tmpSD)
                # increment, added one block to segment
                indx += 1
            # else we are either adding to an existing segment
            # or we are creating a new segment
            else:
                # pass the current segment list and determine wheter or not
                # the current block belongs in segment or not
                currMin = tmpY - tmpSD
                currMax = tmpY + tmpSD
                inRange = checkPrevInList(tmpYList,tmpSDList,currMin,currMax,checkPrev)
                # if it does overlap all of the previous checkPrev blocks in segment
                if(inRange):
                    # add block to current segment
                    tmpXList.append(tmpX)
                    tmpYList.append(tmpY)
                    tmpSDList.append(tmpSD)
                    # if this is the last group to be added
                    if(tsIndx == endIndx-1):
                        copyX = tmpXList.copy()
                        avgXList.append(copyX)
                        copyY = tmpYList.copy()
                        avgYList.append(copyY)
                        copySD = tmpSDList.copy()
                        stdevYList.append(copySD)
                    indx += 1
                # else add old segment to overall list and create new segment 
                else:
                    # create a copy of prev segment, create a new segment
                    copyX = tmpXList.copy()
                    tmpXList = []
                    tmpXList.append(tmpX)
                    # add old segment to overall list
                    avgXList.append(copyX)
                    
                    # do the same for Y
                    copyY = tmpYList.copy()
                    tmpYList = []
                    tmpYList.append(tmpY)
                    avgYList.append(copyY)

                    copySD = tmpSDList.copy()
                    tmpSDList = []
                    tmpSDList.append(tmpSD)
                    stdevYList.append(copySD) 
                    # if this is the last group to be added
                    if(tsIndx == endIndx-1):
                        copyX = tmpXList.copy()
                        avgXList.append(copyX)
                        copyY = tmpYList.copy()
                        avgYList.append(copyY)
                        copySD = tmpSDList.copy()
                        stdevYList.append(copySD)
                    # reset indx to 1 because we are creating a new segment
                    indx = 1

        #print(len(avgXList))  
        # implement a color change every time a new segment begins 
        plt.figure(fileIndx+1)
        for pltIndx in range(len(avgXList)):
            if(pltIndx%2 == 0):
                plt.errorbar(avgXList[pltIndx],avgYList[pltIndx],yerr = stdevYList[pltIndx], ecolor = 'b',color = 'k',ms = 3,capsize = 3)
            else:
                plt.errorbar(avgXList[pltIndx],avgYList[pltIndx],yerr = stdevYList[pltIndx], ecolor = 'r',color = 'k',ms = 3,capsize = 3)
        plt.xlabel(xAxisLabels[fileIndx])
        plt.ylabel(yAxisLabels[fileIndx])
        plt.legend([legends[fileIndx]])
        plt.show()



main()	

