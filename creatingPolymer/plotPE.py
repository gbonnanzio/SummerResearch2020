#! Program Files (x86)/PythonDownload
import matplotlib.pyplot as plt
import statistics as stat 
import math
from readFiles import readMyPeOUTFile


def checkInRange(min1,max1,min2,max2):
    check1 = max2 > min1 and max2 < max1
    check2 = min2 < max1 and min2 > min1
    check3 = max1 < max2 and min1 > min2
    finalCheck = check1 or check2 or check3
    return finalCheck 

def checkWholeList(peList,sdList,myMin,myMax,numToCheck):
    lenList = len(peList)
    if(lenList<numToCheck):
        numToCheck = lenList
    for i in range(numToCheck):
        backIndx = -1-i
        tmpMin = peList[backIndx] - sdList[backIndx]
        tmpMax = peList[backIndx] + sdList[backIndx]
        testFlag = checkInRange(tmpMin,tmpMax, myMin,myMax)
        if(testFlag == False):
            return False
    return True

def main():
    fileList = ['runLammps/PE.OUT']
    colors = ['bo']
    labels = ['Calculated Potential Energy']
    numBlocks = 100
    checkPrev = 5
    
    for fileIndx in range(len(fileList)):
        avgTSList = []
        tmpTSList = []
        avgPEList = []
        tmpPEList = []
        stdevPEList = []
        tmpSDList = []
        ts, pe = readMyPeOUTFile(fileList[fileIndx])
        indx = 0
        endIndx = math.ceil(len(ts)/numBlocks)-1
        for tsIndx in range(endIndx):
            first = tsIndx*numBlocks
            last = (tsIndx+1)*numBlocks-1
            tmpTS = stat.mean(ts[first:last])
            tmpPE = stat.mean(pe[first:last])
            tmpPEstdev = stat.stdev(pe[first:last])
            if(indx == 0):
                tmpTSList.append(tmpTS)
                tmpPEList.append(tmpPE)
                tmpSDList.append(tmpPEstdev)
                indx += 1
            else:
                prevMin = tmpPEList[indx-1] - tmpSDList[indx-1]
                prevMax = tmpPEList[indx-1] + tmpSDList[indx-1]
                currMin = tmpPE - tmpPEstdev
                currMax = tmpPE + tmpPEstdev
                inRange = checkWholeList(tmpPEList,tmpSDList,currMin,currMax,checkPrev)
                if(inRange):
                    tmpTSList.append(tmpTS)
                    tmpPEList.append(tmpPE)
                    tmpSDList.append(tmpPEstdev)
                    if(tsIndx == endIndx-1):
                        copyTS = tmpTSList.copy()
                        avgTSList.append(copyTS)
                        copyPE = tmpPEList.copy()
                        avgPEList.append(copyPE)
                        copySD = tmpSDList.copy()
                        stdevPEList.append(copySD)
                    indx += 1
                else:
                    copyTS = tmpTSList.copy()
                    tmpTSList = []
                    tmpTSList.append(tmpTS)
                    avgTSList.append(copyTS)

                    copyPE = tmpPEList.copy()
                    tmpPEList = []
                    tmpPEList.append(tmpPE)
                    avgPEList.append(copyPE)

                    copySD = tmpSDList.copy()
                    tmpSDList = []
                    tmpSDList.append(tmpPEstdev)
                    stdevPEList.append(copySD) 

                    if(tsIndx == endIndx-1):
                        copyTS = tmpTSList.copy()
                        avgTSList.append(copyTS)
                        copyPE = tmpPEList.copy()
                        avgPEList.append(copyPE)
                        copySD = tmpSDList.copy()
                        stdevPEList.append(copySD)
                    indx = 1

        #print(len(avgTSList))  
        plt.figure(fileIndx+1)
        for pltIndx in range(len(avgTSList)):
            if(pltIndx%2 == 0):
                plt.errorbar(avgTSList[pltIndx],avgPEList[pltIndx],yerr = stdevPEList[pltIndx], ecolor = 'b',color = 'k',ms = 3,capsize = 3)
            else:
                plt.errorbar(avgTSList[pltIndx],avgPEList[pltIndx],yerr = stdevPEList[pltIndx], ecolor = 'r',color = 'k',ms = 3,capsize = 3)
        plt.xlabel('Timestep')
        plt.ylabel('Potential Energy (KT)')
        plt.legend(['Avg. Potential Energy'])
        plt.show()



main()	

