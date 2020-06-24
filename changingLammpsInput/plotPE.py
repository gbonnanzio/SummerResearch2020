#! Program Files (x86)/PythonDownload
import matplotlib.pyplot as plt
import statistics as stat 
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
    fileList = ['thermo.125.OUT','thermo.25.OUT','thermo.50.OUT']
    colors = ['ro','bs','g^']
    labels = ['ρ = 0.125','ρ = 0.25','ρ = 0.50']
    blockAvg = 25
    avgTSList = []
    tmpTSList = []
    avgPEList = []
    tmpPEList = []
    stdevPEList = []
    tmpSDList = []
    checkLast = 5
    for fileIndx in range(1):
        ts, pe = readMyPeOUTFile(fileList[fileIndx])
        first = 0
        last = blockAvg
        indx = 0
        for tsIndx in range(len(ts)//blockAvg):
            tmpTS = stat.mean(ts[tsIndx*blockAvg:(tsIndx+1)*blockAvg-1])
            tmpPE = stat.mean(pe[tsIndx*blockAvg:(tsIndx+1)*blockAvg-1])
            tmpPEstdev = stat.stdev(pe[tsIndx*blockAvg:(tsIndx+1)*blockAvg-1])
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
                inRange = checkWholeList(tmpPEList,tmpSDList,currMin,currMax,checkLast)
                if(inRange):
                    tmpTSList.append(tmpTS)
                    tmpPEList.append(tmpPE)
                    tmpSDList.append(tmpPEstdev)
                    if(tsIndx == len(ts)//blockAvg-1):
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

                    if(tsIndx == len(ts)//blockAvg-1):
                        copyTS = tmpTSList.copy()
                        avgTSList.append(copyTS)
                        copyPE = tmpPEList.copy()
                        avgPEList.append(copyPE)
                        copySD = tmpSDList.copy()
                        stdevPEList.append(copySD)
                    indx = 1

        #print(len(avgTSList))  
        for pltIndx in range(len(avgTSList)):
            if(pltIndx%2 == 0):
                plt.errorbar(avgTSList[pltIndx],avgPEList[pltIndx],yerr = stdevPEList[pltIndx], ecolor = 'b',color = 'k',ms = 3,capsize = 3,label = labels[fileIndx])
            else:
                plt.errorbar(avgTSList[pltIndx],avgPEList[pltIndx],yerr = stdevPEList[pltIndx], ecolor = 'r',color = 'k',ms = 3,capsize = 3,label = labels[fileIndx])
    plt.xlabel('Timestep')
    plt.ylabel('Potential Energy (KT)')
    #plt.legend()
    plt.show()

main()	

