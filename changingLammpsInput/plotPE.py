#! Program Files (x86)/PythonDownload
import matplotlib.pyplot as plt
import statistics as stat 
from readFiles import readMyPeOUTFile

def checkInRange(min1,max1,min2,max2):
    check1 = min1 < max2 and min1 > min2
    check2 = max1 < max2 and max1 > min2
    finalCheck = check1 or check2
    return finalCheck 

def main():
    fileList = ['thermo.125.OUT','thermo.25.OUT','thermo.50.OUT']
    colors = ['ro','bs','g^']
    labels = ['ρ = 0.125','ρ = 0.25','ρ = 0.50']
    blockAvg = 10
    avgTSList = []
    tmpTSList = []
    avgPEList = []
    tmpPEList = []
    stdevPEList = []
    tmpSDList = []
    for fileIndx in range(1):
        ts, pe = readMyPeOUTFile(fileList[fileIndx])
        print(len(ts))
        first = 0
        last = blockAvg
        indx = 0
        for tsIndx in range(len(ts)//blockAvg):
            tmpTS = stat.mean(ts[tsIndx*blockAvg:(tsIndx+1)*blockAvg])
            tmpPE = stat.mean(pe[tsIndx*blockAvg:(tsIndx+1)*blockAvg])
            tmpPEstdev = stat.stdev(pe[tsIndx*blockAvg:(tsIndx+1)*blockAvg])
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
                inRange = checkInRange(prevMin,prevMax,currMin,currMax)  
                if(inRange):
                    tmpTSList.append(tmpTS)
                    tmpPEList.append(tmpPE)
                    tmpSDList.append(tmpPEstdev)
                    indx += 1
                elif(tsIndx == len(ts)-1):
                    copyTS = tmpTSList.copy()
                    tmpTSList = []
                    avgTSList.append(copyTS)
                    copyPE = tmpPEList.copy()
                    tmpPEList = []
                    avgPEList.append(copyPE)
                    copySD = tmpSDList.copy()
                    tmpSDList = []
                    stdevPEList.append(copySD)
                else:
                    copyTS = tmpTSList.copy()
                    tmpTSList = []
                    avgTSList.append(copyTS)
                    copyPE = tmpPEList.copy()
                    tmpPEList = []
                    avgPEList.append(copyPE)
                    copySD = tmpSDList.copy()
                    tmpSDList = []
                    stdevPEList.append(copySD)
                    indx = 0

            first += blockAvg
            last += blockAvg

            
        print(avgTSList)
        for pltIndx in range(len(avgTSList)):
            if(pltIndx%2 == 0):
                plt.errorbar(avgTSList[pltIndx],avgPEList[pltIndx],yerr = stdevPEList[pltIndx], ecolor = 'b',ms = 3,label = labels[fileIndx])
            else:
                plt.errorbar(avgTSList[pltIndx],avgPEList[pltIndx],yerr = stdevPEList[pltIndx], ecolor = 'g',ms = 3,label = labels[fileIndx])
    plt.xlabel('Timestep')
    plt.ylabel('Potential Energy (KT)')
    #plt.legend()
    plt.show()

main()	

