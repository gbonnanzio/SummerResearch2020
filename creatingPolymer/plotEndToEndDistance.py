#! Program Files (x86)/PythonDownload/python
import statistics as stat 
import os
import matplotlib.pyplot as plt
import numpy as np

#os.chdir("/Users/gbonn/Summer_Research_2020")
from readFiles import readMyDumpFileForAvgVal
#os.chdir("/Users/gbonn/Summer_Research_2020/creatingPolymer")



def binEndToEndDistance():
    
    # takes the coordinates of a certain number of atoms 
    # determines the distance between every possible pair and the Lennard Jones potential 
    # and returns the radii with the pair (assumes periodic boundary conditions)
    
    bondFile = 'bonds.lmpdump'
    angleFile = 'angles.lmpdump'
    totBonds = 58
    totMols = 2
    allBonds = readMyDumpFileForAvgVal(bondFile,9,totBonds)
    #print(allBonds)
    allREE = []
    currRee = 0
    for tsList in allBonds:
        first = 0
        last = int(totBonds/totMols - 1)
        for polNum in range(totMols):
            tmpList = (tsList[first:last])
            tmpRee = (sum(tmpList)**2)/len(tmpList)
            first += last
            last += last
            allREE.append(tmpRee)
            

    #we have a list allREE of all the lengths of the polymer at every timestep
    maxR = max(allREE)
    print(maxR)
    minR = min(allREE)
    print(minR)
    nBins = 100
    rRange = np.linspace(minR,maxR,nBins)
    rStep = rRange[1] - rRange[0]
    binnedREE = [0]*len(rRange)
    for ts in range(len(allREE)):
        binnedIndx = int((allREE[ts]-minR)//rStep)
        binnedREE[binnedIndx] += 1 
        
    binnedR = [x+0.5*rStep for x in rRange]
    totalBinned = float(len(allREE))
    freqREE = [num/totalBinned for num in binnedREE]
    print(sum(freqREE))
    #print(binnedR)
    plt.bar(binnedR,freqREE,width = rStep)
    plt.xlabel('Bond Length')
    plt.ylabel('Frequency')
    #plt.ylim(0,.15)
    #print(bins)
    plt.show()


def main():
    binEndToEndDistance()

main()
    
