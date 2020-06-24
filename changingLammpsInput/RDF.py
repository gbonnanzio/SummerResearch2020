#! Program Files (x86)/PythonDownload/python
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
from readFiles import readMyRdfFile
from readFiles import readMyDumpFileForXYZ

################## determinePE ##############################################################

def determineR(coordsTS,lengthBox,allBinned,binWidth):
    
    # takes the coordinates of a certain number of atoms 
    # determines the distance between every possible pair and the Lennard Jones potential 
    # and returns the radii with the pair (assumes periodic boundary conditions)
    
    halfBox = lengthBox/2.0
    nAtoms = len(coordsTS)
    allR = []
    allPE = []
    
    # j is adjusted so that you never compare any two atoms twice 
    for i in range(nAtoms-1):
        x1 = coordsTS[i][0]
        y1 = coordsTS[i][1]
        z1 = coordsTS[i][2]
        for j in range(i+1,nAtoms):
            
            x2 = coordsTS[j][0]
            y2 = coordsTS[j][1]
            z2 = coordsTS[j][2]

            # implement PBC for box length 10 
            diffX = x2-x1
            if(diffX>halfBox):
                diffX = diffX - lengthBox
            if(diffX<-halfBox):
                diffX = diffX + lengthBox
            diffY = y2-y1
            if(diffY>halfBox):
                diffY = diffY - lengthBox
            if(diffY<-halfBox):
                diffY = diffY + lengthBox
            diffZ = z2-z1
            if(diffZ>halfBox):
                diffZ = diffZ - lengthBox
            if(diffZ<-halfBox):
                diffZ = diffZ + lengthBox

            # calculate r and potential energy 
            rSq = (diffX**2)+(diffY**2)+(diffZ**2)
            r = math.sqrt(rSq)
            binIndx =  int(r//binWidth)
            allBinned[binIndx] = allBinned[binIndx] + 2.0
            #tempPE = 4.0*(((1/rSq)**6)-((1/rSq)**3))
            #allPE.append(tempPE)
        
    return allBinned


######################## binning ##############################################################

def binning(binSpace, currRValues,allBinned):
    
    # takes the current r and PE lists for the current time step 
    # as well as a list of all the past binned r and PE values 
    # and the number of bins 
    # finds what bin r is in and averages in PE to running total
    # returns compiled bin and the specified 
   
    # if we are creating the final bins for the first time
    
    # bin the current (r,PE) pairs   
    for atomNum in range(len(currRValues)):
        tempR = currRValues[atomNum]
        # iterate from back of list starting at the largest r until
        # we find spot for tempR in bin 
        binIndx = int(tempR//binSpace)
        # if initializing this bin
        allBinned[binIndx] = allBinned[binIndx] + 2.0
            
    return allBinned

######################### radialDistFun ##########################################################

def radialDistFun(binWidth,allRad,allCts,runs,totAtoms,lenBox):
    
    # takes a binSize, all of the radii from all time steps, the bin count for each
    # radius and the number of time steps performed
    # creates the radial distribution function g(r)
    radDistVal = []
    n = []
    nID = []
    meanDens = (totAtoms)/(lenBox**3)
    for i in range(len(allRad)):
        n.append(allCts[i]/(totAtoms*runs))
        tempNId = 4.0*math.pi*(allRad[i]**2)*binWidth*meanDens
        nID.append(tempNId)
        radDistVal.append(n[i]/nID[i])

    return radDistVal



################################### MAIN ##########################################################

def plotRDF(fileLocs,numFiles,lengthOfBoxes):

    #specify where files are located
    #coordFile = fileLoc #contains all the coords of the atoms at every time stamp
    
    # initialize necessary variables
    numAtoms = 500
    totBins = 150
    finalR = []
    finalG = []
    for plotIndx in range(numFiles):
        
        boxLen = lengthOfBoxes[plotIndx]
        minR = 0.0
        maxR = (boxLen/2.0)*math.sqrt(3) # for box length 10 with PBC
        binSize = (maxR-minR)/(totBins - 1)
        rRange = np.linspace(minR, maxR, totBins)
        
        binnedPts = [0 for i in rRange] # how many pairs we added to each bin
    
        sys.stdout.write('\nSTARTING DUMPFILE NUMBER ' + str(plotIndx) + '\n')
        sys.stdout.flush()
        allCoords = readMyDumpFileForXYZ(fileLocs[plotIndx],9, numAtoms,0)
        currTs = 1
        for ts in range(len(allCoords)):
            
            binnedPts = determineR(allCoords[ts],boxLen,binnedPts,binSize)
        
            sys.stdout.write('DUMP: ' + str(plotIndx) + ' TS: ' + str(currTs)+'\n')
            sys.stdout.flush()
            currTs += 1
            if(currTs > 1000):
                break
        #print(sum(binnedPts))
        tmpRRange = [i+0.5*binSize for i in rRange]
        finalR.append(tmpRRange)
        # generate radial distribution function 
        finalG.append(radialDistFun(binSize,tmpRRange,binnedPts, currTs,numAtoms,boxLen))


    # plot the RDF 
    #calculated
    colors = ['ro','bs','g^']
    labels = ['ρ = 0.125','ρ = 0.25','ρ = 0.50']
    for plotIndx in range(numFiles):
        plt.plot(finalR[plotIndx],finalG[plotIndx],colors[plotIndx],ms = 3, label = labels[plotIndx])
        plt.plot(finalR[plotIndx],finalG[plotIndx],colors[plotIndx][0],lw = 1,label = None)
    
    onesList = [1]*(totBins)
    plt.plot(finalR[0],onesList,'k',ls = '--',label = None)
    plt.xlim(0,5)
    #plt.ylim(0, 3.0)
    plt.xlabel("Distance (d)")
    plt.ylabel("g(r)")
    plt.legend()
    #plt.legend(["Calculated", "VMD","LAMMPS"])
    plt.show()
    
    
def main():
    fileLocList = ['coords.125.lmpdump','coords.25.lmpdump','coords.50.lmpdump']
    boxLengthList = [15.874010519681994,12.599210498948732,10.0]
    plotRDF(fileLocList,len(fileLocList),boxLengthList)

#run functions    
main()


