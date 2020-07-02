#! Program Files (x86)/PythonDownload/python
import matplotlib.pyplot as plt
import math 
import numpy as np
import statistics as stats
import sys
import glob
#os.chdir("/Users/gbonn/Summer_Research_2020")
from readFiles import (readIndividualDump, readFormattedDump,readMyDumpFileForXYZ, formatDumpFileForXYZ, readMyVMDradGyrData, readMyCOMFile)
#os.chdir("/Users/gbonn/Summer_Research_2020/creatingPolymer")

def determineDis(atom1, atom2,boxLen):
    # determines the distance between two atoms in a PBC
    # cubic box. Usually used to determine end to end distance
    # atom1/atom2 -> list, [x,y,z]
    # boxLen -> float, length of cubic box

    halfBoxLen = boxLen/2.0

    x1 = atom1[0]
    x2 = atom2[0]
    y1 = atom1[1]
    y2 = atom2[1]
    z1 = atom1[2]
    z2 = atom2[2]
            
    diffX = x2-x1
    if(diffX>halfBoxLen):
        diffX = diffX - boxLen
    if(diffX<-halfBoxLen):
        diffX = diffX + boxLen
    diffY = y2-y1
    if(diffY>halfBoxLen):
        diffY = diffY - boxLen
    if(diffY<-halfBoxLen):
        diffY = diffY + boxLen
    diffZ = z2-z1
    if(diffZ>halfBoxLen):
        diffZ = diffZ - boxLen
    if(diffZ<-halfBoxLen):
        diffZ = diffZ + boxLen
 
    rSq = (diffX**2)+(diffY**2)+(diffZ**2)
    finalR = math.sqrt(rSq)

    return finalR


def determineRadGyr(xCoords,yCoords,zCoords,COM,boxLen):
    # determine the center of mass for polymer with equally
    # weighted monomers and then find the average distance
    # of atoms away from this (distance^2)
    halfBoxLen = boxLen/2.0
    numAtoms = len(xCoords)
    #print(numAtoms)
    cmX = COM[0]%boxLen
    cmY = COM[1]%boxLen
    cmZ = COM[2]%boxLen
    allRSq = []

    for i in range(numAtoms):
        x = xCoords[i]
        y = yCoords[i]
        z = zCoords[i]

        diffX = x-cmX
        if(diffX>halfBoxLen):
            diffX = diffX - boxLen
        if(diffX<-halfBoxLen):
            diffX = diffX + boxLen
        diffY = y-cmY
        if(diffY>halfBoxLen):
            diffY = diffY - boxLen
        if(diffY<-halfBoxLen):
            diffY = diffY + boxLen 
        diffZ = z-cmZ
        if(diffZ>halfBoxLen):
            diffZ = diffZ - boxLen
        if(diffZ<-halfBoxLen):
            diffZ = diffZ + boxLen

        rSq = (diffX**2)+(diffY**2)+(diffZ**2)
        allRSq.append(rSq)
    
    avgRSq = stats.mean(allRSq)
    
    return avgRSq




def binEndToEndDistance(numAtomsInMol,totMols,boxWidth,nBins,numTStoSkip,allFiles,avgOutFile):
    # calculates the end to end distance (Ree) of every polymer at every timestep and creates
    # a frequency distribution of all of these lengths. Also outputs the avg and SD of all Ree at
    # each timestep to the file avgOutFile
    # totAtoms -> int TOTAL number of atoms in system, totMols -> int number of polymer chains in system    
    # boxWidth -> float length of cubic box, nBins -> number of bins your histogram uses
    # numTStoSkip -> skips this many timesteps in analysis
    # allCoords -> list of all atoms in format [[[x,y,z],[x,y,z]...n],[[x,y,z],[x,y,z]...n]...ts]
    # allCoords is a list of all n atom coords at each time step (atoms are not broken down by respective polymer)
    # avgOutFile -> string of file name to output avg Ree at each ts

    # every Ree at every timestep
    allRee = []
    # avg and SD of Ree at every timestep
    avgPolRee = []
    SDPolRee = []

    # outer for loop iterates through timesteps
    for ts in range(numTStoSkip,len(allFiles)):
        sys.stdout.write('Current TS: ' + str(ts) + '\n')
        sys.stdout.flush()
        # list of polymers at ts
        tsList = readIndividualDump(allFiles[ts],numAtomsInMol,totMols)
    
        # Ree of each polymer at given timestep (gets reset to null each ts)
        tmpReeList = [] 
        
        # iterate through number of polymers in system
        for polNum in range(len(tsList)):
            currPol = tsList[polNum]
            # determine current Ree
            tmpR = determineDis(currPol[0],currPol[-1],boxWidth)
            allRee.append(tmpR)
            tmpReeList.append(tmpR)
            # update which polymer we're looking at
        
        # finding avg and stdev of Ree among all polymers at given timestep
        avgPolRee.append(stats.mean(tmpReeList))
        SDPolRee.append(stats.stdev(tmpReeList))
        
    # output avg and stdev to outputfile
    fileObj = open(avgOutFile,'w') 
    sys.stdout.write('WRITING AVERAGES TO FILE\n')
    sys.stdout.flush()
    for ts in range(len(avgPolRee)):
        fileObj.write(str(avgPolRee[ts])+ ' ' + str(SDPolRee[ts]) +'\n')
    fileObj.close()

    sys.stdout.write('PLOTTING DATA\n')
    sys.stdout.flush()
    
    # bin Ree to histogram
    maxR = max(allRee)
    minR = 0
    rRange = np.linspace(minR,maxR,nBins)
    rStep = rRange[1] - rRange[0]
    
    # initialize all bins to 0
    binnedRee = [0]*len(rRange)
    for ts in range(len(allRee)):
        # determine bin indx and add 1 to current size
        binnedIndx = int((allRee[ts]-minR)//rStep)
        binnedRee[binnedIndx] += 1 
    # adjust R to plot in center of bin
    binnedR = [x+0.5*rStep for x in rRange]
    #find the frequency of each bin 
    totalBinned = float(len(allRee))
    freqRee = [num/totalBinned for num in binnedRee]
    
    # plot end to end distance distribution
    plt.plot(binnedR,freqRee,'.k',label = 'P($R_{ee}$)')
    plt.plot(binnedR,freqRee,'k',lw = 1,label = None)
    plt.xlabel('$R_{ee}$ [d]')
    plt.ylabel('P($R_{ee}$)') 
    plt.legend()
    plt.show()
    
def binRadGyr(numAtomsInMol,totMols,centersOfMass,boxWidth,nBins,numTStoSkip,allFiles,avgOutFile):
    # calculates the radius of gyration (Rg) of every polymer at every timestep and creates
    # a frequency distribution of all of these lengths. Also outputs the avg and SD of all Rg at
    # each timestep to the file avgOutFile
    # totAtoms -> int TOTAL number of atoms in system, totMols -> int number of polymer chains in system    
    # boxWidth -> float length of cubic box, nBins -> number of bins your histogram uses
    # numTStoSkip -> skips this many timesteps in analysis
    # allCoords -> list of all atoms in format [[[x,y,z],[x,y,z]...n],[[x,y,z],[x,y,z]...n]...ts]
    # allCoords is a list of all n atom coords at each time step (atoms are not broken down by respective polymer)
    # avgOutFile -> string of file name to output avg Rg at each ts

    # every radius of gyration of every polymer at every timestep (one long list)
    allRadGyr = []
    # list of avg and SD of polymer Rg at each timestep
    avgPolRg = [] 
    SDPolRg = []
    
    # iterates through all timesteps
    for ts in range(numTStoSkip,len(allFiles)):
        sys.stdout.write('Current TS: ' + str(ts) + '\n')
        sys.stdout.flush()
        # list of polymers at ts
        tsList = readIndividualDump(allFiles[ts],numAtomsInMol,totMols)
        
        #radius of gyration of each polymer at a given time step (gets reset to null each ts)
        tmpRgList = [] 
        
        #iterates through all polymers at a given ts
        for polNum in range(len(tsList)):
            # get and arrange coords of ONE polymer at a ts (last index is exclusive)
            tmpMolList = tsList[polNum]
            # sort coords into 3 lists of x, y, z
            currX = []
            currY = []
            currZ = []
            for atom in tmpMolList:
                currX.append(atom[0])
                currY.append(atom[1])
                currZ.append(atom[2])
            # get center of mass of that polymer (out file from lammps)
            tmpCM = centersOfMass[ts][polNum]
            # Rg of ONE polymer at a ts
            currRg = determineRadGyr(currX,currY,currZ,tmpCM,boxWidth)
            allRadGyr.append(currRg)
            tmpRgList.append(currRg)   
            # go to next molecule indices
            #first += numIndxs 
            #last += numIndxs
        
        # determine avg and SD of all polymer Rg's at a given ts
        avgPolRg.append(stats.mean(tmpRgList))
        SDPolRg.append(stats.stdev(tmpRgList))
    
    # write avg and SD to output file 
    fileObj = open(avgOutFile,'w')
    sys.stdout.write('WRITING AVERAGES TO FILE\n')
    sys.stdout.flush()
    for ts in range(len(avgPolRg)):
        fileObj.write(str(avgPolRg[ts]) + ' ' + str(SDPolRg[ts]) + '\n')
    fileObj.close()
    sys.stdout.write('PLOTTING DATA\n')
    sys.stdout.flush()
    # bin all Rg
    maxR = max(allRadGyr)
    minR = 0
    rRange = np.linspace(minR,maxR,nBins)
    rStep = rRange[1] - rRange[0]
    # initialize counts to 0
    binnedRadGyr = [0]*len(rRange) 
    for ts in range(len(allRadGyr)):
        # determine indx and increment count by 1
        binnedIndx = int((allRadGyr[ts]-minR)//rStep)
        binnedRadGyr[binnedIndx] += 1 
    
    # adjust R to plot in center of bin
    binnedR = [x+0.5*rStep for x in rRange] #plot the pt in center of each bin
    # determine frequency of each bin
    totalBinned = float(len(allRadGyr)) 
    freqGyr = [num/totalBinned for num in binnedRadGyr] 
    
    # plot radius of gyration
    plt.plot(binnedR,freqGyr,'.k',label = 'P($R_g^2$)')
    plt.plot(binnedR,freqGyr,'k',label = None)
    plt.xlabel('Radius of Gyration <$R_g^2$> [$d^2$]')
    plt.xlim(0,boxWidth/2)
    plt.ylabel('P($R_g^2$)') 
    plt.legend()
    plt.show()
    

def radGyrVMD(VMDFile,numBins):
    
    #plot radius of gyration from data from VMD
    vmdRadGyrList = readMyVMDradGyrData(VMDFile)
    maxR = max(vmdRadGyrList)
    minR = 0
    rRange = np.linspace(minR,maxR,numBins)
    rStep = rRange[1] - rRange[0]
    binnedRadGyrVMD = [0]*len(rRange)
    for ts in range(len(vmdRadGyrList)):
        binnedIndx = int((vmdRadGyrList[ts]-minR)//rStep)
        binnedRadGyrVMD[binnedIndx] += 1 
        
    binnedR = [x+0.5*rStep for x in rRange]
    totalBinned = float(len(vmdRadGyrList))
    freqGyrVMD = [num/totalBinned for num in binnedRadGyrVMD]
    print(sum(freqGyrVMD))
    # plot radius of gyration
    plt.plot(binnedR,freqGyrVMD,'sb',label = 'P(r) generated by VMD')
    plt.plot(binnedR,freqGyrVMD,'b',lw = 1, label = None)


def main():
   
    # Files to get lammps data from
    polymerDumpXYZFile = 'polymerCoordsXYZ.lmpdump'
    centerOfMassFile = 'runLammps/com.out'
    
    # Files to output avg Ree and avg Rg at each timestep to
    ReeOutFile = 'Calc_Ree.txt'
    RgOutFile = 'Calc_Rg.txt'

    # simulation variables to change
    numMols = 50
    atomsPerMol = 20
    numAtoms = numMols*atomsPerMol
    boxLength = 30.0
    totalBins = 45
    tsToSkip = 0
    
    # get list of coordinate files from directory
    fileList = glob.glob('./runLammps/dump.coords*')
    
    # get centers of mass of each polymer at each timestep
    comList = readMyCOMFile(centerOfMassFile,numMols)
    
    # create a probability distribution of end to end distances of polymers
    sys.stdout.write('BINNING END TO END DISTANCE\n')
    sys.stdout.flush()
    binEndToEndDistance(atomsPerMol,numMols,boxLength,totalBins,tsToSkip,fileList,ReeOutFile)
    
    # create a probability distribution of the radius of gyration of polymers
    sys.stdout.write('BINNING RADIUS OF GYRATION\n')
    sys.stdout.flush()
    binRadGyr(atomsPerMol,numMols,comList,boxLength,totalBins,tsToSkip,fileList,RgOutFile)
    
main()
    

