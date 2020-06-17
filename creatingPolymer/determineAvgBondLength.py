#! Program Files (x86)/PythonDownload/python
import statistics as stat 
import os

#os.chdir("/Users/gbonn/Summer_Research_2020")
from readFiles import readMyDumpFileForAvgVal
#os.chdir("/Users/gbonn/Summer_Research_2020/creatingPolymer")



def determineAvgBondAndAngle():
    
    # takes the coordinates of a certain number of atoms 
    # determines the distance between every possible pair and the Lennard Jones potential 
    # and returns the radii with the pair (assumes periodic boundary conditions)
    
    bondFile = 'bonds.lmpdump'
    angleFile = 'angles.lmpdump'
    allBonds = readMyDumpFileForAvgVal(bondFile,9,116)
    allAngles = readMyDumpFileForAvgVal(angleFile,9,112)
    #print(allBonds,allAngles)
    print(stat.mean(allBonds),stat.mean(allAngles))
    print(stat.stdev(allBonds),stat.stdev(allAngles))
    #numCounted = 1
    #currBondMean = 0
    #for bondList in allBonds:
     #   currBondMean = (mean(bondList) + numCounted*currBondMean)/numCounted
      #  numCounted += 1

    #numCounted = 1
    #currAngleMean = 0
    #for angleList in allAngles:
     #   currAngleMean = (mean(angleList) + numCounted*currAngleMean)/numCounted
      #  numCounted += 1
    
    #print(currBondMean,currAngleMean)

def main():
    determineAvgBondAndAngle()

main()
    
