#! Program Files (x86)/PythonDownload
import numpy as py
import math

def createPolymers(numPolymers,bondAngle,bondLength,numAtomsInChain):
    
    pi = math.pi
    #carbon carbon carbon angle
    angleUsed = (bondAngle/2)*(pi/180)
    
    ccBL = 1.1 #arbitrary change later
    allAtoms = []
    allBonds = []
    allAngles = []
    atomsInChain = [] # type -> 1 for C [type, x, y, z] coords 
    bondsInChain = [] # type -> 1 for C-C, [type, atom#, atom#]
    anglesInChain = [] # type -> 1 for C-C-C [type, end, middle, end] (these are atom numbers)
   
    currAtom = 0 # num indicates atom you are about to create
    currBond = 0
    currAngle = 0
    for numChain in range(numPolymers):

        #create first C at origin
        currAtom += 1
        atomsInChain.append([1,0,0,numChain*2*bondLength])

        while(len(atomsInChain) < numAtomsInChain):
            
            currAtom += 1

            if(currAtom%2 == 0):
                cxAdd = bondLength*math.cos(angleUsed)
            else:
                cxAdd = -bondLength*math.cos(angleUsed)
            
            cyAdd = bondLength*math.sin(angleUsed)
            
            #create a new carbon

            #get coords of the last carbon
            newAtomCoords = atomsInChain[currAtom%numAtomsInChain - 2].copy()
            # no change [0], all same type
            newAtomCoords[1] = newAtomCoords[1] + cxAdd 
            newAtomCoords[2] = newAtomCoords[2] + cyAdd
            #no change in z (all in xy plane)
            atomsInChain.append(newAtomCoords)
            if(currAtom%(numAtomsInChain) != 1):
                bondsInChain.append([1,currAtom-1,currAtom])
            
            # if more than 2 C, create c-c-c angle
            if(currAtom>2 and (currAtom%numAtomsInChain == 0 or currAtom%numAtomsInChain > 2)):
                anglesInChain.append([1,currAtom-2,currAtom-1,currAtom])
        tmpAtomsCopy = atomsInChain.copy()
        allAtoms.append(tmpAtomsCopy)
        atomsInChain = []
        tmpBondsCopy = bondsInChain.copy()
        allBonds.append(tmpBondsCopy)
        bondsInChain = []
        tmpAnglesCopy = anglesInChain.copy()
        allAngles.append(tmpAnglesCopy)
        anglesInChain = []

    
    return allAtoms, allBonds, allAngles


def main():
    outFile = open('polymerData.data','w')
    numPolymers = 4
    atomsList, bondsList, anglesList = createPolymers(numPolymers,109.5,1.1,30)
    #print(atomsList)
    outFile.write('# Polymer Data file\n')
    outFile.write(str(len(atomsList)*len(atomsList[0])) + '\t atoms\n')
    outFile.write(str(len(bondsList)*len(bondsList[0])) + '\t bonds\n')
    outFile.write(str(len(anglesList)*len(anglesList[0])) + '\t angles\n\n')
    outFile.write('1\tatom types\n1\tbond types\n1\tangle types\n\n')
    outFile.write('-10.0000\t40.0000 xlo xhi\n')
    outFile.write('-10.0000\t40.0000 ylo yhi\n')
    outFile.write('-10.0000\t40.0000 zlo zhi\n\n')
    outFile.write('Masses\n\n1\t1.05\n')


    outFile.write('\nAtoms\n\n')
    totAtoms = 1
    for polymerIndx in range(len(atomsList)):
        tmpAtomList = atomsList[polymerIndx] #atoms in one polymer
        for atomIndx in range(len(tmpAtomList)):
            outFile.write(str(totAtoms) + '\t')
            totAtoms += 1
            for coordIndx in range(len(tmpAtomList[atomIndx])):
                if coordIndx == 0:
                    outFile.write(str(polymerIndx+1) + '\t'+ str(tmpAtomList[atomIndx][coordIndx])+'\t')
                else:
                    tmpCoord = "{:.9f}".format(tmpAtomList[atomIndx][coordIndx])
                    outFile.write(tmpCoord +'\t')
            outFile.write('\n')
    
    outFile.write('\nBonds\n\n')
    totBonds = 1
    for polymerIndx in range(len(bondsList)):
        tmpBondList = bondsList[polymerIndx] #atoms in one polymer
        for bondIndx in range(len(tmpBondList)):
            outFile.write(str(totBonds) + '\t')
            totBonds += 1
            for element in tmpBondList[bondIndx]:
                outFile.write(str(element)+'\t')
            outFile.write('\n')


    outFile.write('\nAngles\n\n')
    totAngles = 1
    for polymerIndx in range(len(anglesList)):
        tmpAngleList = anglesList[polymerIndx] #atoms in one polymer
        for angleIndx in range(len(tmpAngleList)):
            outFile.write(str(totAngles) + '\t')
            totAngles += 1
            for element in tmpAngleList[angleIndx]:
                outFile.write(str(element)+'\t')
            outFile.write('\n')


main()

