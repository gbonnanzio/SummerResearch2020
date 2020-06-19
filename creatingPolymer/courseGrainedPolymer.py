#! Program Files (x86)/PythonDownload
import numpy as py
import math

def createPolymers(numPolymers,bondAngle,bondLength,numAtomsInChain):
    
    # cg polymer is being built along the y axis in the xy plane
    # additional molecules are built upwards in the z direction

    pi = math.pi
    # angle to be used for calculations
    angleUsed = (bondAngle/2)*(pi/180)
     
    # list of all atom info broken down by molecule [[[type,x,y,z],[type,x,y,z]...][[type,x,y,z],...]]
    allAtoms = []
    # same format-ish as allAtoms
    allBonds = []
    allAngles = []
    atomsInChain = [] #  [type, x, y, z] coords 
    bondsInChain = [] #  [type, atom#, atom#]
    anglesInChain = [] # [type, end, middle, end] (these are atom numbers)
   
    currAtom = 0 # indicates atom you are about to create
    for numChain in range(numPolymers):

        # create first bead
        currAtom += 1
        # started at 1 so not on boundary of sim box
        atomsInChain.append([1,1.0,1.0,1.0+numChain*2*bondLength])

        while(len(atomsInChain) < numAtomsInChain):
            
            currAtom += 1
            #if you are about to create an even numbered atom
            if(currAtom%2 == 0):
                cxAdd = bondLength*math.cos(angleUsed)
            else:
                cxAdd = -bondLength*math.cos(angleUsed)
            
            cyAdd = bondLength*math.sin(angleUsed)
            
            # create a new bead
            
            # get coords of the last bead added
            if(currAtom%numAtomsInChain == 0):
                newAtomCoords = atomsInChain[numAtomsInChain-2].copy()
            else:
                newAtomCoords = atomsInChain[(currAtom%numAtomsInChain)-2].copy()
            # no change [0], all same type
            newAtomCoords[1] = newAtomCoords[1] + cxAdd 
            newAtomCoords[2] = newAtomCoords[2] + cyAdd
            #no change in z (all in xy plane)
            atomsInChain.append(newAtomCoords)

            # excludes the first atom created in each molecule
            if(currAtom%(numAtomsInChain) != 1):
                bondsInChain.append([1,currAtom-1,currAtom])
            
            # excludes the first 2 atoms in each chain (0 for last atom in chain)
            if(currAtom%numAtomsInChain == 0 or currAtom%numAtomsInChain > 2):
                anglesInChain.append([1,currAtom-2,currAtom-1,currAtom])
        
        # add data to overall lists and reset
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
    
    # variables to change
    outFile = open('polymerData.data','w') # where you want to store the coords of the polymer
    numPolymers = 2 # how many polymers you're creating
    numAtomsPerMol = 30 # number of beads in each polymer
    polymerBL = 1.0 # bond length between beads
    polymerAngle = 109.5 # angle between three beads
    cutoffLength = 3.5 # cutoff distance for LJ Pot. to be calculated

    # atomsList has coords of every atom, bondsList has bonds between every atom, anglesList has angles between atoms
    # coords are in their own lists by molecule
    atomsList, bondsList, anglesList = createPolymers(numPolymers,polymerAngle,polymerBL,numAtomsPerMol)
    lengthBox = numAtomsPerMol + cutoffLength 

    #write to our file all the info we need
    outFile.write('# Polymer Data file\n')
    outFile.write(str(len(atomsList)*len(atomsList[0])) + '\t atoms\n')
    outFile.write(str(len(bondsList)*len(bondsList[0])) + '\t bonds\n')
    outFile.write(str(len(anglesList)*len(anglesList[0])) + '\t angles\n\n')
    outFile.write('1\tatom types\n1\tbond types\n1\tangle types\n\n')
    outFile.write('0.0000\t' + str(lengthBox)+' xlo xhi\n')
    outFile.write('0.0000\t' + str(lengthBox)+' ylo yhi\n')
    outFile.write('0.0000\t' + str(lengthBox)+' zlo zhi\n')
    outFile.write('Masses\n\n1\t1.0\n')

    
    # write atom coords in format atom# mol# type# x y z
    outFile.write('\nAtoms\n\n')
    totAtoms = 1 # number of atoms we've created
    for polymerIndx in range(len(atomsList)):
        tmpAtomList = atomsList[polymerIndx] # all atoms in one polymer
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
    
    # write bond info in format bond# type# atom1 atom2 
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

    
    # write angle info in format angle# type# atom1 atom2 atom3
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

# call functions
main()

