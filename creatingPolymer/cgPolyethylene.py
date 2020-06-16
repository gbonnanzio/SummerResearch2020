#! Program Files (x86)/PythonDownload
import numpy as py
import math

def createPolymer():
    
    pi = math.pi
    #carbon carbon carbon angle
    ccca = pi/3
    
    # rotation matrix
    #R = [[math.cos(ccca), math.sin(ccca), 0],[-math.sin(ccca),math.cos(ccca),0],[0,0,1]]

    #mirMatrix = [[0,1,0],[-1, 0, 0],[0,0,1]]

    #RH = [[1,0,0],[0,math.cos(pi/4),math.sin(pi/4)],[0,-math.sin(pi/4),math.cos(pi/4)]]

    ccBond = 2
    chBond = 1

    atoms = [] # atom x y z coords
    types = [] # 1 for C, 2 for H
    bonds = [] # type -> 1 for C-C, 2 for C-H. [type, atom#, atom#]
    angles = [] # type -> [end, middle, end] (atom number)
    dihedrals = []
    #create first carbon (type 1) at origin
    atoms.append([0,0,0])
    types.append(1)
    
    #create H in yz plane
    atoms.append([-chBond*math.cos(pi/4), 0, chBond*math.sin(pi/4)])
    types.append(2)
    bonds.append([2,1,2])
    
    #create another H in yz plane but now -z
    atoms.append([-chBond*math.cos(pi/4),0, -chBond*math.sin(pi/4)])
    types.append(2)
    bonds.append([2,1,3])
    angles.append([2,1,3])

    currAtom = 3
    while len(atoms) < 30:
        currAtom += 1
        if(currAtom%2 == 0):
            cxAdd = ccBond*math.cos(pi/3)
            hxAdd = chBond*math.cos(pi/4)
        else:
            cxAdd = -ccBond*math.cos(pi/3)
            hxAdd = -chBond*math.cos(pi/4)
        
        cyAdd = ccBond*math.sin(pi/3)
        
        #create a new carbon
        prevAtomCoords = atoms[currAtom - 4].copy()
        prevAtomCoords[0] = prevAtomCoords[0] + cxAdd 
        prevAtomCoords[1] = prevAtomCoords[1] + cyAdd
        atoms.append(prevAtomCoords)
        types.append(1)
        bonds.append([1,currAtom-3,currAtom])
        
        #create c-c-c angle
        if(currAtom > 6):
            angles.append([currAtom,currAtom-3,currAtom-6])
        #create dihedral
        if(currAtom > 9):
            dihedrals.append([currAtom-9,currAtom-6,currAtom-3,currAtom])
        # add first H to new C
        currAtom += 1
        currHCoords = [prevAtomCoords[0]+hxAdd,prevAtomCoords[1],chBond*math.sin(pi/4)]
        atoms.append(currHCoords)
        types.append(2)
        bonds.append([2,currAtom,currAtom-1])
        # add second H to new C
        currAtom += 1
        currHCoords[2] = -currHCoords[2]
        atoms.append(currHCoords)
        types.append(2)
        bonds.append([2,currAtom,currAtom-2])
        angles.append([currAtom,currAtom-2,currAtom -1])
    
    return atoms, types, bonds, angles, dihedrals


def main():
    outFile = open('polymerOutput.txt','w')
    atomsList, typesList, bondsList, anglesList, dihedralsList = createPolymer()
    outFile.write(str(len(atomsList)) + '\t atoms\n')
    outFile.write(str(len(bondsList)) + '\t bonds\n')
    outFile.write(str(len(anglesList)) + '\t angles\n')
    outFile.write(str(len(dihedralsList)) + '\t dihedrals\n')
    
    outFile.write('\nAtoms\n\n')
    for indx in range(len(atomsList)):
        outFile.write(str(indx+1) + '\t' + str(typesList[indx]) + '\t')
        for coord in atomsList[indx]:
            coord = round(coord,10)
            outFile.write(str(coord)+'\t')
        outFile.write('\n')
    
    outFile.write('\nBonds\n\n')
    for indx in range(len(bondsList)):
        outFile.write(str(indx+1) + '\t')
        for element in bondsList[indx]:
            #round(element,5)
            outFile.write(str(element)+'\t')
        outFile.write('\n')
 
    outFile.write('\nAngles\n\n')
    for indx in range(len(anglesList)):
        outFile.write(str(indx+1) + '\t')
        for element in anglesList[indx]:
            #round(element,5)
            outFile.write(str(element)+'\t')
        outFile.write('\n')

    outFile.write('\nDihedrals\n\n')
    for indx in range(len(dihedralsList)):
        outFile.write(str(indx+1) + '\t')
        for element in dihedralsList[indx]:
            #round(element,5)
            outFile.write(str(element)+'\t')
        outFile.write('\n')

main()

