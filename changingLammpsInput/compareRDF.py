#! /Program Files (x86)/PythonDownload
from readFiles import (readMyRdfFile, readMyRdfFileOneAtaTime)
import matplotlib.pyplot as plt

def collectiveRDF():

    # rdf files to read

    fileBeginning = 'rdf_D_'
    fileEnd = '.rdf'
    numList = ['0.50','0.25','0.20','0.125']
    colors = ['b.','g^','mx','rs']
    for indx in range(len(numList)):
        fullFile = fileBeginning + str(numList[indx]) + fileEnd
        tmpR, tmpGr = readMyRdfFile(fullFile,3,1,200,100)
        plt.plot(tmpR,tmpGr,colors[indx],ms = 3,label = 'ρ = ' + numList[indx])
        plt.plot(tmpR,tmpGr,colors[indx][0],lw = 1,label = None)
    
    oneLine = [1.0 for i in tmpR]
    plt.plot(tmpR,oneLine,'k',ls = '--',label = None)
    plt.xlim(0,5)
    plt.ylim(0,3)
    plt.xlabel("Distance (d)")
    plt.ylabel("g(r)")
    plt.legend()
    #plt.legend(['T = 1.0','T = 3.0','T = 5.0','T = 7.0','T = 9.0'])
    plt.show()

def singleRDF():
    
    rdfFile = 'rdf_D_0.25.rdf'
    rList, tsRDFLists = readMyRdfFileOneAtaTime(rdfFile,3,1,200)
    colors = ['b','g','r','c','m','y','k']
    for indx in range(500,507):
        plt.plot(rList,tsRDFLists[indx],colors[indx%7],label = str(indx*100))
    
    oneLine = [1.0 for i in rList]
    plt.plot(rList,oneLine,'k',ls = '--',label = None)
    plt.legend()
    plt.xlabel("Distance (d)")
    plt.ylabel('g(r)')
    plt.show()


def main():
    #singleRDF()
    collectiveRDF()

main()
