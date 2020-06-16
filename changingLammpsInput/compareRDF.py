#! /Program Files (x86)/PythonDownload
from readFiles import readMyRdfFile
import matplotlib.pyplot as plt

def main():

    # rdf files to read

    numList = range(1,10,2)
    fileBeginning = 'rdf_T_'
    fileEnd = '.rdf'
    for element in numList:
        fullFile = fileBeginning + str(element) + fileEnd
        tmpR, tmpGr = readMyRdfFile(fullFile,3,1,100,1000000/100)
        plt.plot(tmpR,tmpGr)

    plt.xlim(0,4)
    plt.ylim(0,2.5)
    plt.xlabel("Distance (d)")
    plt.ylabel("g(r)")
    plt.legend(['T = 1.0','T = 3.0','T = 5.0','T = 7.0','T = 9.0'])
    plt.show()

main()
