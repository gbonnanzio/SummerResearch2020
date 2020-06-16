#! Program Files (x86)/PythonDownload
import matplotlib.pyplot as plt
from readFiles import readMyPeOUTFile

def main():
	myFile = 'potentialEnergies.out'
	ts, pe = readMyPeOUTFile(myFile)
	plt.plot(ts,pe)
	plt.xlabel('Time step')
	plt.ylabel('Potential Energy (kT)')
	plt.legend(['Potential Energy']) 
	plt.show()

main()	
