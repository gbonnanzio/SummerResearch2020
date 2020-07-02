import glob
def main():
    myList = glob.glob('./runLammps/dump.coords*')
    print(myList[0])

main()


