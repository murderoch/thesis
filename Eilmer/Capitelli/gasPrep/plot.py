import matplotlib.pyplot as plt 


TList = []
CpList = []

with open('dataFile.dat', 'r') as inFile:
    next(inFile)
    for row in inFile:
        data = row.split(' ')
        TList.append(float(data[1]))
        CpList.append(float(data[2]))

plt.plot(TList, CpList)
plt.show()

