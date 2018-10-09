import matplotlib.pyplot as plt

cells = [30, 60, 120]

plotCells = 120
species = ['N', 'O']
rR = '2_5e-5'

folderName = 'nonReac/'

CapFileName = species[0] + '_' + species[1] + '-Cap.dat'
CEAFileName = species[0] + '_' +  species[1] + '-CEA.dat'

#capFiles = ['Cap-pR2_5e-5.dat', 'Cap-pR1_0e-4.dat', 'Cap-pR5_0e-6.dat']
#CEAFiles = 


CEA = {}
with open(folderName + CEAFileName, 'r') as CEAFile:
    for row in CEAFile:
        data = row.split(' ')
        ID = (data[0], int(data[2]))
        CEA[ID] = float(data[8])

Cap = {}
#for fileName in capFiles:
with open(folderName + CapFileName, 'r') as CapFile:
    for row in CapFile:
        data = row.split(' ')
        ID = (data[0], int(data[2]))
        Cap[ID] = float(data[8])


plt.figure(2)
plt.title('pR5.0e-6')
plt.figure(3)
plt.title('pR1.0e-4')
plt.figure(1)
plt.xlabel('Velocity (km/s)')
plt.ylabel('Proportional Shock Standoff (D/d)')
plt.title('Lobb-Sphere Shock Standoff Comparison for O+N Monatoms')


CEAPoints = [[],[]]
CapPoints = [[],[]]

for keys, results in CEA.items():
    uIdx = keys[0].index('u')
    vel = float(keys[0][uIdx+1:])
    print(vel, results)
    if keys[0][:8] == 'pR2.5e-5':
        plt.figure(1)
    elif keys[0][:8] == 'pR5.0e-6':
        plt.figure(2)
    else:
        plt.figure(3)

    if keys[1] == plotCells or plotCells == None:
        CEAPoints[0].append(vel)
        CEAPoints[1].append(results)
        #plt.plot(vel, results, 'ko', markersize=7, label='CEA Coefficients')
    elif keys[1] == plotCells or plotCells == None:
        plt.plot(vel, results, 'g+')
    elif keys[1] == plotCells or plotCells == None:
        plt.plot(vel, results, 'r+')

for keys, results in Cap.items():
    uIdx = keys[0].index('u')
    vel = float(keys[0][uIdx+1:])
    if keys[0][:8] == 'pR2.5e-5':
        plt.figure(1)
    elif keys[0][:8] == 'pR5.0e-6':
        plt.figure(2)
    else:
        plt.figure(3)

    if keys[1] == plotCells or plotCells == None:
        CapPoints[0].append(vel)
        CapPoints[1].append(results)
        #plt.plot(vel, results, 'r+', markersize = 10, label='Updated Coefficients')
    elif keys[1] == plotCells or plotCells == None:
        plt.plot(vel, results, 'go')
    elif keys[1] == plotCells or plotCells == None:
        plt.plot(vel, results, 'ro')

Lobb = [[4.600, 4.900, 5.100, 5.450, 5.700, 6.400], [0.048, 0.055, 0.045, 0.05, 0.052, 0.049]]
Zandar = [[8.66, 9.7], [0.036, 0.0335]]

plt.figure(1)
plt.plot(CEAPoints[0], CEAPoints[1], 'ko', markersize=7, label='CEA Coefficients')
plt.plot(CapPoints[0], CapPoints[1], 'r+', markersize=10, label='Updated Coefficients')
#plt.plot(Lobb[0], Lobb[1], 'b^', markersize=7, label='Lobb Experimental Data')
#plt.plot(Zandar[0], Zandar[1], 'bs', markersize=7, label='Zandar Experimental Data (Mean)')
plt.ylim([0.0, 0.13])
plt.legend(loc = 4)
plt.grid()
plt.show()