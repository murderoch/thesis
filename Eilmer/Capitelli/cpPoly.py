import csv
from scipy import optimize, interpolate
import matplotlib.pyplot as plt
import numpy as np

class CoefficientList:
    def __init__(self, species):
        self.name = species
        self.tempList = {}
        self.luaString = 'db[\'' + self.name + '\'].thermoCoeffs = {\n'

    def setCoefficients(self, tempRange, coefficients):
        self.tempList[tempRange] = coefficients


temps = []
speciesList = ['O', 'O+', 'O2', 'O2+', 'N', 'N+', 'N2', 'N2+', 'NO', 'NO+', 'e']

#speciesList = ['O+', 'O2', 'O2+', 'N', 'N+', 'N2', 'N2+', 'NO', 'e']

CpList = {species: [] for i,species in enumerate(speciesList)}

R = 8.3144598

with open('thermoData.csv', 'r') as readFile:
    readData = csv.reader(readFile, delimiter=',')
    next(readData)
    next(readData)
    for row in readData:
        temps.append(float(row[0]))
        for i, species in enumerate(speciesList):
            CpList[species].append(float(row[i + 1]) / R)


plotSpecies = ''
regions = {'O':   [0, 6000, 14000, 22000, 50000],   #good
           'O+':  [0, 14000, 32000, 50000],         #good
           'N':   [0, 6000, 14000, 22000, 50000],   #good
           'N+':  [0, 14000, 30000, 50000],         #good
           'O2':  [0, 600, 4000, 10000, 50000],     #good
           'O2+': [0, 600, 4000, 14000, 50000],     #good
           'N2':  [0, 600, 10000, 20000, 50000],    #good
           'N2+': [0, 600, 5000, 12000, 50000],     #good
           'NO':  [0, 600, 4000, 15000, 50000],     #good
           'NO+': [0, 600, 12000, 20000, 50000],    #good
           'e':   [0, 50000]}


def func(x, a, b, c, d, e, f, g):
    y = a*x**-2. + b*x**-1 + c + d*x + e*x**2. + f*x**3. + g*x**4. + g*x**5. + g*x**6.
    return y

allSpeciesList = []

with open('output.dat', 'w') as outputFile:

    for species in speciesList:

        speciesObj = CoefficientList(species)

        outputFile.writelines('\n\n~~~~~~~' + species + '~~~~~~~\n')
        
        if species == plotSpecies or plotSpecies == 'all':
            plt.figure()
        speciesRegions = regions[species]
        CpSpecies = CpList[species]
        for i in range(len(speciesRegions[:-1])):

            tempsMin = min([x for x in temps if x >= speciesRegions[i]])

            tempsMax = max([x for x in temps if x <= speciesRegions[i+1]])

            tempsMinIdx = temps.index(tempsMin)
            tempsMaxIdx = temps.index(tempsMax) + 1

            tempRegion = temps[tempsMinIdx:tempsMaxIdx]
            cpRegion = CpSpecies[tempsMinIdx:tempsMaxIdx]

            p, e = optimize.curve_fit(func, tempRegion, cpRegion)
            
            xd = np.linspace(min(tempRegion), max(tempRegion), 1000)

            if species == plotSpecies or plotSpecies == 'all':
                plt.plot(xd, func(xd, *p))
                plt.title(species)
                plt.plot(temps, CpSpecies, 'kx')

            coeffStr = str(tempsMin) + '-' + str(tempsMax) 
            spaces = 16 - len(coeffStr)
            for i in range(spaces):
                coeffStr += ' '
            coeffStr += '| \n'
            outputFile.writelines(coeffStr)

            coeffStr = ''
            for coeff in p:
                if coeff > 0:
                    coeffStr += "       {:.9e}".format(coeff) + ',\n'
                else:
                    coeffStr += "      {:.9e}".format(coeff) + ',\n'
            outputFile.writelines(coeffStr)

            speciesObj.setCoefficients((tempsMin, tempsMax), p)
        
        allSpeciesList.append(speciesObj)

#plt.show()



def getCoeffs(temp, species):
    for key, value in species.tempList.items():
        if temp < key[1] and temp >= key[0]:
            return value
    return 0


plotTemps = range(300, 50000, 100)
yplot = []

species = allSpeciesList[1]
for temp in plotTemps:
    coeff = getCoeffs(temp, species)
    print(temp, coeff)
    yplot.append(func(temp, coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5], coeff[6]))

plt.figure()
plt.plot(plotTemps, yplot)
plt.show()


