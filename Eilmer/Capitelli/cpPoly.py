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


plotSpecies = 'N2'

'''
regions = {'O':   [0, 6000, 14000, 22000, 50000],   #good
           'O+':  [0, 14000, 32000, 50000],         #good
           'N':   [0, 6000, 14000, 22000, 50000],   #good
           'N+':  [0, 6000, 24000, 38000, 50000],   #good
           'O2':  [0, 600, 4000, 10000, 50000],     #good
           'O2+': [0, 600, 4000, 14000, 50000],     #good
           'N2':  [0, 600, 10000, 20000, 50000],    #good
           'N2+': [0, 600, 5000, 12000, 50000],     #good
           'NO':  [0, 600, 4000, 15000, 50000],     #good
           'NO+': [0, 600, 12000, 20000, 50000],    #good
           'e':   [0, 50000]}
'''
regions = {'O':   [0, 6000, 14000, 22000, 50000],   #good
           'O+':  [0, 14000, 32000, 50000],         #good
           'N':   [0, 6000, 14000, 22000, 50000],   #good
           'N+':  [0, 6000, 24000, 38000, 50000],   #good
           'O2':  [0, 600, 4000, 10000, 50000],     #good
           'O2+': [0, 600, 4000, 14000, 50000],     #good
           'N2':  [100, 600, 10000, 20000, 50000],   #good
           'N2+': [0, 600, 5000, 12000, 50000],     #good
           'NO':  [0, 600, 4000, 15000, 50000],     #good
           'NO+': [0, 600, 12000, 20000, 50000],    #good
           'e':   [0, 50000]}

def func(x, a, b, c, d, e, f, g):
    y = a*x**-2. + b*x**-1 + c + d*x + e*x**2. + f*x**3. + g*x**4.
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

            cpWeight = np.empty(len(cpRegion))
            cpWeight.fill(10)
            cpWeight[0] = cpWeight[-1] = 0.1

            p, e = optimize.curve_fit(func, tempRegion, cpRegion, sigma = cpWeight)
            
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

def getCoeffs(temp, tempList):
    for key, value in tempList:
        if temp < key[1] and temp >= key[0]:
            return value
    return 0

plt.figure()


plotTemps = range(200, 20000, 10)
yplot = []

#species = allSpeciesList[6]
species = { (100, 600):  [ -1.357998069e+03,
                        4.074344664e+01,
                        3.009094920e+00,
                        3.036860470e-03,
                       -9.939853784e-06,
                        1.575896162e-08,
                       -8.493317618e-12,],
    
        (600, 1000): [  3.561032511e+04,
                        -4.823315234e+02,
                        4.176549841e+00,
                        2.797898299e-04,
                        -5.559406916e-08,
                        3.471667200e-12,
                        5.866562848e-17,],
        (1000, 50000): [ 1.000000000e+00,
                        3.948581822e+03,
                        5.389505726e+01,
                        -1.520063528e-02,
                        1.654607432e-06,
                        -7.406707299e-11,
                        1.174825602e-15,]}

for temp in plotTemps:
    #coeff = getCoeffs(temp, species.tempList.items())
    coeff = getCoeffs(temp, species)
    yplot.append(func(temp, coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5], coeff[6]))

plt.plot(plotTemps, yplot, label='mine')

CEA = {(200, 1000): [2.210371497e+04,
                -3.818461820e+02,
                6.082738360e+00,
                -8.530914410e-03,
                1.384646189e-05,
                -9.625793620e-09,
                2.519705809e-12,
                7.108460860e+02,
                -1.076003744e+01],
    
       (1000, 6000): [5.877124060e+05,
                    -2.239249073e+03,
                    6.066949220e+00,
                    -6.139685500e-04,
                    1.491806679e-07,
                    -1.923105485e-11,
                    1.061954386e-15,
                    1.283210415e+04,
                    -1.586640027e+01,],
        (6000, 50000): [8.310139160e+08,
                    -6.420733540e+05,
                    2.020264635e+02,
                    -3.065092046e-02,
                    2.486903333e-06,
                    -9.705954110e-11,
                    1.437538881e-15,
                    4.938707040e+06,
                    -1.672099740e+03,]}


yplot = []
for temp in plotTemps:
    coeff = getCoeffs(temp, CEA.items())
    yplot.append(func(temp, coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5], coeff[6]))

plt.plot(plotTemps, yplot, label = 'CEA')
plt.legend()
plt.show()