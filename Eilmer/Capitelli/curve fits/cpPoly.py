import CEAData
import outputDat
import plotting

import csv
from scipy import optimize, integrate
import matplotlib.pyplot as plt
import numpy as np
import math

class CoefficientList:
    def __init__(self, species):
        self.name = species
        self.cpList = {}
        self.enthalpyList = {}
        self.entropyList = {}
        self.luaString = 'db[\'' + self.name + '\'].thermoCoeffs = {\n'

    def setCp(self, tempRange, coefficients):
        self.cpList[tempRange] = coefficients

    def setEnthalpy(self, tempRange, coefficients):
        self.enthalpyList[tempRange] = coefficients

    def setEntropy(self, tempRange, coefficients):
        self.entropyList[tempRange] = coefficients

CpTemps = []
enthalpyTemps = []
entropyTemps = []

speciesList = ['O', 'O+', 'O2', 'O2+', 'N', 'N+', 'N2', 'N2+', 'NO', 'NO+']
heatOfFormation = [249.175, 1568.787, 0.00, 1171.828, 462.68, 1882.128, 0.00, 1509.508, 91.271, 990.810]

CpList = {species: [] for i,species in enumerate(speciesList)}
enthalpyList = {species: [] for i,species in enumerate(speciesList)}
entropyList = {species: [] for i,species in enumerate(speciesList)}

R = 8.3144598

with open('thermoData.csv', 'r') as readFile:
    readData = csv.reader(readFile, delimiter=',')
    next(readData)
    next(readData)
    for row in readData:
        CpTemps.append(float(row[0]))
        for i, species in enumerate(speciesList):
            CpList[species].append(float(row[i + 1]) / R)

with open('enthalpyData.csv', 'r') as readFile:
    readData = csv.reader(readFile, delimiter=',')
    next(readData)
    next(readData)
    for row in readData:
        enthalpyTemps.append(float(row[0]))
        for i, species in enumerate(speciesList):
            if row[i + 1] != '':
                enthalpy = (float(row[i + 1]) - heatOfFormation[i])/ (R * float(row[0]))
                enthalpyList[species].append(enthalpy)

with open('entropyData.csv', 'r') as readFile:
    readData = csv.reader(readFile, delimiter=',')
    next(readData)
    next(readData)
    for row in readData:
        entropyTemps.append(float(row[0]))
        for i, species in enumerate(speciesList):
            if row[i + 1] != '':
                entropyList[species].append(float(row[i + 1]) / R)

plotSpecies = ''
regions = {'O':   [0, 6000, 14000, 22000, 50000],   #good
           'O+':  [0, 14000, 32000, 50000],         #good
           'N':   [0, 6000, 12000, 22000, 50000],   #good
           'N+':  [0, 6000, 24000, 38000, 50000],   #good
           'O2':  [0, 600, 4000, 10000, 50000],     #good
           'O2+': [0, 600, 4000, 14000, 50000],     #good
           'N2':  [0, 600, 10000, 20000, 50000],    #good
           'N2+': [0, 600, 5000, 12000, 50000],     #good
           'NO':  [0, 600, 4000, 15000, 50000],     #good
           'NO+': [0, 600, 12000, 20000, 50000],    #good
           'e':   [0, 50000]}                       #good
           

def CpFunc(x, a, b, c, d, e, f, g):
    y = a*x**-2. + b*x**-1 + c + d*x + e*x**2. + f*x**3. + g*x**4.
    return y

def makeEnthalpy(a, b, c, d, e, f, g):
    def enthalpy(x, intCoeff):
        y = -a*x**-2. + b*np.log(x)/x + c + d*x/2. + e*x**2./3. + f*x**3./4. + g*x**4./5. + intCoeff/x
        return y/x
    return enthalpy

def enthFunc(x, a, b, c, d, e, f, g, intCoeff):
    y = -a*x**-2. + b*np.log(x)/x + c + d*x/2. + e*x**2./3. + f*x**3./4. + g*x**4./5. + intCoeff/x
    return y/x

def makeEntropy(a, b, c, d, e, f, g):
    def entropy(x, intCoeff):
        y = (-a*x**-2.)/2. - b*x**-1. + c*np.log(x) + d*x + (e*x**2.)/2. + (f*x**3.)/3. + (g*x**4.)/4. + intCoeff
        return y
    return entropy

def entrFunc(x, a, b, c, d, e, f, g, intCoeff):
    y = (-a*x**-2.)/2. - b*x**-1. + c*np.log(x) + d*x + (e*x**2.)/2. + (f*x**3.)/3. + (g*x**4.)/4. + intCoeff
    return y

allSpeciesList = {}
for species in speciesList:
    speciesObj = CoefficientList(species)

    speciesRegions = sorted(regions[species])

    CpSpecies = CpList[species]

    enthalpySpecies = sorted(enthalpyList[species])
    entropySpecies = sorted(entropyList[species])

    for i in range(len(speciesRegions[:-1])):

        tempsMin = min([x for x in CpTemps if x >= speciesRegions[i]])
        tempsMax = max([x for x in CpTemps if x <= speciesRegions[i+1]])

        tempsMinIdx = CpTemps.index(tempsMin)
        tempsMaxIdx = CpTemps.index(tempsMax) + 1

        tempRegion = CpTemps[tempsMinIdx:tempsMaxIdx]
        cpRegion = CpSpecies[tempsMinIdx:tempsMaxIdx]
        

        cpWeight = np.empty(len(cpRegion))
        cpWeight.fill(1000)
        cpWeight[0] = cpWeight[-1] = 1e-5

        p, e = optimize.curve_fit(CpFunc, tempRegion, cpRegion, sigma = cpWeight)
        
        speciesObj.setCp((tempsMin, tempsMax), p)

        ######### Enthalpy ###########
        lenTemps = len(enthalpyTemps[:tempsMaxIdx])

        if i == 0:
            boundaryTempEnth = enthalpyTemps[int(lenTemps/4)]
            boundaryEnth = enthalpySpecies[int(lenTemps/4)]
            Ep, Ee = optimize.curve_fit(makeEnthalpy(*p), boundaryTempEnth, boundaryEnth)
            prevEnth = enthFunc(tempsMax, *p, *Ep)
        else:
            
            boundaryEnth = prevEnth
            Ep, Ee = optimize.curve_fit(makeEnthalpy(*p), tempsMin, boundaryEnth)
            prevEnth = enthFunc(tempsMax, *p, *Ep)
        speciesObj.setEnthalpy((tempsMin, tempsMax), Ep)
        
        ######### Entropy ###########
        
        if i == 0:
            boundaryTempEntr = entropyTemps[int(lenTemps/3)]
            boundaryEntr = entropySpecies[int(lenTemps/3)]
            Sp, Se = optimize.curve_fit(makeEntropy(*p), boundaryTempEntr, boundaryEntr)
            prevEntr = entrFunc(tempsMax, *p, *Sp)
            
        else:
            boundaryEntr = prevEntr
            Sp, Se = optimize.curve_fit(makeEntropy(*p), tempsMin, boundaryEntr)
            prevEntr = entrFunc(tempsMax, *p, *Sp)
        
        speciesObj.setEntropy((tempsMin, tempsMax), Sp)
        
    allSpeciesList[species] = speciesObj


outputDat.outputDat(allSpeciesList)

CEA = CEAData.getCEA()
CEA = None

'''
for species in speciesList:
    plt.figure()
    plt.title(species + ' Cp')
    Capitelli = [CpTemps, CpList[species]]
    plotTemps = range(200, 50000, 100)
    plotting.plotCp(species, plotTemps, allSpeciesList, CEA, Capitelli)
    plt.ylim([min(Capitelli[1])*0.8, max(Capitelli[1])*1.2])

    plt.figure()
    plt.title(species + ' Enthalpy')
    Capitelli = [enthalpyTemps, enthalpyList[species]]
    plotTemps = range(200, 50000, 100)
    plotting.plotEnthalpy(species, plotTemps, allSpeciesList, CEA, Capitelli)

    plt.figure()
    plt.title(species + ' Entropy')
    Capitelli = [entropyTemps, entropyList[species]]
    plotTemps = range(200, 50000, 100)
    plotting.plotEntropy(species, plotTemps, allSpeciesList, CEA, Capitelli)

    plt.show()
'''

