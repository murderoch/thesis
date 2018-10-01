#import operator
from math import sqrt

import util
import species

#constants = util.Constants()

class Configuration:
    def __init__(self, core, excitedState, term):
        self.core = core
        self.excitedState = excitedState
        self.term = term

        if self.excitedState:
            self.ID = self.core.atomicNotation + '.' + self.excitedState.subShell.getSubShellString() + '--' + term.getTermString()
        else:
            self.ID = self.core.atomicNotation + '--' + term.getTermString()

class Core:
    def __init__(self, configuration):
        self.configuration = configuration
        self.noElectrons = 0    
        self.atomicNotation = ''
        self.maxN = 0
        self.maxL = 0

        for subShell in configuration:
            self.noElectrons += subShell.noElectrons
            self.atomicNotation = self.atomicNotation + '.' + subShell.getSubShellString()
            if subShell.n > self.maxN:
                self.maxN = subShell.n
            if subShell.L > self.maxL:
                self.maxL = subShell.L 

        self.atomicNotation = self.atomicNotation[1:]

        self.maxS = self.getSMax()
    
    def getSMax(self):
        Ms = 0
        for subShell in self.configuration:
            n = subShell.n
            L = subShell.L
            noElectrons = subShell.noElectrons
            maxElec = 4 * L + 2

            if noElectrons < maxElec * 0.5:
                Ms += noElectrons * 0.5
            else:
                Ms += (maxElec - noElectrons) * 0.5
        return Ms
            
class SubShell:
    def __init__(self, n, L, noElectrons):
        self.n = n
        self.L = L
        self.noElectrons = noElectrons

    def getSubShell(self):
        return [self.n, self.L, self.noElectrons]

    def getSubShellString(self):
        return str(self.n) + str(util.subShellMap[self.L]) + str(self.noElectrons)

class ExcitedState:
    def __init__(self, subShell):
        self.subShell = subShell
        self.S = 0.5
        self.L = self.subShell.L


class Term:
    def __init__(self, L, S, n):
        self.n = n
        self.L = L
        self.S = S
        self.J = self.getJ()
    
        self.levels = []
        for J in self.J:
            self.levels.append(Level(J, None))

    def getJ(self):
        JMax = self.L + self.S
        JMin = abs(self.L - self.S)

        i = JMax
        J = []
        while i >= JMin:
            J.append(i)
            i = i - 1
        return J

    def getTerm(self):
        return [self.L, self.S, self.J]

    def getTermString(self):
        return str(2.*self.S + 1.) + util.termMap[self.L] + str(self.J)


class Level:
    def __init__(self, J, energy):
        self.J = J
        self.energy = energy

    def setEnergy(self, energy):
        self.energy = energy


class CalcEnergy:
    def __init__(self, nMax, species, NIST, theory, allSpecies):
        self.nMax = nMax
        self.species = species
        self.NIST = NIST
        self.theory = theory
        self.allSpecies = allSpecies
        self.constants = util.Constants()

        self.ions = []
        for ionSpecies in self.allSpecies:
            speciesObj = ionSpecies[0]
            if speciesObj.noElectrons == self.species.noElectrons and speciesObj.name != self.species.name:
                self.ions.append(ionSpecies)

    def populateTheory(self):        
        allLevels = self.theory

        for config in allLevels:
            #print(config.ID)
            NISTConfig = next((x for x in self.NIST if x.ID == config.ID), None)
            #print('NIST', NISTConfig)
            if NISTConfig:
                #print(NISTConfig.ID)
                #print(NISTConfig.term.levels)
                for level in config.term.levels:
                    #print(level.J)
                    NISTlevel = next((x for x in NISTConfig.term.levels if x.J == level.J), None)
                    if NISTlevel:
                        level.setEnergy(NISTlevel.energy)
                    else:
                        #print(config.ID)
                        energy = self.getEnergy(config, level)
                        level.setEnergy(energy)
            else:
                for level in config.term.levels:
                    #print(config.ID)
                    energy = self.getEnergy(config, level)
                    level.setEnergy(energy)

        return allLevels

    def getEnergy(self, configuration, level):      
        #self.energy = self.RitzRydberg(self.configuration, self.level)
        
        ionNISTLevels = []
        if configuration.term.n < 20:
            for ion in self.ions:
                speciesObj = ion[0]
                speciesNIST = ion[1]
                speciesNISTConfig = next((x for x in speciesNIST if x.ID == configuration.ID), None)

                if speciesNISTConfig:
                    #print(speciesNISTConfig.ID)
                    speciesNISTLevel = next((x for x in speciesNISTConfig.term.levels if x.J == level.J), None)
                    if speciesNISTLevel:
                        ionNISTLevels.append([speciesObj, speciesNISTLevel, speciesNISTConfig])

        if len(ionNISTLevels) == 0:
            energy = self.hydrogenicApprox(configuration.term.n)
        else:
            energy = self.ritzRydberg(configuration, level, ionNISTLevels, configuration.term.n)

        return energy

    def hydrogenicApprox(self, n):
        IH = self.constants.IH
        Io = self.species.Io

        energy = Io - IH/n**2
        return energy


    def ritzRydberg(self, configuration, level, ionNISTLevels, n):
        IH = self.constants.IH
        I = configuration.core.energy
        z = self.species.charge
        print(I)
        A, B = self.getCoefficients(ionNISTLevels)
        energy = I - IH*(z + 1.)**2. / (n + A + B/n**2.)**2.
        #print(energy)

        return energy

    def getCoefficients(self, ionNISTLevels):
        if len(ionNISTLevels) == 1:
            speciesObj = ionNISTLevels[0][0]
            speciesNISTLevel = ionNISTLevels[0][1]
            speciesNISTConfig = ionNISTLevels[0][2]

            print(self.constants.IH, speciesObj.Io)

            B = 0
            A = sqrt(self.constants.IH * (speciesObj.charge + 1.)**2. / (speciesObj.Io - speciesNISTLevel.energy)) - speciesNISTConfig.term.n
        else:
            A = 0
            B = 0
        
        print(A, B)
        return A, B



def clebschGordon(x, y):
    i = x + y
    iMin = abs(x - y)
    QNumberArray = []
    while i >= iMin:
        QNumberArray.append(i)
        i = i - 1
    return QNumberArray

def readNISTSpectra(species):
    NIST = []
    filename = util.getSpectralDataDir() + species.name + '.dat'
    
    def isNumber(string):
        try:
            int(string)
            return True
        except ValueError:
            return False


    def formatConfig(string):
        nStr = ''
        LStr = ''
        noElectronsStr = ''
        LIdx = len(string) + 1
        for i, char in enumerate(string):
            if isNumber(char) and i < LIdx:
                nStr += char
            elif not isNumber(char):
                LStr += char
                LIdx = i
            if isNumber(char) and i >= LIdx:
                noElectronsStr += char

        n = int(nStr)
        L = util.subShellMap.index(LStr)
        if noElectronsStr != '':
            noElectrons = int(noElectronsStr)
        else:
            noElectrons = 1
        return n, L, noElectrons


    def readConfiguration(inputStr):
        configs = inputStr.split('.')
        coreShells = [SubShell(1, 0, 2)]
        for config in configs[:-1]: 
            n, L, noElectrons = formatConfig(config)
            coreShells.append(SubShell(n, L, noElectrons))

        core = Core(coreShells)
        n, L, noElectrons = formatConfig(configs[-1])
        excitedShell = SubShell(n, L, noElectrons)
        excitedState = ExcitedState(excitedShell)

        return core, excitedState

    def readTerm(inputStr):
        LStr = ''
        SStr = ''
        LIdx = len(inputStr) + 1
        for i, char in enumerate(inputStr):
            if isNumber(char) and i < LIdx:
                SStr += char
            elif not isNumber(char):
                LStr += char
                LIdx = i
        L = util.termMap.index(LStr)
        S = (float(SStr) - 1.)/2.
        return L, S

    with open(filename) as dataFile:
        next(dataFile)
        IDs = []
        for row in dataFile:
            if not row.strip():
                continue

            values = row.rstrip('\n').split(', ')
            core, excitedState = readConfiguration(values[0])
            n = excitedState.subShell.n

            L, S = readTerm(values[1])
            term = Term(L, S, n)
            configuration = Configuration(core, excitedState, term)
            if configuration.ID not in IDs:
                IDs.append(configuration.ID)
                NIST.append(configuration)
            else:
                configuration = NIST[IDs.index(configuration.ID)]
            
            for level in configuration.term.levels:
                if configuration.term.L == L and configuration.term.S == S and level.J == float(values[2]):
                    level.setEnergy(float(values[3]))
    return NIST


def calculateExpectedStates(species, nMax):
    states = []

    cores = []
    for core in species.cores:
        coreShells = []
        for subShell in core[0]:
            coreShells.append(SubShell(subShell[0], subShell[1], subShell[2]))
        cores.append(Core(coreShells))
    for core in cores:

        if core.noElectrons == species.noElectrons:
            n = core.maxN
            LArr = range(0, core.maxL+1, 1)
            SArr = []

            i = core.maxS

            while i >= 0:
                SArr.append(i)
                i = i - 1
        
            for L in LArr:
                for S in SArr:
                    term = Term(L, S, n)
                    configuration = Configuration(core, None, term)
                    states.append(configuration)
        
        else:
            for n in range(core.maxN+1, nMax+1, 1):
                for L in range(0, n, 1):
                    excitedShell = SubShell(n, L, 1)
                    excitedState = ExcitedState(excitedShell)
                    SArr = clebschGordon(core.maxS, excitedState.S)
                    n = excitedState.subShell.n

                    for S in SArr:
                        term = Term(L, S, n)
                        configuration = Configuration(core, excitedState, term)
                        states.append(configuration)
                        
    return states