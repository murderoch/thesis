from math import sqrt
from itertools import repeat, product, combinations
from collections import Counter
import numpy as np
from scipy import optimize

import matplotlib.pyplot as plt

import util
import species

#constants = util.Constants()

class Configuration:
    def __init__(self, core, excitedState, term):
        self.core = core
        self.excitedState = excitedState
        self.term = term

        if core.term:
            coreTermStr = '(' + str(core.term.S * 2 + 1.) + str(util.termMap[core.term.L]) + ')'
        elif self.excitedState:
            coreTermStr = ''


        if self.excitedState:
            self.ID = self.core.atomicNotation + coreTermStr + '.' + self.excitedState.subShell.getSubShellString() + '--' + term.getTermString()
        else:
            self.ID = self.core.atomicNotation + '--' + term.getTermString()

        self.method = None

    def setMethod(self, method):
        self.method = method
        
class Core:
    def __init__(self, subShells, term):
        self.subShells = subShells
        self.term = term

        self.noElectrons = 0    
        self.atomicNotation = ''
        self.maxN = 0
        self.maxL = 0

        for subShell in subShells:
            self.noElectrons += subShell.noElectrons
            self.atomicNotation = self.atomicNotation + '.' + subShell.getSubShellString()
            if subShell.n > self.maxN:
                self.maxN = subShell.n
            if subShell.L > self.maxL:
                self.maxL = subShell.L 

        self.atomicNotation = self.atomicNotation[1:]
        self.maxS = self.getSMax()

        if term:
            self.ID = self.getCoreID()

    def getCoreID(self):
        ID = self.atomicNotation + '--' + self.term.getTermString()
        return ID

    def getSMax(self):
        Ms = 0
        for subShell in self.subShells:
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

    def checkLevelEnergies(self):
        self.maxEnergy = 0
        for level in self.levels:
            if level.energy == None:
                level.energy = 0
            if level.energy > self.maxEnergy:
                self.maxEnergy = level.energy

        emptyLevels = []
        for level in self.levels:
            if level.energy == None:
                emptyLevels.append(level)
        if len(emptyLevels) == 0:
            return None
        else:
            return emptyLevels

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
        self.constants = util.Constants()

        self.speciesIon = None

        with open('hi.dat', 'w') as outFile:
            for config in self.NIST:
                outFile.write(config.ID + '\n')

        self.comprableIons = []
        for ionSpecies in allSpecies:
            speciesObj = ionSpecies[0]
            if speciesObj.noElectrons == self.species.noElectrons and speciesObj.name != self.species.name:
                self.comprableIons.append(ionSpecies)

            elif speciesObj.atomicNumber == self.species.atomicNumber and speciesObj.noElectrons == self.species.noElectrons - 1:
                self.speciesIon = ionSpecies

        if len(self.comprableIons) == 0:
            print('missing ion data for ', self.species.name)


    def populateTheory(self):
        allLevels = self.theory
        calculatedLevels = []
        for config in allLevels:
            NISTConfig = next((x for x in self.NIST if x.ID == config.ID), None)
            if NISTConfig:
                for level in config.term.levels:
                    NISTlevel = next((x for x in NISTConfig.term.levels if x.J == level.J and x.energy != None), None)
                    if NISTlevel:
                        level.setEnergy(NISTlevel.energy)
                    else:
                        energy = self.getEnergy(allLevels, config, level)
                        level.setEnergy(energy)
            
            else:
                for level in config.term.levels:
                    energy = self.getEnergy(allLevels, config, level)
                    level.setEnergy(energy)
                calculatedLevels.append(config)                   
            

        return allLevels, calculatedLevels

    def getEnergy(self, allLevels, configuration, level):      
        energy = None
        if configuration.term.n <= 20 and self.speciesIon != None:
            
            if configuration.excitedState == None:
                configuration.setMethod('coreIonExtrap')
                energy = self.coreIonEnergyExtrapolation(configuration, level)

            else:
                matchedCores = []
                for config in self.NIST:
                    if config.excitedState:
                        print(config.ID)
                        if config.core.ID == configuration.core.ID:
                            matchedCores.append(config)
                            
                if len(matchedCores) > 0:
                    matchedTerms = []
                    for config in matchedCores:
                        if config.term.getTerm() == configuration.term.getTerm():
                            matchedTerms.append(config)

                    if len(matchedTerms) > 0:
                        coreEnergy = None
                        speciesIonNISTConfig = next((x for x in self.speciesIon[1] if x.ID == configuration.core.ID))
                        if speciesIonNISTConfig:
                            energyList = [x.energy for x in speciesIonNISTConfig.term.levels]
                            coreEnergy = sum(energyList)/len(energyList)

                        if coreEnergy != None:
                            configuration.setMethod('ritzRydberg')
                            energy = self.ritzRydberg(configuration, level, coreEnergy, matchedTerms)
                            
                    else:
                        configuration.setMethod('azimuthalQuantExtrap')
                        energy = self.azimuthalQuantumExtrapolation(configuration, level, matchedCores)
                        
        if energy == None:
            configuration.setMethod('hydrogenic')
            energy = self.hydrogenicApprox(configuration.term.n)
            

        return energy

    def hydrogenicApprox(self, n):
        IH = self.constants.IH
        Io = self.species.Io
        energy = Io - IH/n**2
        return energy

    def coreIonEnergyExtrapolation(self, configuration, level):
        ionNISTLevels = []
        for ion in self.comprableIons:
            speciesObj = ion[0]
            speciesNIST = ion[1]
            speciesNISTConfig = next((x for x in speciesNIST if x.ID == configuration.ID), None)
            if speciesNISTConfig:
                speciesNISTLevel = next((x for x in speciesNISTConfig.term.levels if x.J == level.J), None)
                if speciesNISTLevel:
                    ionNISTLevels.append([speciesObj.charge, speciesNISTLevel.energy])
        
        if len(ionNISTLevels) > 1:
            # fit linear trend to comparable ions (same electronic config)
            m = (ionNISTLevels[0][1] - ionNISTLevels[1][1]) / (ionNISTLevels[0][0] - ionNISTLevels[1][0])
            c = ionNISTLevels[0][1] - m * ionNISTLevels[0][0]
            
            # extrapolate to species' charge
            energy = m * self.species.charge + c         
            return energy
        else:
            return None

    def azimuthalQuantumExtrapolation(self, configuration, level, matchedCores):       
        fitEnergyList = []
        fitLList = []
        for config in matchedCores:
            if config.excitedState.subShell.n == configuration.excitedState.subShell.n:
                
                energyList = [x.energy for x in config.term.levels if x.energy != None]
                avgEnergy = sum(energyList)/len(energyList)
                fitLList.append(config.excitedState.L)
                fitEnergyList.append(avgEnergy)

        if len(fitEnergyList) == 0 or len(set(fitLList)) <= 1:
            energy = self.hydrogenicApprox(configuration.term.n)
            configuration.setMethod('hydrogenic')
        else:
            fit = np.polyfit(fitLList, fitEnergyList, 1)
            energy = fit[0] * configuration.excitedState.L + fit[1]
        
        return energy


    def ritzRydberg(self, configuration, level, coreEnergy, matchedTerms):
        IH = self.constants.IH
        z = self.species.atomicNumber
        n = configuration.term.n
        I = self.species.Io + coreEnergy

        A, B = self.getCoefficients(coreEnergy, matchedTerms)

        if A != None and B != None:
            energy = I - IH*(z + 1.)**2. / (n + A + B/n**2.)**2.
        else:
            energy = None


        if 0:
            def func(n, A, B):
                energy = I - IH*(z + 1.)**2. / (n + A + B/n**2.)**2.
                return energy

            if configuration.ID == '1s2.2s2.2p3(4.0S).19p1--3.0P[2.0, 1.0, 0.0]':
                plt.figure()
                xPlot = range(2,10)
                yPlot = [func(x, A, B) for x in xPlot]
                print('\n', configuration.ID, 'No Matches:', len(matchedTerms))
                for config in matchedTerms:
                    energyList = [x.energy for x in config.term.levels]
                    Anergy = sum(energyList)/len(energyList)
                    plt.plot(config.term.n, Anergy, 'bo')
                    print('matched to: ', config.ID)
                    print(config.term.n, Anergy)    
                plt.plot(xPlot, yPlot, 'b--')
                plt.plot(n, energy, 'rx')
                plt.show()
        return energy


    def getCoefficients(self, coreEnergy, matchedTerms):
        # fit curve and calculate coefficients
        IH = self.constants.IH
        z = self.species.atomicNumber
        I = self.species.Io + coreEnergy
        #Ry = self.constants.Rydberg * (self.species.atomicNumber - self.speciesIon[0].noElectrons)**2.

        if len(matchedTerms) == 1:
            config = matchedTerms[0]
            n = config.term.n
            energyList = [x.energy for x in config.term.levels if x.energy != None]
            energy = sum(energyList)/len(energyList)

            if energy > I:
                return None, None

            B = 0
            A = sqrt(self.constants.IH * (z + 1.)**2. / (I - energy)) - n
            #A = sqrt(Ry / (I - energy)) - n
        
        elif len(matchedTerms) >= 2:
            fitnList = []
            fitEnergyList = []
            for config in matchedTerms:
                n = config.term.n
                energyList = [x.energy for x in config.term.levels if x.energy != None]
                energy = sum(energyList)/len(energyList)
                fitnList.append(n)
                fitEnergyList.append(energy) 
            if energy > I:
                return None, None

            def func(n, A, B):
                energy = I - IH*(z + 1.)**2. / (n + A + B/n**2.)**2.
                #energy = I - Ry /  (n + A + B/n**2.)**2.
                return energy

            #sigma = [1 for x in range(len(fitnList))]
            #sigma[-1] = 0.2
            p, _ = optimize.curve_fit(func, fitnList, fitEnergyList)   
            A = p[0]
            B = p[1]
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
        coreJ = None

        _, _, finalSubShellNoElectrons = formatConfig(configs[-1])

        if finalSubShellNoElectrons == 1:
            coreTerm = None
            for config in configs[:-1]: 
                if config[0] == '(':

                    if '<' in config:
                        openIdx = config.index('<')
                        closeIdx = config.index('>')
                        if '/' in config:
                            slashIdx = config.index('/')
                            num = float(config[openIdx+1:slashIdx])
                            denom = float(config[slashIdx+1:closeIdx])
                            coreJ = num/denom
                        else:
                            coreJ = float(config[openIdx+1:closeIdx])

                    coreS = (float(config[1]) - 1.) / 2.
                    coreL = util.termMap.index(config[2])
                    coreTerm = Term(coreL, coreS, coreShells[-1].n)


                else:
                    n, L, noElectrons = formatConfig(config)
                    coreShells.append(SubShell(n, L, noElectrons))

            core = Core(coreShells, coreTerm)
            n, L, noElectrons = formatConfig(configs[-1])
            excitedShell = SubShell(n, L, noElectrons)
            excitedState = ExcitedState(excitedShell)

        else:
            for config in configs:
                n, L, noElectrons = formatConfig(config)
                coreShells.append(SubShell(n, L, noElectrons))
                core = Core(coreShells, None)
                excitedState = None

        return core, excitedState, coreJ

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
            core, excitedState, coreJ = readConfiguration(values[0])
            L, S = readTerm(values[1])
            if excitedState:
                n = excitedState.subShell.n
            else:
                n = core.maxN

            term = Term(L, S, n)
            configuration = Configuration(core, excitedState, term)
            if configuration.ID not in IDs:
                IDs.append(configuration.ID)
                NIST.append(configuration)
            else:
                configuration = NIST[IDs.index(configuration.ID)]
            
            for level in configuration.term.levels:
                if configuration.term.L == L and configuration.term.S == S and level.J == float(values[2]):
                    if values[3] != 'None':
                        level.setEnergy(float(values[3]))
    return NIST

def calculateExpectedStates(species, nMax):
    states = []
    cores = []
    for core in species.cores:
        coreShells = []

        for subShell in core:
            
            n = subShell[0]
            L = subShell[1]
            noElectrons = subShell[2]

            coreShells.append(SubShell(n, L, noElectrons))

            maxElectrons = 2.*(L * 2. + 1.)
            if noElectrons < maxElectrons:
                coreSLs = getTerms(L, noElectrons)

                for coreSL in coreSLs:
                    termObj = Term(coreSL[0], coreSL[1], n)
                    cores.append(Core(coreShells, termObj))
                    
    for core in cores:
        if core.noElectrons == species.noElectrons:
            n = core.maxN

            termObj = core.term
            configuration = Configuration(core, None, termObj)
            states.append(configuration)
        
        else:
            for n in range(core.maxN+1, nMax+1, 1):
                for L in range(0, n, 1):
                    
                    #This gets us to all of the possible excited electron states

                    #now need to consider interaction with coreState

                    #excitedState is good
                    excitedShell = SubShell(n, L, 1)    
                    excitedState = ExcitedState(excitedShell)


                    #Good, all possible S states from combination of exc and core S
                    #added all possible L states from combination of exc and core L
                    SList = clebschGordon(core.term.S, excitedState.S)
                    LList = clebschGordon(core.term.L, excitedState.L)

                    for Lcomb in LList:
                        for Scomb in SList:
                            term = Term(Lcomb, Scomb, n)
                            configuration = Configuration(core, excitedState, term)
                            states.append(configuration)

    return states



def getTerms(L, noElectrons):
    """Return a list of term symbols for the configuration l^noElectrons."""

    # Total number of (ml, ms) pairs for this subshell.
    n = (2*L+1)*2

    # All possible values of ml = -l, -l+1, ..., l-1, l.
    ml = list(range(-L,L+1))

    # All possible values of 2ms = -1, 1. That is, ms = -1/2, +1/2. We work
    # with 2ms instead of ms so that we can handle integers only.
    ms2 = [-1,1]

    # All possible (ml, 2ms) pairs for this subshell.
    ml_ms2 = list(product(ml, ms2))

    # All possible microstates for r electrons in this subshell.
    microstates = list(combinations(range(n), noElectrons))

    # The totals ML = sum(ml) and MS2 = sum(2ms) for each microstate
    ML = [sum([ml_ms2[microstate[j]][0] for j in range(noElectrons)]) for microstate in microstates]
    MS2 = [sum([ml_ms2[microstate[j]][1] for j in range(noElectrons)]) for microstate in microstates]

    # Count the microstates (MS, ML). Store them this way round so we can
    # pick off the ground state term (maximum S) first.
    MS2_ML = Counter(zip(MS2,ML))
    N = len(microstates)

    # Extract the term symbols by starting at the minimum (ML, MS) value and
    # removing microstates corresponding to the (L, S) term it belongs to.
    # Repeat until we're out of microstates.
    terms = []
    while N>0:
        S, L = min(MS2_ML)
        terms.append([-L, -S/2.])
        for ML in range(L, -L+1):
            for MS in range(S, -S+1,2):
                MS2_ML[MS,ML] -= 1
                if MS2_ML[MS,ML] == 0:
                    del MS2_ML[MS,ML]
                N -= 1
    return terms


def sortSpectra(levelList, outputFileName):
    for config in levelList:
        if config.term.checkLevelEnergies() != None:
            print('Missing energy for:', config.ID, config.core.term.getTermString())
            for level in config.term:
                level.setEnergy(0.)

    levelList.sort(key=lambda x: x.term.maxEnergy)
 
    with open(outputFileName + '.dat', 'w') as outputFile:
        for config in levelList:
            outputFile.write('\n')
            for level in config.term.levels:
                writeStr = config.ID + ', ' + str(level.J) + ', ' + str(level.energy)
                
                if outputFileName == 'calculated':
                    writeStr += ' ' + config.method + '\n'
                else:
                    writeStr += '\n'
                
                outputFile.writelines(writeStr)

    return levelList




if __name__ == '__main__':
    O = species.Species('O')

    #allSpecies = [allSpecies]

    NIST = readNISTSpectra(O)
    theory = calculateExpectedStates(O, 5)
    calcEnergy = CalcEnergy(5, O, NIST, theory, O)

    allLevels = calcEnergy.populateTheory()

    for config in allLevels:
        if config.term.checkLevelEnergies() != None:
            print(config.ID, config.core.term.getTermString())