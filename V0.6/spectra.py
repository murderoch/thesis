from math import sqrt
from itertools import repeat, product, combinations
from collections import Counter

import util
import species

#constants = util.Constants()

class Configuration:
    def __init__(self, core, excitedState, term):
        self.core = core
        self.excitedState = excitedState
        self.term = term

        '''
        if self.excitedState:
            coreS = self.term.S - 0.5
            coreL = self.term.L - self.excitedState.L
        else:
            coreS = self.term.S
            coreL = self.term.L
        coreN = self.core.maxN

        coreTerms = []
        for subShell in self.core.subShells:
            maxElectrons = 2.*(subShell.L * 2. + 1.)
            if subShell.noElectrons < maxElectrons:
                coreTerms += getTerms(subShell.L, subShell.noElectrons)

        self.coreTerm = Term(coreL, coreS, coreN)
        '''

        if self.excitedState:
            self.ID = self.core.atomicNotation + '.' + self.excitedState.subShell.getSubShellString() + '--' + term.getTermString()
        else:
            self.ID = self.core.atomicNotation + '--' + term.getTermString()

        
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
        
        for ionSpecies in self.allSpecies:
            speciesObj = ionSpecies[0]
            if speciesObj.atomicNumber == self.species.atomicNumber and speciesObj.noElectrons == self.species.noElectrons - 1:
                self.speciesIon = ionSpecies
                break

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
            
            if energy == None:
                energy = self.hydrogenicApprox(configuration.term.n)

        #print(configuration.ID, energy)

        return energy

    def hydrogenicApprox(self, n):
        IH = self.constants.IH
        Io = self.species.Io

        energy = Io - IH/n**2
        return energy


    def ritzRydberg(self, configuration, level, ionNISTLevels, n):
        IH = self.constants.IH
        z = self.species.charge

        #print('\n', configuration.ID)

        def getCoreEnergy(configuration, ionNISTLevels):
            for ionConfig in self.speciesIon[1]:
                #print(ionConfig.ID.split('--')[0], ionConfig.term.getTermString())
                '''
                if ionConfig.ID.split('--')[0] == configuration.core.atomicNotation:
                    print('matched config is', ionConfig.ID, configuration.core.atomicNotation, configuration.core.term.getTermString())
                    print('check term is', ionConfig.term.getTermString(), configuration.core.term.getTermString())
                    if ionConfig.term.getTerm() == configuration.core.term.getTerm():
                        print('matched term is ', ionConfig.ID)
                        for ionLevel in ionConfig.term.levels:
                            print('level is', ionLevel.J, 'want', configuration.core.term.J)
                            if ionLevel.J in configuration.core.term.J:
                                print('matched level is', ionConfig.ID, configuration.core.term.J, 'with energy', ionLevel.energy)
                                #print('matched config is', ionConfig.ID, configuration.core.atomicNotation, configuration.coreTerm.getTermString())
                                
                                coreEnergy = ionLevel.energy
                                return coreEnergy
                '''
            return 0

        I = getCoreEnergy(configuration, ionNISTLevels)

        if I == 0:
            return self.hydrogenicApprox(n)
        else:
            A, B = self.getCoefficients(ionNISTLevels)
            energy = I - IH*(z + 1.)**2. / (n + A + B/n**2.)**2.
            return energy

    def getCoefficients(self, ionNISTLevels):
        if len(ionNISTLevels) == 1:
            speciesObj = ionNISTLevels[0][0]
            speciesNISTLevel = ionNISTLevels[0][1]
            speciesNISTConfig = ionNISTLevels[0][2]

            #print(self.constants.IH, speciesObj.Io)

            B = 0
            A = sqrt(self.constants.IH * (speciesObj.charge + 1.)**2. / (speciesObj.Io - speciesNISTLevel.energy)) - speciesNISTConfig.term.n
        else:
            A = 0
            B = 0
        
        #print(A, B)
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
        #coreTerms = []
        #cores = []
        _, _, finalSubShellNoElectrons = formatConfig(configs[-1])


        if finalSubShellNoElectrons == 1:
            for config in configs[:-1]: 
                n, L, noElectrons = formatConfig(config)
                coreShells.append(SubShell(n, L, noElectrons))
                '''
                maxElectrons = 2.*(L * 2. + 1.)
                if noElectrons < maxElectrons:
                    coreSLs = getTerms(L, noElectrons)
                    for coreSL in coreSLs:
                        termObj = Term(coreSL[0], coreSL[1], n)
                        cores.append(Core(coreShells, termObj))
                '''
            core = Core(coreShells, None)
            n, L, noElectrons = formatConfig(configs[-1])
            excitedShell = SubShell(n, L, noElectrons)
            excitedState = ExcitedState(excitedShell)

        else:
            for config in configs:
                n, L, noElectrons = formatConfig(config)
                coreShells.append(SubShell(n, L, noElectrons))
                core = Core(coreShells, None)
                excitedState = None

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
            L, S = readTerm(values[1])

            if excitedState:
                n = excitedState.subShell.n
                coreS = [S - 0.5, S + 0.5]
                coreL = clebschGordon(L, excitedState.L)

                print(core.atomicNotation, excitedState.subShell.getSubShellString(), coreS, coreL)

                subShellTerms = []
                for subShell in core.subShells:
                    if subShell.noElectrons != (4. * subShell.L + 2.):
                        subShellTerms.append(getTerms(subShell.L, subShell.noElectrons))
                
                print(subShellTerms)
                coreShellTerms = subShellTerms[0]   #Will need to change this to combine multiple shells
                for L in coreL:
                    for S in coreS:
                        print([L, S])
                        if [L, S] in coreShellTerms:
                            print('matched ', [L, S])

                            matchedL = L
                            matchedS = S
                            break



                
                '''
                combinedTerms = [[0], [0]]
                for subShell in subShellTerms:
                    LArr = [x[0] for x in subShell]
                    SArr = [x[1] for x in subShell]

                    maxL = max([x + y for x in LArr for y in combinedTerms[0]])
                    minL = min([abs(x - y) for x in LArr for y in combinedTerms[0]])

                    maxS = max([x + y for x in SArr for y in combinedTerms[1]])
                    minS = min([abs(x - y) for x in SArr for y in combinedTerms[1]])
                    
                    def fRange(x, y):
                        while x <= y:
                            yield x
                            x += 1.

                    if minL != maxL:
                        combinedTerms[0] = list(fRange(minL, maxL))
                    else:
                        combinedTerms[0] = list(minS)

                    if minS != maxS:
                        combinedTerms[1] = list(fRange(minS, maxS))
                    else:
                        combinedTerms[0] = list(minS)

                print(combinedTerms)
                
                for L in coreL:
                    for S in coreS:
                        print([L, S])
                        if L in combinedTerms[0] and S in combinedTerms[1]:
                            print('matched', [L, S])
                            break
                    
                '''
                coreN = core.maxN
                core.term = Term(matchedL, matchedS, coreN)
            else:
                n = core.maxN
                core.term = Term(L, S, n)
            


            term = Term(L, S, n)
            if excitedState:
                print(excitedState.subShell.getSubShellString(), term.getTermString(), 'core - ', core.atomicNotation, core.term.getTermString())
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
                    excitedShell = SubShell(n, L, 1)
                    excitedState = ExcitedState(excitedShell)
                    SArr = clebschGordon(core.maxS, excitedState.S)
                    n = excitedState.subShell.n

                    for S in SArr:
                        term = Term(L, S, n)
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


if __name__ == '__main__':
    O = species.Species('O')
    NIST = readNISTSpectra(O)