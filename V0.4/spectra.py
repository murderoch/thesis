import util
import species

#constants = util.Constants()


class Configuration:
    def __init__(self, core, excitedState, term):
        self.core = core
        self.excitedState = excitedState
        self.term = term

        #self.S = self.ClebschGordon(self.core.SMax, self.excitedState.S)
        #self.L = self.ClebschGordon(self.core.LMax, self.excitedState.L)

        '''
        self.terms = []
        for L in self.L:
            for S in self.S:
                self.terms.append(Term(L, S, n))
        '''
        if self.excitedState:
            self.ID = self.core.atomicNotation + '.' + self.excitedState.subShell.getSubShellString() + '--' + term.getTermString()
        else:
            self.ID = self.core.atomicNotation + '--' + term.getTermString()
    '''
    def ClebschGordon(self, x, y):
        i = x + y
        iMin = abs(x - y)
        QNumberArray = []
        while i >= iMin:
            QNumberArray.append(i)
            i = i - 1
        return QNumberArray
    '''

class Core:
    def __init__(self, configuration, energy):
        
        self.configuration = configuration
        self.noElectrons = 0    
        self.atomicNotation = ''
        self.maxN = 0
        self.maxL = 0
        self.energy = energy

        for subShell in configuration:
            self.noElectrons += subShell.noElectrons
            self.atomicNotation = self.atomicNotation + '.' + subShell.getSubShellString()
            if subShell.n > self.maxN:
                self.maxN = subShell.n
            if subShell.L > self.maxL:
                self.maxL = subShell.L 
        #print(self.LMax, '\n\n')
        self.atomicNotation = self.atomicNotation[1:]
        #self.SMax = self.noElectrons * 0.5
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
        '''
        self.S = self.getS()

        self.terms = []
        for S in self.S:
            self.terms.append(Term(self.L, S, self.n))
        '''
    def getSubShell(self):
        return [self.n, self.L, self.noElectrons]

    def getSubShellString(self):
        return str(self.n) + str(util.subShellMap[self.L]) + str(self.noElectrons)
    
    '''
    def getS(self):
        S = []
        i = self.noElectrons * 0.5
        while i >= 0:
            S.append(i)
            i = i - 1
        return S
    '''

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
        #self.J = self.L + self.S
    
        self.levels = []
        for J in self.J:
            self.levels.append(Level(J, None))


    def getJ(self):
        JMax = self.L + self.S
        JMin = abs(self.L - self.S)

        #print(JMax, JMin)

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
        #self.L = L
        #self.S = S
        #self.n = n
        #self.energy = self.hydrogenicApprox()
        #self.energy = self.ritzRydberg()

    def setEnergy(self, energy):
        self.energy = energy


class CalcEnergy:
    def __init__(self, species, configuration, level):
        self.species = species
        self.configuration = configuration
        self.level = level
        self.constants = util.Constants()

        self.energy = self.hydrogenicApprox()

    def getEnergy(self):
        return self.energy

    def hydrogenicApprox(self):
        IH = self.constants.IH
        Io = self.species.Io
        if self.configuration.excitedState:
            n = self.configuration.excitedState.subShell.n
        else:
            n = self.configuration.core.maxN

        #print(IH, Io)
        energy = Io - IH/n**2
        #print(self.n, energy)#/constants.Cm_1ToJoules)
        return energy#/constants.Cm_1ToJoules

    '''
    def ritzRydberg(self):
        IH = constants.IH
        Io = species.Io
        z = species.charge
        self.A = 0
        self.B = 0
        energy = Io - (IH*(z + 1.)**2.)/(self.n + self.A + self.B/(self.n**2.))**2.0
        return energy
    '''

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

        core = Core(coreShells, 0)
        #print(core.SMax)
        n, L, noElectrons = formatConfig(configs[-1])
        excitedShell = SubShell(n, L, noElectrons)
        excitedState = ExcitedState(excitedShell)
        #print(core.atomicNotation, excitedState.subShell.getSubShellString())

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
            elif row[0:5] == 'Limit':
                ionizationEnergy = float(row.rstrip('\n').split(', ')[1])
                break

            values = row.rstrip('\n').split(', ')
            #print(values)
            core, excitedState = readConfiguration(values[0])
            #SArr = clebschGordon(core.maxS, excitedState.S)
            #LArr = range(0, max(clebschGordon(core.maxL, excitedState.L))+1)
            n = excitedState.subShell.n
            '''
            for L in LArr:
                #print(L)
                for S in SArr:
                    #print(L, S)
                    termL, termS = readTerm(values[1])
                    term = Term(termL, termS, n)
                    #print(term.getTermString())
                    configuration = Configuration(core, excitedState, term)
                    
                    if configuration.ID not in IDs:
                        IDs.append(configuration.ID)
                        NIST.append(configuration)
                    else:
                        configuration = NIST[IDs.index(configuration.ID)]

                    L, S = readTerm(values[1])
                    for level in configuration.term.levels:
                        if configuration.term.L == L and configuration.term.S == S and level.J == float(values[2]):
                            level.setEnergy(float(values[3]))
            '''

            L, S = readTerm(values[1])
            term = Term(L, S, n)
            configuration = Configuration(core, excitedState, term)
            if configuration.ID not in IDs:
                IDs.append(configuration.ID)
                NIST.append(configuration)
            else:
                configuration = NIST[IDs.index(configuration.ID)]
            #print(configuration.term.levels)
            for level in configuration.term.levels:
                #print(level)
                if configuration.term.L == L and configuration.term.S == S and level.J == float(values[2]):
                    #print(configuration.ID, level.J)
                    level.setEnergy(float(values[3]))
                    #print(level.energy)
    return NIST

def calculateExpectedStates(species, nMax):
    states = []

    cores = []
    for core in species.cores:
        coreShells = []
        for subShell in core[0]:
            coreShells.append(SubShell(subShell[0], subShell[1], subShell[2]))
        cores.append(Core(coreShells, core[1]))

    for core in cores:

        if core.noElectrons == species.noElectrons:
            n = core.maxN
            LArr = range(0, core.maxL+1, 1)
            SArr = []

            i = core.maxS
            #print(core.atomicNotation, core.maxS)
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
                    #LArr = clebschGordon(core.LMax, excitedState.L)
                    n = excitedState.subShell.n

                    #print(core.SMax, excitedState.S)
                    #for L in LArr:
                    for S in SArr:
                        term = Term(L, S, n)
                        #print(S)
                        #rint(term.getTermString())
                        configuration = Configuration(core, excitedState, term)
                        #print(excitedState.L)
                        #print(excitedState.subShell.getSubShellString())
                        states.append(configuration)
        
    return states

def calculateEnergies(nMax, species, NIST, theory):

    for config in theory:
        #print(config.ID)
        NISTConfig = next((x for x in NIST if x.ID == config.ID), None)
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
                    energy = CalcEnergy(species, config, level)
                    level.setEnergy(energy.getEnergy())
        else:
            for level in config.term.levels:
                #print(config.ID)
                energy = CalcEnergy(species, config, level)
                level.setEnergy(energy.getEnergy())
    
    return theory


nMax = 10
O = species.Species('O')
#He = species.Species('He')

use = O

NIST = readNISTSpectra(use)
theory = calculateExpectedStates(use, nMax)


'''
print('')
for config in NIST:
    for level in config.term.levels:
        print('NIST', config.ID, level.J, level.energy)

print('')
for config in theory:
    print('Thry', config.ID)
'''


Olevels = calculateEnergies(nMax, use, NIST, theory)
#energies = [[[x for x in level.energy] for level in config.term.levels] for config in Olevels]

#for config in Olevels:
    #for level in config.term.levels:
        #print(level.energy)
        #print(type(level.energy))

energies = [level.energy for config in Olevels for level in config.term.levels]

configs = [config for config in Olevels for config in config.term.levels]

#asdf = zip(configs, energies)
#asdf = list(zip())
#asdf = [x for _,x in *sorted(zip(configs, energies))]


for x in energies:
    if type(x) != float:
        print(type(x))

'''
values = [[level for level in config.term.levels] for config in Olevels]
energies = [[level.energy for level in config.term.levels] for config in Olevels]

energy = [energy for energy in energies]
print(energy)

zipped = zip(values, energies)
for level in zipped:
    print(level[1])


#print(sorting)
#Olevels = sorting
#sortedLevels = sorted(zip(Olevels, [x for x in Olevels.configuration.terms.levels.energy]))
#Olevels = NIST
#Olevels = theory
'''



with open('output.dat', 'w') as output:
    for config in Olevels:
        for level in config.term.levels:
            string = config.ID + ', ' + str(level.J) + ', ' + str(level.energy) + '\n'
            output.write(string)


'''
for config in NIST:
    for level in config.term.levels:
        #print(config.ID, level.energy)
        continue

for config in theory:
    NISTConfig = next((x for x in NIST if x.ID == config.ID), None)
    if NISTConfig:
        for level in config.term.levels:
            NISTlevel = next((x for x in NISTConfig.term.levels if x.J == level.J), None)
            if NISTlevel:
                level.setEnergy(NISTlevel.energy)
            else:
                print(config.ID)
                energy = CalcEnergy(He, config, level)
                level.setEnergy(energy.getEnergy())
    else:
        for level in config.term.levels:
            print(config.ID)
            energy = CalcEnergy(He, config, level)
            level.setEnergy(energy.getEnergy())
'''


