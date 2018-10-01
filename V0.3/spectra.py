import util
import species

#constants = util.Constants()


class Configuration:
    def __init__(self, core, excitedState):
        self.core = core
        self.excitedState = excitedState
        
        n = self.excitedState.subShell.n

        self.S = self.ClebschGordon(self.core.SMax, self.excitedState.S)
        self.L = self.ClebschGordon(self.core.LMax, self.excitedState.L)

        self.terms = []
        for L in self.L:
            for S in self.S:
                self.terms.append(Term(L, S, n))

        self.ID = []
        for term in self.terms:
            self.ID.append(self.core.atomicNotation + '-' + self.excitedState.subShell.getSubShellString() + '--' + term.getTermString())

    def ClebschGordon(self, x, y):
        i = x + y
        iMin = abs(x - y)
        QNumberArray = []
        while i >= iMin:
            QNumberArray.append(i)
            i = i - 1
        return QNumberArray


class Core:
    def __init__(self, configuration, energy):

        self.noElectrons = 0    
        self.atomicNotation = ''
        self.maxN = 0
        self.LMax = 0
        self.energy = energy

        for subShell in configuration:
            self.noElectrons += subShell.noElectrons
            self.atomicNotation = self.atomicNotation + '.' + subShell.getSubShellString()
            if subShell.n > self.maxN:
                self.maxN = subShell.n
            if subShell.L > self.LMax:
                self.LMax = subShell.L 
        
        self.atomicNotation = self.atomicNotation[1:]
        self.SMax = self.noElectrons * 0.5
        

class SubShell:
    def __init__(self, n, L, noElectrons):
        self.n = n
        self.L = L
        self.noElectrons = noElectrons
        self.S = self.getS()

        self.terms = []
        for S in self.S:
            self.terms.append(Term(self.L, S, self.n))

    def getSubShell(self):
        return [self.n, self.L, self.noElectrons]

    def getSubShellString(self):
        return str(self.n) + str(util.subShellMap[self.L]) + str(self.noElectrons)
    
    def getS(self):
        S = []
        i = self.noElectrons * 0.5
        while i >= 0:
            S.append(i)
            i = i - 1
        return S


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
        #self.L = L
        #self.S = S
        #self.n = n
        #self.energy = self.hydrogenicApprox()
        #self.energy = self.ritzRydberg()

    def setEnergy(self, energy):
        self.energy = energy

'''
    def hydrogenicApprox(self):
        IH = constants.IH
        Io = species.Io
        energy = Io - IH/self.n**2
        #print(self.n, energy)#/constants.Cm_1ToJoules)
        return energy#/constants.Cm_1ToJoules

    def ritzRydberg(self):
        IH = constants.IH
        Io = species.Io
        z = species.charge
        self.A = 0
        self.B = 0
        energy = Io - (IH*(z + 1.)**2.)/(self.n + self.A + self.B/(self.n**2.))**2.0
        return energy
'''


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
        coreShells = []
        for config in configs[:-1]: 
            n, L, noElectrons = formatConfig(config)
            coreShells.append(SubShell(n, L, noElectrons))
        
        core = Core(coreShells, 0)
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
            core, excitedState = readConfiguration(values[0])
            configuration = Configuration(core, excitedState)

            if configuration.ID not in IDs:
                IDs.append(configuration.ID)
                NIST.append(configuration)
            else:
                configuration = NIST[IDs.index(configuration.ID)]

            for term in configuration.terms:
                for level in term.levels:
                    #print(term.getTermString())
                    #print(values[2])#, values[2], values[3])
                    L, S = readTerm(values[1])
                    #print(L, S, '\n')

                    #term.J = values[2]
                    '''
                    if term.n < 3:
                        print('')
                        print(configuration.ID)
                        print(term.getTermString(), values[1])
                        print(term.L, L)
                        print(term.S, S)
                        print(level.J, values[2])
                    '''
                    if term.L == L and term.S == S and level.J == float(values[2]):
                        #print('matched')
                        #print(term.getTermString(), values)
                        level.setEnergy(float(values[3]))
                        #if term.n < 3:
                            #print('match')
                            #print(level.energy)
                        #print(term.energy)
                    #print('')
                #level = Level(int(values[2]), float(values[3]))
                #term = Term()



            #self.j.append(float(values[2]))
            #self.epsilon.append(float(values[3])*constants.Cm_1ToJoules)

            # create new energy state for imported values
            #state = EnergyState(values[0], values[1], values[2], values[3])
            # append to dict organised by core configuration
            #self.states[state.config].append(state)

    return NIST

def calculateExpectedStates(species, nMax):
    states = []
    #core = Core([SubShell(1, 0, 1)], 0)
    #print(species.cores[0])
    #core = species.cores[0]
    #print(species.cores)
    #print('asdf')
    for core in species.cores:
        #print(core)
        coreShells = []
        for subShell in core[0]:
            #print(subShell)
            coreShells.append(SubShell(subShell[0], subShell[1], subShell[2]))
        core = Core(coreShells, core[1])

    for n in range(core.maxN+1, nMax+1, 1):
        for L in range(0, nMax, 1):
            excitedShell = SubShell(n, L, 1)
            excitedState = ExcitedState(excitedShell)
            configuration = Configuration(core, excitedState)
            states.append(configuration)

    return states

nMax = 10
He = species.Species('He')
O = species.Species('O')
NIST = readNISTSpectra(He)
theory = calculateExpectedStates(He, 10)

NISTStateIDs = []
for state in NIST:
    for IDterm in state.ID:
        NISTStateIDs.append(IDterm)

for state in theory:
    for termID in state.ID:
        if termID in NISTStateIDs:


'''
if __name__ == '__main__':
    constants = util.Constants()
    species = species.Species('He')

    #N Core
    #core = Core([SubShell(1, 0, 2), SubShell(2, 0, 2), SubShell(2, 1, 2)], 0)
    #core1 = Core([SubShell(1, 0, 1), SubShell(2, 0, 1), SubShell(2, 1, 4)])
    #core1 = Core([SubShell(1, 0, 1), SubShell(2, 1, 5)])

    #He Core
    core = Core([SubShell(1, 0, 1)], 0)

    nMax = 20

    configurations = []

    for n in range(core.maxN+1, nMax+1, 1):
        for L in range(0, nMax, 1):
            excitedShell = SubShell(n, L, 1)
            excitedState = ExcitedState(excitedShell)
            configuration = Configuration(core, excitedState, n, species)
            configurations.append(configuration)

    #print('')

    with open('output.dat', 'w') as output:
        output.write('Core Config, Excited Config, Term, J, Energy\n')
        for configuration in configurations:
            terms = []
            for i, term in enumerate(configuration.terms):
                #terms.append([term.getTermString()])
                for level in term.levels:
                    #terms[i].append(level.energy)
            
                    outputString = configuration.core.atomicNotation + ', ' + configuration.excitedState.subShell.getSubShellString() + ', '\
                                + term.getTermString() + ', ' + str(level.J) + ', ' + str(level.energy) + '\n'
                    output.write(outputString)
            #print(configuration.core.atomicNotation, '+', configuration.excitedState.subShell.getSubShellString(), terms)
                #print(configuration.excitedState.subShell.getSubShellString(),configuration.excitedState.term.getTermString())



'''