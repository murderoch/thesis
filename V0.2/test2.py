class Constants:
    def __init__(self):
        self.kB = 1.38064852E-23                        #
        self.R = 8.3144598                              #
        self.Cm_1ToJoules = 1.60217E-19 * 1.23986E-4    #
        self.eVToCm_1 = 8065.544005                     #
        self.IH = 13.59844 * self.eVToCm_1              #CRC Handbook of Chemistry and Physics (2003)

class Species:
    def __init__(self, Io):
        self.constants = Constants()
        self.Io = Io
        self.charge = 0

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

        self.subShellMap = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', \
                            'l', 'm', 'n', 'o', 'q', 'r', 't', 'u', 'v', \
                            'w', 'x', 'y', 'z']
        
        self.terms = []
        for S in self.S:
            self.terms.append(Term(self.L, S, self.n))

    def getSubShell(self):
        return [self.n, self.L, self.noElectrons]

    def getSubShellString(self):
        return str(self.n) + str(self.subShellMap[self.L]) + str(self.noElectrons)
    
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
        #self.term = Term(self.subShell.L, 0.5, subShell.n)

class Configuration:
    def __init__(self, core, excitedState, n, species):
        self.core = core
        self.excitedState = excitedState

        self.S = self.ClebschGordon(self.core.SMax, self.excitedState.S)
        self.L = self.ClebschGordon(self.core.LMax, self.excitedState.L)

        self.terms = []
        for L in self.L:
            for S in self.S:
                self.terms.append(Term(L, S, n))

    def ClebschGordon(self, x, y):
        i = x + y
        iMin = abs(x - y)
        QNumberArray = []
        while i >= iMin:
            QNumberArray.append(i)
            i = i - 1
        return QNumberArray


class Term:
    def __init__(self, L, S, n):
        self.n = n
        self.L = L
        self.S = S
        self.J = self.getJ()

        self.termMap = ['S', 'P', 'D', 'F', 'G', 'H', 'I', 'J', 'K', \
                        'L', 'M', 'N', 'O', 'Q', 'R', 'T', 'U', 'V', \
                        'W', 'X', 'Y', 'Z']

        self.levels = []
        for J in self.J:
            self.levels.append(Level(J, self.L, self.S, self.n))

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
        return str(2.*self.S + 1.) + self.termMap[self.L] + str(self.J)

class Level:
    def __init__(self, J, L, S, n):
        self.J = J
        self.L = L
        self.S = S
        self.n = n
        #self.energy = self.hydrogenicApprox()
        self.energy = self.ritzRydberg()

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


constants = Constants()
species = Species(24.58741 * constants.eVToCm_1)

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



