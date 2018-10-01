import itertools

class Configuration:
    def __init__(self, subShells):
        
        self.noElectrons = 0
        self.LMax = 0
        self.nMax = 0
        self.subShells = subShells
        self.atomicNotation = ''

        for subShell in subShells:
            self.noElectrons += subShell.noElectrons
            self.LMax += subShell.noElectrons * subShell.L
            if subShell.n > self.nMax:
                self.nMax = subShell.n

            self.atomicNotation = self.atomicNotation + '.' + subShell.getSubShellString()

        self.L = self.getL()
        self.S = self.getS()

        self.terms = []

        for L in self.L:
            for S in self.S:
                self.terms.append(Term(L, S))

    def getS(self):
        S = []
        i = self.noElectrons * 0.5
        while i >= 0:
            S.append(i)
            i = i - 1
        return S

    def getL(self):
        L = []
        i = self.LMax
        while i >= 0:
            L.append(i)
            i = i - 1
        return L

    def ClebschGordon(self, x, y=None):
        i = sum(x)
        if y == None:
            iMin = abs(x)
        QNumberArray = []
        while i >= iMin:
            QNumberArray.append(i)
            i = i - 1
        return QNumberArray

class SubShell():
    def __init__(self, n, L, noElectrons):
        self.n = n
        self.L = L
        self.noElectrons = noElectrons
        self.subShellMap = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', \
                            'l', 'm', 'n', 'o', 'q', 'r', 't', 'u', 'v', \
                            'w', 'x', 'y', 'z']

    def getSubShell(self):
        return [self.n, self.L, self.noElectrons]

    def getSubShellString(self):
        return str(self.n) + str(self.subShellMap[self.L]) + str(self.noElectrons)

class Term:
    def __init__(self, L, S):
        self.L = L
        self.S = S

        self.J = self.getJ()

        self.termMap = ['S', 'P', 'D', 'F', 'G', 'H', 'I', 'J', 'K', \
                        'L', 'M', 'N', 'O', 'Q', 'R', 'T', 'U', 'V', \
                        'W', 'X', 'Y', 'Z']
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

nMax = 2
noElectrons = 1

subShells = []
for E in range(1, noElectrons+1, 1):
    for L in range(0, nMax, 1):
        for n in range(1, nMax+1, 1):
            if E > (4. * L + 2):
                break
            print([n, L, E])
            subShells.append([n, L, E])       

def addShell(oldShell, E, subShells, noElectrons):
    if E >= noElectrons:
        return oldShell
    else:
        for subShell in subShells:
            if not subShell in oldShell:
                if E + subShell[2] <= noElectrons:
                    oldShell.append(subShell)
                    E = E + subShell[2]
                    oldShell = addShell(oldShell, E, subShells, noElectrons)
                    return oldShell

'''
newShells = []
for i in range(len(subShells)):
    for subShell in modSubShells:
        newShells.append(addShell([subShell], subShell[2], subShells, noElectrons))

for shell in newShells:
    print(shell)


for subShell in subShells:
    print(subShell)
'''

'''
def getLArray(n):
    LMax = n-1
    L = []
    i = LMax
    while i >= 0:
        L.append(i)
        i = i - 1
    return L

def getElectronArray(L):
    EMax = 4.*L + 2.
    E = []
    i = EMax
    while i >= 0:
        E.append(i)
        i = i - 1
    return E

configurations = []

for n in range(1, nMax+1, 1):
    totalElectrons = 0
    LArr = getLArray(n)
    for L in LArr:
        EArr = getElectronArray(L)
        subShell = []
        for E in EArr:
            #print(n, L, E)
            subShell.append(SubShell(n, L, noElectrons))
        configuration = Configuration(subShell)
        configurations.append(configuration)

for configuration in configurations:
    for term in configuration.terms:
        for subShell in configuration.subShells:
            #print(subShell.getSubShellString(), term.getTermString())
            #print(subShell.getSubShellString())
            pass
    #print(configuration.atomicNotation)   
'''