
def getQuantumNumberArray(x, y=0):
    QNumber = []
    if y != 0:
        i = x + y
        iMin = abs(x - y)
    else:
        i = x
        iMin = 0

    while i >= iMin:
        QNumber.append(i)
        i = i - 1
    return QNumber


ElectronConfigMap = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'q', 'r', 't', 'u', 'v', 'w', 'x', 'y', 'z']
TermMap = [item.upper() for item in ElectronConfigMap]


'''
noElectrons = 0

for shell in shellConfigs:
    noElectrons += int(shell[-1])
'''
#nMax = int(shellConfigs[-1][0])
#SArr = getQuantumNumberArray(noElectrons * 0.5)
#LArr = getQuantumNumberArray(nMax - 1)

class SubShell:
    def __init__(self, input):
        self.n = input[0]
        self.L = input[1]
        self.NoElectrons = input[-1]

        self.calculateState
    
    def calculateState(self):
        if self.NoElectrons > 1 and self.L == 0:
            self.S = 2

    def maxElectrons(self):
        return self.L * 4. + 2.

class Shell:
    def __init__(self, input):
        self.n = input[0]
    
    def populateSubshells(self):
        pass


config = '1s1.2s1'      # get from atom species

#shellConfigs = config.split('.')
shellConfigs = []
for shellString in config.split('.'):
    shellConfigs.append(SubShell(shellString))

for shell in shellConfigs:
    print(shell.n, shell.L, shell.NoElectrons)

LMin = ElectronConfigMap.index(shellConfigs[0][1])
LMax = ElectronConfigMap.index(shellConfigs[-1][1])

n = 2
noElectrons = 1

SMax = 1
SMin = 0

LArr = getQuantumNumberArray(LMax, LMin)
SArr = getQuantumNumberArray(SMax, SMin)

print(LArr, SArr)

terms = []

for L in LArr:
    for S in SArr:
        JArr = getQuantumNumberArray(L, S)
        for J in JArr:
            term = str(2*S+1) + str(LTermMap[L]) + str(J)
            terms.append(term)
print(terms)

'''
excitedConfigs = []

for i in range(n):
    term = str(n) + LTermMap[i]
    excitedConfigs.append(term)


SCore = (float(coreconfig[0])- 1.0)/2.
LCoreElectronConfig.index(coreconfig[1])



#print(SCore)

SExcited = 0.5
LExcited = LTermMap.index(excitedConfig[1])

SArr = getQuantumNumberArray(SCore, SExcited)
LArr = getQuantumNumberArray(LCore, LExcited)

#print(SArr, LArr)   
'''