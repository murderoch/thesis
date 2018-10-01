import matplotlib.pyplot as plt 
from math import exp

from collections import defaultdict

import util

class EnergyState:
    
    def __init__(self, config, term, j, energy):
        self.j = float(j)
        self.energy = float(energy)*constants.Cm_1ToJoules
        self.config = config
        self.term = term

        if len(config.split(')')) == 2:
            self.coreConfig = config.split('(')[0][0:-1]
            self.coreTerm = '(' + config.split('(')[1].split(')')[0] + ')'
            self.excitedConfig = config.split(')')[1].replace('.', '')
        else:
            self.coreConfig = None
            self.coreTerm = None
            self.excitedConfig = None
    
        self.__name__ = self.config + self.term
        
class Species:

    def __init__(self, species, cutoffCm):
        self.charge = species.count('i') - 1
        self.species = species
        self.readSpectra()
        #self.ionizationEnergy = self.epsilon[-1]
        self.cutoffEnergy = cutoffCm * constants.Cm_1ToJoules
        self.cutoffIdx = self.calcCutoffIdx(-1)
        
        #del self.epsilon[-1]
        #del self.j[-1]

        self.calcSpectra()

    def readSpectra(self):
        self.j = []
        self.epsilon = []
        self.states = defaultdict(list)

        filename = util.getSpectralDataDir() + self.species + '.dat'

        with open(filename) as dataFile:
            next(dataFile)
            for row in dataFile:
                if not row.strip():
                    continue
                elif row[0:5] == 'Limit':
                    self.ionizationEnergy = float(row.rstrip('\n').split(', ')[1])
                    break

                values = row.rstrip('\n').split(', ')
                
                self.j.append(float(values[2]))
                self.epsilon.append(float(values[3])*constants.Cm_1ToJoules)

                # create new energy state for imported values
                state = EnergyState(values[0], values[1], values[2], values[3])
                # append to dict organised by core configuration
                self.states[state.config].append(state)

    def calcSpectra(self):
        #for state in coreState for coreState in list(self.states.values()):
        if 0:
            for coreState in self.states.values(): 
                for state in coreState:
                    print(state.config)
                    print(state.coreConfig, state.coreTerm)
        
        test = self.states['2s2.2p2']

        def calcRitzRydberg(self, state):

            A = 2
            B = 3

            en = self.ionizationEnergy - constants.IH * (self.charge + 1.)**2 \
                / (state.n + A + B/(state.n**2))**2

            return en

    def calcSpecificHeat(self, T):
        # function which evaluates the specific heat (J / mol K))
        # T - temperature at which to evaluate the specific heat (K)
        # cutoff - index of maximum energy level   

        #print(constants.kB*T/self.cm_1ToJoules, 
        #self.epsilon[self.cutoffIdx]/self.cm_1ToJoules)

        #cutoffIdx = self.calcCutoffIdx(T)
        self.cutoffIdx = self.calcCutoffIdx(T)
        #print(self.cutoffIdx)

        def Q():
            Q = 0.
            for m in range(self.cutoffIdx):
                Q += (2.*self.j[m]+1.) * exp(-self.epsilon[m] /(constants.kB * T))
            return Q
        
        def Qdot():
            Qdot = 0.
            for m in range(self.cutoffIdx):
                Qdot += (2.*self.j[m]+1.) * self.epsilon[m] \
                        * exp(-self.epsilon[m] /(constants.kB*T))

            return 1.0/(constants.kB * T**2.0) * Qdot

        def Qdot2():
            Qdot2 = 0.
            for m in range(self.cutoffIdx):
                Qdot2 += (2.*self.j[m]+1.) * self.epsilon[m] \
                        * (self.epsilon[m] - 2.*constants.kB*T) \
                        * exp(-self.epsilon[m] / (constants.kB*T))
                        
            return Qdot2 / (constants.kB**2.0 * T**4.0) 

        Cp = T**2*Qdot2()/Q() - (T*Qdot()/Q())**2 + 2.0*T*Qdot()/Q() + 5./2.

        return Cp * constants.R


    def calcCutoffIdx(self, T):

        # Capitelli style fixed cutoffEnergy method
        cutoffJoules = self.ionizationEnergy - self.cutoffEnergy

        # Gordon and McBride 'TEMPER' Method
        #cutoffJoules = self.ionizationEnergy - T * constants.kB #/ self.cm_1ToJoules

        cutoff = 0 
        for e, energy in enumerate(self.epsilon):
            #print(energy)
            if energy > cutoffJoules:
                cutoff = e - 1
                break
        #print(cutoff)
        return cutoff


    def calcCpRange(self, temps):
        # calculates a range of Cp values (J / mol K)) for a range of temperatures
        # returns a 2D array of temperatures and values

        # temps - array of input temperatures (K)
        self.CpArray = [[],[]]

        for T in temps:
            self.CpArray[0].append(T)
            self.CpArray[1].append(self.calcSpecificHeat(T))
            
        return self.CpArray

if __name__ == '__main__':
    Capitelli2005oi = [[100, 200, 500, 700, 1000, 2000, 3000, 5000, 10000, 12000, 13000, 14000,
                    15000, 16000, 17000, 18000, 19000, 20000, 22000, 23000, 24000, 25000, 26000, 27000, 28000,
                    30000, 34000, 40000, 44000, 50000],
                [23.70, 22.74, 21.26, 21.04, 20.91, 20.83, 20.94, 21.80, 23.29, 24.35, 25.74, 28.15,
                    31.89, 37.19, 44.04, 52.19, 61.11, 70.00, 84.37, 88.54, 90.36, 89.94, 87.66, 84.00, 79.46,
                    69.37, 51.52, 35.92, 30.6, 26.30]]

    Gordon1999 = [[100, 200, 500, 700, 1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000],
                [23.703, 22.734, 21.257, 21.040, 20.915, 20.827, 20.937, 21.799, 22.708, 23.150, 23.868, 24.722]]

    constants = util.Constants()

    #oi = Species('oi', 1000)
    #oii = Species('oii', 1000)
    #ni = Species('ni', 1000)
    #nii = Species('nii', 1000)
    #ci = Species('ci', 1000)

    #oi.calcCpRange(range(1000, 51000, 1000))
    #oii.calcCpRange(range(1000, 51000, 1000))
    #ni.calcCpRange(range(1000, 51000, 1000))
    #nii.calcCpRange(range(1000, 51000, 1000))

    he = Species('he', 1000)

    for level in he.states:
        print(level)

    #plt.figure()

    #plt.plot(oi.CpArray[0], oi.CpArray[1], label='o-i')
    #plt.plot(oii.CpArray[0], oii.CpArray[1], label='o-ii')
    #plt.plot(ni.CpArray[0], ni.CpArray[1], label='n-i')
    #plt.plot(nii.CpArray[0], nii.CpArray[1], label='n-ii')

    #plt.plot(Gordon1999[0], Gordon1999[1], label='Gordon Coeffs')
    plt.plot(Capitelli2005oi[0], Capitelli2005oi[1], label='Capitelli o-i')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Cp (J/mol/K)')

    plt.legend()
    #plt.show()