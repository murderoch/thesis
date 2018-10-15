import util

import matplotlib.pyplot as plt 
from math import exp
from collections import defaultdict

class Thermo:

    def __init__(self, species, completeLevels, tempRange, cutoffEnergy):
        self.species = species
        self.completeLevels = completeLevels
        self.tempRange = tempRange
        self.cutoffEnergy = cutoffEnergy

        self.constants = util.Constants()

        self.epsilon = []
        self.J = []

        for config in self.completeLevels:
            for level in config.term.levels:
                if level.energy is not None:
                    self.epsilon.append(level.energy * self.constants.Cm_1ToJoules)
                    self.J.append(level.J)

    
    def getCpRange(self):
        CpRange = self.calcCpRange(self.tempRange)
        return CpRange

    def calcSpecificHeat(self, T):
        # function which evaluates the specific heat (J / mol K))
        # T - temperature at which to evaluate the specific heat (K)
        # cutoff - index of maximum energy level   

        self.cutoffIdx = self.calcCutoffIdx(T)
        #print(self.cutoffIdx)

        def Q():
            Q = 0.
            for m in range(self.cutoffIdx):
                Q += (2.*self.J[m]+1.) * exp(-self.epsilon[m] /(self.constants.kB * T))
            return Q
        
        def Qdot():
            Qdot = 0.
            for m in range(self.cutoffIdx):
                Qdot += (2.*self.J[m]+1.) * self.epsilon[m] \
                        * exp(-self.epsilon[m] /(self.constants.kB*T))

            return 1.0/(self.constants.kB * T**2.0) * Qdot

        def Qdot2():
            Qdot2 = 0.
            for m in range(self.cutoffIdx):
                Qdot2 += (2.*self.J[m]+1.) * self.epsilon[m] \
                        * (self.epsilon[m] - 2.*self.constants.kB*T) \
                        * exp(-self.epsilon[m] / (self.constants.kB*T))
                        
            return Qdot2 / (self.constants.kB**2.0 * T**4.0) 

        Cp = T**2*Qdot2()/Q() - (T*Qdot()/Q())**2 + 2.0*T*Qdot()/Q() + 5./2.
        return Cp * self.constants.R


    def calcCutoffIdx(self, T):
        # Capitelli style fixed cutoffEnergy method
        cutoffJoules = (self.species.Io - self.cutoffEnergy) * self.constants.Cm_1ToJoules

        # Gordon and McBride 'TEMPER' Method
        #cutoffJoules = self.species.Io*self.constants.Cm_1ToJoules - T * self.constants.kB #/ self.cm_1ToJoules
        
        cutoff = len(self.epsilon)

        #print(cutoffJoules)
        for e, energy in enumerate(self.epsilon):
            #print(energy, cutoffJoules)
            if energy > cutoffJoules:
                cutoff = e - 1
                break
        if cutoff < 1:
            cutoff = 1
        
        #print(cutoff)
        return cutoff


    def calcCpRange(self, temps):
        # calculates a range of Cp values (J / mol K)) for a range of temperatures
        # returns a 2D array of temperatures and values

        # temps - array of input temperatures (K)
        self.CpArray = []

        for T in temps:
            self.CpArray.append(self.calcSpecificHeat(T))
            
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

    #he = Species('he', 1000)

    #for level in he.states:
    #    print(level)

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
    plt.show()