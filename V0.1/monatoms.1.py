import matplotlib.pyplot as plt 
from math import *

class species:

    def __init__(self, species, cutoffCm):
        self.kB = 1.38064852E-23
        self.R = 8.3144598
        self.cm_1ToJoules = 1.23986E-4 * 1.60217E-19

        self.g = []
        self.j = []
        self.epsilon = []

        filename = species + '.dat'

        with open(filename) as dataFile:
            next(dataFile)
            for row in dataFile:
                if not row.strip():
                    continue
 
                values = row.rstrip('\n').split(', ')
                self.j.append(float(values[0]))
                #self.g.append(float(values[0])*2+1)
                self.epsilon.append(float(values[1])*self.cm_1ToJoules)

        self.ionizationEnergy = self.epsilon[-1]
        self.cutoffEnergy = cutoffCm * self.cm_1ToJoules
        self.cutoffIdx = self.calcCutoffIdx()
        
        del self.epsilon[-1]
        del self.j[-1]

        

    def calcSpecificHeat(self, T):
        # function which evaluates the specific heat (cal/(g mol K))

        # T - temperature at which to evaluate the specific heat (K)
        # cutoff - index of maximum energy level

        def Q():
            Q = 0.
            for m in range(self.cutoffIdx):
                Q += (2.*self.j[m]+1.) * exp(-self.epsilon[m] /(self.kB * T))
            return Q
        
        def Qdot():
            Qdot = 0.
            for m in range(self.cutoffIdx):
                Qdot += (2.*self.j[m]+1.) * self.epsilon[m] \
                        * exp(-self.epsilon[m] /(self.kB*T))

            return 1.0/(self.kB * T**2.0) * Qdot

        def Qdot2():
            Qdot2 = 0.
            for m in range(self.cutoffIdx):
                Qdot2 += (2.*self.j[m]+1.) * self.epsilon[m] \
                        * (self.epsilon[m] - 2.*self.kB*T) \
                        * exp(-self.epsilon[m] / (self.kB*T))
                        
            return Qdot2 / (self.kB**2.0 * T**4.0)     

        Q = Q()
        Qdot = Qdot()
        Qdot2 = Qdot2()
        
        Cp = T**2*Qdot2/Q - (T*Qdot/Q)**2 + 2.0*T*Qdot/Q + 5./2.

        return Cp * self.R

    def calcCutoffIdx(self):
        cutoffJoules = self.ionizationEnergy - self.cutoffEnergy
        cutoff = 0

        for e, energy in enumerate(self.epsilon):
            #print(energy, cutoffJoules)
            if energy > cutoffJoules:
                cutoff = e - 1
                break

        return cutoff

    def calcCpRange(self, temps):
        # calculates a range of Cp values (cal/(g mol K)) for a range of temperatures
        # returns a 2D array of temperatures and values

        # temps - array of input temperatures (K)
        self.CpArray = [[],[]]

        for T in temps:
            self.CpArray[0].append(T)
            self.CpArray[1].append(self.calcSpecificHeat(T))
            
        return self.CpArray

Capitelli2005 = [[100, 200, 500, 700, 1000, 2000, 3000, 5000, 10000, 12000, 13000, 14000,
                  15000, 16000, 17000, 18000, 19000, 20000, 22000, 23000, 24000, 25000, 26000, 27000, 28000,
                  30000, 34000, 40000, 44000, 50000],
               [23.70, 22.74, 21.26, 21.04, 20.91, 20.83, 20.94, 21.80, 23.29, 24.35, 25.74, 28.15,
                31.89, 37.19, 44.04, 52.19, 61.11, 70.00, 84.37, 88.54, 90.36, 89.94, 87.66, 84.00, 79.46,
                69.37, 51.52, 35.92, 30.6, 26.30]]


#oxygen = species('o', 1000)
nitrogen = species('n', 1000)

#oxygen.calcCpRange(range(500, 50000))
nitrogen.calcCpRange(range(500, 40000, 10))

print(nitrogen.calcSpecificHeat(21000))

plt.figure()
#plt.plot(Capitelli2005[0], Capitelli2005[1])
#plt.plot(oxygen.CpArray[0], oxygen.CpArray[1])
plt.plot(nitrogen.CpArray[0], nitrogen.CpArray[1])
plt.show()