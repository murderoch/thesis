import numpy as np 
from matplotlib import pyplot as plt 
from math import *

kB = 1.38064852E-23
R = 8.3144598
cm_1ToJoules = 1.23986E-4 * 1.60217E-19

J = [2, 1, 0, 2, 0, 2, 1]
g = [5, 3, 1, 5, 1, 5, 3]
epsilon = [0.0, 158.5, 226.5, 15867.7, 33792.4, 73767.81, 76794.69] #cm-1

                #g  J  E(cm-1)
spectralData = np.array([
                [2, 0.0],
                [1, 158.5],
                [0, 226.5],
                
                [2, 15867.7],
                
                [0, 33792.4],
                
                [2, 73767.81],
                
                [1, 76794.69],
                
                [1, 86625.35],
                [2, 86627.37],
                [3, 86631.04],
                
                [2, 88630.84],
                [1, 88630.30],
                [0, 88631.00],
             
                [2, 95476.43],
                [1, 96225.50],
                
                [4, 97420.24],
                [3, 97420.37],
                [2, 97420.37],
                [2, 97420.50],
                [1, 97420.50],
                [0, 97420.50],

                [3, 97488.14],
                [2, 97488.14],
                [1, 97488.14],

                [2, 99690.40],
                [1, 99690.40],
                [0, 99690.40],

                [3, 101135.04],
                [2, 101147.21],
                [1, 101155.10],

                [2, 102411.21],

                [1, 102411.65],

                [2, 102661.64],

                [4, 102865.09],

                [3, 102908.14],
                [2, 102908.14],
                [1, 102908.14],

                [2, 103869.4],
                [1, 103869.4],
                [0, 103869.4],

                [2, 105019.0],

                [1, 105164.9],

                [4, 105385.3],
                [3, 105385.3],
                [2, 105385.3],
                [1, 105385.3],
                [0, 105385.3],

                [3, 105408.58],
                [2, 105408.58],
                [1, 105408.58],

                [2, 105911.3],
                [1, 105911.3],
                [0, 105911.3],

                [2, 106545.1],

                [1, 106627.9],

                [4, 106751.2],
                [3, 106751.2],
                [2, 106751.2],
                [1, 106751.2],
                [0, 106751.2],
                ])
                
def specificHeat(T, L):
   
    def Q():
        Q = 0.
        for m in range(L):
            Q += (2.*spectralData[m,0]+1.) * exp(-spectralData[m,1]*cm_1ToJoules /(kB * T))
        return Q
    
    def Qdot():
        Qdot = 0.
        for m in range(L):
            Qdot += (2.*spectralData[m,0]+1.) * spectralData[m,1]*cm_1ToJoules \
                    * exp(-spectralData[m,1]*cm_1ToJoules /(kB*T))
        return 1.0/(kB * T**2.0) * Qdot

    def Qdot2():
        Qdot2 = 0.
        for m in range(L):
            Qdot2 += (2.*spectralData[m,0]+1.) * spectralData[m,1]*cm_1ToJoules \
                     * (spectralData[m,1]*cm_1ToJoules - 2.*kB*T) \
                     * exp(-spectralData[m,1]*cm_1ToJoules / (kB*T))
        return Qdot2 / (kB**2.0 * T**4.0)            
    
    Q = Q()
    Qdot = Qdot()
    Qdot2 = Qdot2()

    Cp = T**2*Qdot2/Q - (T*Qdot/Q)**2 + 2.0*T*Qdot/Q + 5./2.
    
    return Cp * R #/ 4.183

temps = np.linspace(500, 50000, 50)
bethe = []
russian = []
temper = []
capitelli = []

for temp in temps:
    nMax = int(round(2.461 * temp ** (1./6.),0))
    
    ionizationPotential = 13.61806 * 1.60218e-19 - kB*temp
    
    if ionizationPotential > max(spectralData[:,1]*cm_1ToJoules):
        loweredEnergyLevel = len(spectralData[:,1])
    else:
        loweredEnergyLevel = next(x[0] for x in enumerate(spectralData[:,1]*cm_1ToJoules) if x[1] > ionizationPotential)
    temper.append(specificHeat(temp, loweredEnergyLevel))
    
    ionizationPotential = 13.61806 * 1.60218e-19 - 1000*cm_1ToJoules
    #print(max(spectralData[:,1]*cm_1ToJoules))
    if ionizationPotential > max(spectralData[:,1]*cm_1ToJoules):
        loweredEnergyLevel = len(spectralData[:,1])
    else:
        loweredEnergyLevel = next(x[0] for x in enumerate(spectralData[:,1]*cm_1ToJoules) if x[1] > ionizationPotential)
    capitelli.append(specificHeat(temp, loweredEnergyLevel))


    russian.append(specificHeat(temp, 5))
    bethe.append(specificHeat(temp, nMax))


print(specificHeat(300, 5))

Gordon1999 = [[100, 200, 500, 700, 1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000],
               [23.703, 22.734, 21.257, 21.040, 20.915, 20.827, 20.937, 21.799, 22.708, 23.150, 23.868, 24.722]]

Capitelli2005 = [[100, 200, 500, 700, 1000, 2000, 3000, 5000, 10000, 12000, 13000, 14000,
                  15000, 16000, 17000, 18000, 19000, 20000, 22000, 23000, 24000, 25000, 26000, 27000, 28000,
                  30000, 34000, 40000, 44000, 50000],
               [23.70, 22.74, 21.26, 21.04, 20.91, 20.83, 20.94, 21.80, 23.29, 24.35, 25.74, 28.15,
                31.89, 37.19, 44.04, 52.19, 61.11, 70.00, 84.37, 88.54, 90.36, 89.94, 87.66, 84.00, 79.46,
                69.37, 51.52, 35.92, 30.6, 26.30]]

plt.plot(temps, russian, 'r')
plt.plot(temps, temper, 'g')
plt.plot(temps, bethe, 'b')
plt.plot(Gordon1999[0], Gordon1999[1], 'k')
plt.plot(temps, capitelli, 'violet')
plt.plot(Capitelli2005[0], Capitelli2005[1], 'k')

plt.legend(['Russian Cuttoff', 'Temper Cuttoff', 'Bethe Cuttoff', 'Gordon 1999', 'Capitelli 2005'], loc=4)
plt.xlabel('Temperature (K)')
plt.ylabel('Cp (cal / g mol K)')
plt.ylim([0,30])
plt.grid()
plt.title('Comparison of Monatomic Cuttoff Functions for Atomic Oxygen')
plt.show(block=True)

'''
eilmer4 = np.genfromtxt('O-thermo.dat', skip_header=1, names=True, dtype=float, delimiter = ' ')
eilmer4ArrT = []
eilmer4ArrCp = []

CpConversion = 4.184/0.032

for temp in eilmer4[:]:
    eilmer4ArrT.append(temp[0])
    eilmer4ArrCp.append(temp[1]/1.987/CpConversion)
''' 
