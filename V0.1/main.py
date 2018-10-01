import numpy as np 
from matplotlib import pyplot as plt 
from math import *

kB = 1.38E-23
R = 8.3144598
cm_1ToJoules = 1.23986E-4 * 1.60217E-19

J = [2, 1, 0, 2, 0, 2, 1]
g = [5, 3, 1, 5, 1, 5, 3]
epsilon = [0.0, 158.5, 226.5, 15867.7, 33792.4, 73767.81, 76794.69] #cm-1

def specificHeat(T):

    L = 7 - 1 # Remember arrays start at zero!

    def Q():
        Q = 0.
        for m in range(L):
            Q += g[m] * exp(-epsilon[m]*cm_1ToJoules /(kB * T))
        return Q
    
    def Qdot():
        Qdot = 0.
        for m in range(L):
            Qdot += g[m] * epsilon[m]*cm_1ToJoules * exp(-epsilon[m]*cm_1ToJoules /(kB*T))
        return 1.0/(kB * T**2.0) * Qdot

    def Qdot2():
        Qdot2 = 0.
        for m in range(L):
            Qdot2 += g[m] * epsilon[m]*cm_1ToJoules * (epsilon[m]*cm_1ToJoules - 2.*kB*T) \
                    * exp(-epsilon[m]*cm_1ToJoules / (kB*T))
        return Qdot2 / (kB**2.0 * T**4.0)
    
    Q = Q()
    Qdot = Qdot()
    Qdot2 = Qdot2()

    #print("Q = ", Q)
    #print("Qdot = ", Qdot)
    #print("Qdot2 = ", Qdot2)

    #print("term 1 = ", T**2*Qdot2/Q)
    #print("term 2 = ", (T*Qdot/Q)**2)
    #print("term 3 = ", 2.0*T*Qdot/Q)

    Cp = T**2*Qdot2/Q - (T*Qdot/Q)**2 + 2.0*T*Qdot/Q + 5./2.
    return Cp * R / 4.183

print("Cp(500) = ", specificHeat(500))

xplot = np.linspace(500, 6000, 20)
yplot = []

for i in range(len(xplot)):
    yplot.append(specificHeat(xplot[i]))

eilmer4 = np.genfromtxt('O-thermo.dat', skip_header=1, names=True, dtype=float, delimiter = ' ')
eilmer4ArrT = []
eilmer4ArrCp = []

CpConversion = 4.184/0.032

for temp in eilmer4[:]:
    eilmer4ArrT.append(temp[0])
    eilmer4ArrCp.append(temp[1]/1.987/CpConversion)

plt.plot(xplot, yplot)
plt.plot(eilmer4ArrT, eilmer4ArrCp)
plt.legend(['McBride Partition Func', 'Gas Prep Output'], loc=2)
plt.xlabel('Temperature (K)')
plt.ylabel('Cp (cal/K/mol)')
plt.title('Atomic Oxygen Specific Heat')
plt.show(block=True)