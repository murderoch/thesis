import util

import matplotlib.pyplot as plt
from math import pi, exp, sqrt

constants = util.Constants()

Io = 13.618106 * constants.eVToCm_1

z = 0
kB = constants.kB
epsilon0 = constants.vacPermiativity
qE = constants.electronCharge
nE = 0.02504E27

nE = 3E23

nE = 2.5E19

rhoAtm = 1.225 * 0.001
nE = rhoAtm * constants.avagadro / 16

rhoList = [rhoAtm*0.5, rhoAtm, rhoAtm*10, rhoAtm*30, rhoAtm*100]
nameList = ['$0.5 Atm$', '$1 Atm$', '$10 Atm$', '$30 Atm$', '$100 Atm$']


def getCutoff(T, nE):

    '''
    lambdaD = sqrt(epsilon0 * kB * T / (qE**2. * nE))
    loweringEnergy = qE**2. *  (z + 1)/(4*pi * epsilon0 * lambdaD)

    loweringCm = loweringEnergy/constants.Cm_1ToJoules
    cutoffCm = Io - loweringCm
    '''

    lambdaD = sqrt(epsilon0 * kB * T / (qE**2. * nE*1**2))
    kDH = qE**2./(8*pi*epsilon0*lambdaD)
    loweringEnergy = 2*kDH*(z+1)


    loweringCm = loweringEnergy/constants.Cm_1ToJoules
    cutoffCm = Io - loweringCm

    return loweringCm

print(getCutoff(25000, 1E23))

tempRange = range(10, 50000, 10)

for i, rho in enumerate(rhoList):
    nE = rho * constants.avagadro / 16
    plotY = [getCutoff(T, nE) for T in tempRange]
    plt.plot(tempRange, plotY, label=nameList[i])

plt.grid()
plt.xlabel('Temperature $(K)$')
plt.ylabel('Ionization Lowering $cm^{-1}$')
plt.xlim([0, 50000])
plt.ylim([0, 50000])
plt.legend()
plt.show()
