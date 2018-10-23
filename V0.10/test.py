import util

import matplotlib.pyplot as plt
from math import pi, exp, sqrt

constants = util.Constants()

Io = 13.618106 * constants.eVToCm_1

z = 0
kB = constants.kB
epsilon0 = constants.vacPermiativity
qE = constants.electronCharge
#nE = 0.02504E27

nE = 3E23

def getCutoff(T):
    lambdaD = sqrt(epsilon0 * kB * T / (qE**2. * nE))

    loweringEnergy = qE**2. *  (z + 1)/(4*pi * epsilon0 * lambdaD)
    
    loweringCm = loweringEnergy#/constants.Cm_1ToJoules
    
    cutoffCm = Io - loweringCm

    return loweringCm

print(getCutoff(30000))

tempRange = range(200, 50000, 10)

plotY = [getCutoff(T) for T in tempRange]

plt.figure()
plt.plot(tempRange, plotY, 'k.')
plt.show()
