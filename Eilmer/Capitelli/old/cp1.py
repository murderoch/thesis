import csv
from scipy import optimize
import matplotlib.pyplot as plt
import numpy as np

temps = []
speciesList = ['O', 'O+', 'O2', 'O2+', 'N', 'N+', 'N2', 'N2+', 'NO', 'NO+', 'e']
CpList = {species: [] for i,species in enumerate(speciesList)}



with open('thermoData.csv', 'r') as readFile:
    readData = csv.reader(readFile, delimiter=',')
    next(readData)
    next(readData)
    for row in readData:
        temps.append(float(row[0]))
        for i, species in enumerate(speciesList):
            CpList[species].append(float(row[i + 1]))


def piecewise_poly(x,
                   ka1, ka2, ka3, ka4, ka5, ka6, ka7,
                   kb1, kb2, kb3, kb4, kb5, kb6, kb7,
                   kc1, kc2, kc3, kc4, kc5, kc6, kc7,
                   kd1, kd2, kd3, kd4, kd5, kd6, kd7,
                   ke1, ke2, ke3, ke4, ke5, ke6, ke7):

    condlist = [x < ranges[0], x < ranges[1], x <= ranges[2], x <= ranges[3], x >= ranges[3]]
    funclist = [lambda x: ka1*x**-2. + ka2*x**-1. + ka3 + ka4*x + ka5*x**2. + ka6*x**3. + ka7*x**4.,
                lambda x: kb1*x**-2. + kb2*x**-1. + kb3 + kb4*x + kb5*x**2. + kb6*x**3. + kb7*x**4.,
                lambda x: kc1*x**-2. + kc2*x**-1. + kc3 + kc4*x + kc5*x**2. + kc6*x**3. + kc7*x**4.,
                lambda x: kd1*x**-2. + kd2*x**-1. + kd3 + kd4*x + kd5*x**2. + kd6*x**3. + kd7*x**4.,
                lambda x: ke1*x**-2. + ke2*x**-1. + ke3 + ke4*x + ke5*x**2. + ke6*x**3. + ke7*x**4.]
    
    return np.piecewise(x, condlist, funclist)


species = 'O'

#divide = 12000
ranges = [6000, 12000, 15000, 22000]
#tempIdx = temps.index(divide)


test = CpList[species]
#test1 = test[:tempIdx]
#test2 = test[tempIdx:]
#temps1 = temps[:tempIdx]
#temps2 = temps[tempIdx:]

#global trial

#trial = 1
p, e = optimize.curve_fit(piecewise_poly, temps, test)
#trial = 2
#p2, e2 = optimize.curve_fit(piecewise_poly, temps2, test2)

xd = np.linspace(min(temps), max(temps), 1000)
#xd2 = np.linspace(min(temps2)-1000, max(temps2), 1000)

plt.figure()
plt.plot(temps, test, 'kx')
plt.plot(xd, piecewise_poly(xd, *p))
#plt.plot(xd2, piecewise_poly(xd2, *p2))
plt.show()