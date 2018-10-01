import spectra
import thermo
import util
import species

import matplotlib.pyplot as plt

############### SORT LEVELS ##################


constants = util.Constants()

nMax = 5

speciesList = ['O', 'O+', 'F+']

allSpecies = []
for speciesStr in speciesList:
    speciesObj = species.Species(speciesStr)
    allSpecies.append([speciesObj, spectra.readNISTSpectra(speciesObj)])

#print(allSpecies)

O = species.Species('O')
#He = species.Species('He')

use = O

NIST = spectra.readNISTSpectra(use)
theory = spectra.calculateExpectedStates(use, nMax)
calcEnergy = spectra.CalcEnergy(nMax, use, NIST, theory, allSpecies)
completeLevels = calcEnergy.populateTheory()

tempRange = range(300, 50000, 50)

asdf = thermo.Thermo(use, completeLevels, tempRange, 250)
CpRange = asdf.getCpRange()


with open('output.dat', 'w') as outputFile:
    for config in completeLevels:
        for level in config.term.levels:
            outputFile.writelines(config.ID + ', ' + str(level.J) + ', ' + str(level.energy) + '\n')


Capitelli2005oi = [[100, 200, 500, 700, 1000, 2000, 3000, 5000, 10000, 12000, 13000, 14000,
                15000, 16000, 17000, 18000, 19000, 20000, 22000, 23000, 24000, 25000, 26000, 27000, 28000,
                30000, 34000, 40000, 44000, 50000],
            [23.70, 22.74, 21.26, 21.04, 20.91, 20.83, 20.94, 21.80, 23.29, 24.35, 25.74, 28.15,
                31.89, 37.19, 44.04, 52.19, 61.11, 70.00, 84.37, 88.54, 90.36, 89.94, 87.66, 84.00, 79.46,
                69.37, 51.52, 35.92, 30.6, 26.30]]

Gordon1999 = [[100, 200, 500, 700, 1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000],
            [23.703, 22.734, 21.257, 21.040, 20.915, 20.827, 20.937, 21.799, 22.708, 23.150, 23.868, 24.722]]

plt.figure()
plt.plot(tempRange, CpRange, label = 'Mine')
plt.plot(Capitelli2005oi[0], Capitelli2005oi[1], label = 'Capitelli')
plt.plot(Gordon1999[0], Gordon1999[1], label = 'Gordon')
#plt.xlim([0, 10000])
plt.legend()
plt.show()
