import spectra
import thermo
import util
import species

import matplotlib.pyplot as plt

############### SORT LEVELS ##################


constants = util.Constants()

nMax = 30

speciesList = ['O', 'O+', 'F+', 'Ne++']

allSpecies = []
for speciesStr in speciesList:
    speciesObj = species.Species(speciesStr)
    allSpecies.append([speciesObj, spectra.readNISTSpectra(speciesObj)])

O = species.Species('O')

use = O

NIST = spectra.readNISTSpectra(use)
theory = spectra.calculateExpectedStates(use, nMax)

calcEnergy = spectra.CalcEnergy(nMax, use, NIST, theory, allSpecies)
completeLevels, calculatedLevels = calcEnergy.populateTheory()

for config in calculatedLevels:
    if config.term.checkLevelEnergies() != None:
        print('Missing energy for:', config.ID, config.core.term.getTermString())

calculatedLevels.sort(key=lambda x: x.term.maxEnergy)
with open('calculatedLevels.dat', 'w') as outputFile:
    for config in calculatedLevels:
        outputFile.write('\n')
        for level in config.term.levels:
            outputFile.writelines(config.ID + ', ' + str(level.J) + ', ' + str(level.energy) + '\n')



for config in completeLevels:
    if config.term.checkLevelEnergies() != None:
        print('Missing energy for:', config.ID, config.core.term.getTermString())

with open('output.dat2', 'w') as outputFile:
    for config in completeLevels:
        outputFile.write('\n')
        for level in config.term.levels:
            outputFile.writelines(config.ID + ', ' + str(level.J) + ', ' + str(level.energy) + '\n')

completeLevels.sort(key=lambda x: x.term.maxEnergy)

tempRange = range(300, 50000, 50)

asdf = thermo.Thermo(use, completeLevels, tempRange, 250)
CpRange = asdf.getCpRange()


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
