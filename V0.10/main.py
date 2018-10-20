import spectra
import thermo
import util
import species
import plotter

import matplotlib.pyplot as plt

############### SORT LEVELS ##################


constants = util.Constants()

nMax = 10

speciesList = ['C', 'C+', 'N', 'N+', 'N++', 'O', 'O+', 'O++', 'F+', 'F++', 'F+++', 'Ne++', 'Ne+++']

allSpecies = []
for speciesStr in speciesList:
    speciesObj = species.Species(speciesStr)
    allSpecies.append([speciesObj, spectra.readNISTSpectra(speciesObj)])

O = species.Species('O')
N = species.Species('N')
OI = species.Species('O+')
NI = species.Species('N+')

use = OI

NIST = spectra.readNISTSpectra(use)
theory = spectra.calculateExpectedStates(use, nMax)

calcEnergy = spectra.CalcEnergy(nMax, use, NIST, theory, allSpecies)
completeLevels, calculatedLevels = calcEnergy.populateTheory()

calculatedLevels = spectra.sortSpectra(calculatedLevels, None)
completeLevels = spectra.sortSpectra(completeLevels, None)
NISTSorted = spectra.sortSpectra(NIST, None)

tempRange = range(300, 50000, 50)
ionizationLowering = 250
#ionizationLowering = 'DebyeHuckel'
numberDensity = 3E23

asdf = thermo.Thermo(use, completeLevels, tempRange, ionizationLowering, numberDensity)

plotter.plotCp(use, tempRange, asdf, ionizationLowering, True, False)
