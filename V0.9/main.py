import spectra
import thermo
import util
import species
import plotter

import matplotlib.pyplot as plt

############### SORT LEVELS ##################


constants = util.Constants()

nMax = 30

speciesList = ['C', 'C+', 'N', 'N+', 'N++', 'O', 'O+', 'O++', 'F+', 'F++', 'F+++', 'Ne++', 'Ne+++']

allSpecies = []
for speciesStr in speciesList:
    speciesObj = species.Species(speciesStr)
    allSpecies.append([speciesObj, spectra.readNISTSpectra(speciesObj)])

O = species.Species('O')
N = species.Species('N')
OI = species.Species('O+')
NI = species.Species('N+')

use = NI

NIST = spectra.readNISTSpectra(use)
theory = spectra.calculateExpectedStates(use, nMax)

calcEnergy = spectra.CalcEnergy(nMax, use, NIST, theory, allSpecies)
completeLevels, calculatedLevels = calcEnergy.populateTheory()

calculatedLevels = spectra.sortSpectra(calculatedLevels, 'calculated')
completeLevels = spectra.sortSpectra(completeLevels, 'complete')
NISTSorted = spectra.sortSpectra(NIST, 'NIST')

tempRange = range(300, 50000, 50)
ionizationLowering = 250

asdf = thermo.Thermo(use, completeLevels, tempRange, ionizationLowering, 3E23)

plotter.plotCp(use, tempRange, asdf, ionizationLowering, True, False)
