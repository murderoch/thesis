import spectra
import thermo
import util
import species
import plotter

import matplotlib.pyplot as plt

############### SORT LEVELS ##################


constants = util.Constants()

nMax = 20

speciesList = ['C', 'C+', 'N', 'N+', 'N++', 'O', 'O+', 'O++', 'F+', 'F++', 'F+++', 'Ne++', 'Ne+++']
useList = ['O', 'N', 'O+', 'N+']

for useName in useList:
    use = species.Species(useName)
    allSpecies = []
    for speciesStr in speciesList:
        speciesObj = species.Species(speciesStr)
        allSpecies.append([speciesObj, spectra.readNISTSpectra(speciesObj)])

    NIST = spectra.readNISTSpectra(use)
    theory = spectra.calculateExpectedStates(use, nMax)

    calcEnergy = spectra.CalcEnergy(nMax, use, NIST, theory, allSpecies)
    completeLevels, calculatedLevels = calcEnergy.populateTheory()

    calculatedLevels = spectra.sortSpectra(calculatedLevels, None)
    completeLevels = spectra.sortSpectra(completeLevels, None)
    NISTSorted = spectra.sortSpectra(NIST, None)

    tempRange = range(300, 50000, 50)
    ionizationLowering = 500
    #ionizationLowering = 'DebyeHuckel'
    numberDensity = 3E23

    thermoObj = thermo.Thermo(use, completeLevels, tempRange, ionizationLowering, numberDensity)

    #plotter.plotCp(use, tempRange, thermoObj, ionizationLowering, True, False)

    thermoProps = thermoObj.calcThermoPropsRange(tempRange)
    CpArr = thermoProps[0]
    HArr = thermoProps[1]
    SArr = thermoProps[2]
    
    fileName = useName + '-output.dat'
    with open(fileName, 'w') as exportFile:
        exportFile.write('T, Cp, H, S\n')
        for i, temp in enumerate(tempRange):
            string = str(temp) + ', ' + str(CpArr[i]) + ', ' + str(HArr[i]) + ', ' + str(SArr[i]) + '\n'
            exportFile.write(string)
