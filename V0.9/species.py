import json

import util

class Species:
    def __init__(self, species):
        self.name = species
        self.constants = util.Constants()
        self.readProperties()

    def readProperties(self):
        jsonFilePath = util.getSpectralDataDir() + '/atomData.json'
        with open(jsonFilePath, 'r') as jsonFile:
            jsonData = json.load(jsonFile)

            speciesProperties = jsonData[self.name]
            self.atomicNumber = speciesProperties['atomicNumber']
            self.charge = speciesProperties['charge']
            self.noElectrons = speciesProperties['noElectrons']
            self.Io = speciesProperties['Io'] * self.constants.eVToCm_1
            self.cores = []
            for core in speciesProperties['cores'].items():
                self.cores.append(core[1])
    

#He = Species('He')
