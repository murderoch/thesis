import spectra
import thermo
import util
import species

constants = util.Constants()

#species = species.Species('He')
nMax = 20

asdf = species.hi()

he = []

for n in range(species.core.maxN+1, nMax+1, 1):
    for L in range(0, nMax, 1):
        excitedShell = spectra.SubShell(n, L, 1)
        excitedState = spectra.ExcitedState(excitedShell)
        configuration = spectra.Configuration(species.core, excitedState, n, species)
        he.append(configuration)

