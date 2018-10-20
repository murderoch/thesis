import util
import species

import numpy as np
from math import sqrt, exp, pi
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import warnings


class EnergyModes:
    def __init__(self, species):
        self.constants = util.Constants()
        self.species = species

    def getEnergy(self, T):
        
        sigma = self.species.symmetryFactor

        QelecSum = 0
        QvibSum = 0
        QrotSum = 0

        for n in range(0,1):# len(self.species.spectralData['Te'])):
            Eelec, ge = self.Eelec(n)
            vlim = self.vibrational(n)
            QelecSum += ge * exp(- Eelec*self.constants.Cm_1ToJoules /(self.constants.kB * T))

            for v in range(0, vlim):
                Evib = self.Evib(n, v)
                QvibSum += exp(-Evib*self.constants.Cm_1ToJoules / (self.constants.kB * T))
                
                jlim = self.matchVrmPair(n, v, Evib)
                print(v, jlim)
                for J in range(0, jlim):
                    Erot = self.Erot(n, v, J)
                    QrotSum += (2*J + 1) * exp(-Erot*self.constants.Cm_1ToJoules / (self.constants.kB * T))     

        Q = sigma * QelecSum * QvibSum * QrotSum
        
        return Q

    def Eelec(self, n):
        energy = self.species.spectralData['Te'][n]
        ge = self.species.spectralData['gi'][n]
        return energy, ge
    
    def Evib(self, n, v):
        spectralData = self.species.spectralData
        Ediss = spectralData['Ediss'][n]# * self.constants.Cm_1ToJoules
        omegaE = [spectralData['we'][n], spectralData['wexe'][n], spectralData['weye'][n], spectralData['weze'][n], spectralData['weke'][n]]

        omega0 = [omegaE[0] - omegaE[1] + 3./4.*omegaE[2] + 1./8.*omegaE[3] + 3./16.*omegaE[4],
                  omegaE[1] - 3./2.*omegaE[2] - 3./2.*omegaE[3] - 5./4.*omegaE[4],
                  omegaE[2] + 2.*omegaE[3] + 5./2.*omegaE[4],
                  omegaE[3] + 5./4.*omegaE[4],
                  omegaE[4]
                 ]

        vibrationalEnergy_0 = 1./2.*omegaE[0] + 1./4.*omegaE[1] + 1./8.*omegaE[2] + 1./16.*omegaE[3] + 1./32.*omegaE[4]
        return vibrationalEnergy_0 + omega0[0]*v - omega0[1]*v**2. + omega0[2]*v**3. + omega0[3]*v**4. + omega0[4]*v**5.
    

    def Erot(self, n, v, J):
        spectralData = self.species.spectralData  
        Be = spectralData['Be'][n]                  #cm^-1
        De = spectralData['De'][n]                  #cm^-1
        alphaE = spectralData['alphaE'][n]          
        betaE = spectralData['betaE'][n]

        Bv = Be - alphaE*(v + 1./2.)
        Dv = De - betaE*(v + 1./2.)
        energyRot = Bv*J*(J + 1) - Dv*J**2.*(J + 1)**2.

        return energyRot


    def vibrational(self, n):
        spectralData = self.species.spectralData
        Ediss = spectralData['Ediss'][n]# * self.constants.Cm_1ToJoules
        omegaE = [spectralData['we'][n], spectralData['wexe'][n], spectralData['weye'][n], spectralData['weze'][n], spectralData['weke'][n]]

        omega0 = [omegaE[0] - omegaE[1] + 3./4.*omegaE[2] + 1./8.*omegaE[3] + 3./16.*omegaE[4],
                  omegaE[1] - 3./2.*omegaE[2] - 3./2.*omegaE[3] - 5./4.*omegaE[4],
                  omegaE[2] + 2.*omegaE[3] + 5./2.*omegaE[4],
                  omegaE[3] + 5./4.*omegaE[4],
                  omegaE[4]
                 ]

        #This is weird, could well be wrong to use Capitelli's values (they're negative)
        Ediss = abs(Ediss)
        coeffs = [omega0[3], omega0[2], -omega0[1], omega0[0], -Ediss]
        vLim = np.roots(coeffs)
        vLim = vLim.real[(abs(vLim.imag)<1e-5) & (vLim.real >= 0)]

        if len(vLim) == 0:
            vLim = 0
        else:
            vLim = int(round(min(vLim))) + 1
        return vLim


    def matchVrmPair(self, n, v, Evib):
        spectralData = self.species.spectralData  
        omegaE = spectralData['we'][n]
        Ediss = abs(spectralData['Ediss'][n])

        Jmax = 0
        while self.rotational(n, Jmax, omegaE) != False:
            Jmax += 1
        if Jmax == 0:
            return 0
        

        JenergyMax = 0
        run = True
        while run:
            rm, Urm = self.rotational(n, JenergyMax, omegaE)
            if Urm < Ediss:
                JenergyMax += 1
            else:
                run = False

        JList = range(0, Jmax)
        energyList = []
        UrotList = []

        for J in JList:
            rm, Urot = self.rotational(n, J, omegaE)

            UrotList.append(Urot)

            diff = Evib + self.Erot(n, v, J) - Urot
            energyList.append(diff)
        
        rm, U = self.rotational(0, 168, omegaE)
        print(self.Evib(0, 10) + self.Erot(0, 10, 168), U)


        plt.plot(JList, energyList)
        plt.plot(JList, UrotList)
        plt.show()

        minSep = min([abs(energy) for energy in energyList])

        try:
            Jlim = energyList.index(minSep)
            print(Jlim)
        except ValueError:
            Jlim = energyList.index(-minSep)
            print(Jlim)
        return Jlim


    def rotational(self, n, J, omegaE):
        mu = self.species.reducedMass               #kg
        spectralData = self.species.spectralData    
        re = spectralData['re'][n] * 1E-8           #cm
        Ediss = abs(spectralData['Ediss'][n])       #cm^-1

        beta = omegaE*sqrt(2*pi**2.*self.constants.c*mu / (Ediss * 100.* self.constants.h))

        def potential(rm, J):
            U_0 = Ediss * (1. - exp(-beta*(rm - re)))**2.
            U = U_0 + ((self.constants.h*100))*(J*(J + 1.))/(8*pi**2.*self.constants.c*mu*rm**2.)
            return U
        
        def potentialDerivative(rm, J):
            #print(rm)
            first = 2.*Ediss*beta*exp(-beta*(rm - re))*(1. - exp(-beta*(rm- re)))
            second = (2./(rm**3.)) * ( self.constants.h*100.*J*(J + 1.) ) / (8*pi**2.*mu*self.constants.c)
            dU_dr = first - second
            return dU_dr
        
        if 0:
            #JList = [0, 100, 150, 200, 250]
            JList = [0, 0]
            rmX = np.linspace(0.6E-8, 3E-8, num=100)
            xPlot = []
            UPlot = [[] for J in JList]
            dU_drPlot = [[] for J in JList]

            for rm in rmX:
                for i, J in enumerate(JList):
                    UPlot[i].append(potential(rm, J))
                    dU_drPlot[i].append(potentialDerivative(rm, J))
                xPlot.append(rm)      

            dX = rmX[1] - rmX[0]
            dU_drPlot[1] = np.gradient(UPlot[0], dX)
            

            '''
            plt.figure()
            for i, J in enumerate(JList):
                label = 'J = ' + str(J)
                plt.plot(xPlot, UPlot[i], label=label)
            plt.title('Molecular Potential')
            plt.ylim([0, 3.5E5])
            plt.xlim([0.6E-8, 3.0E-8])
            plt.legend()
            '''
                    
            plt.figure()
            for i, J in enumerate(JList):
                label = 'J = ' + str(J)
                plt.plot(xPlot, dU_drPlot[i], label=label)
            plt.title('Molecular Potential Derivative')
            #plt.ylim([-1.2E10, 0.4E10])
            #plt.ylim([-1E13, 3E13])
            plt.xlim([0.6E-8, 3.0E-8])
            plt.legend()
            plt.show()

        warnings.simplefilter('error')
        try:
            rm = fsolve(potentialDerivative, re*1.1, args=J)
            Urm = potential(rm, J)        
            return rm, Urm
        except RuntimeWarning:
            return False
        



N2 = species.Species('O2')
#print(N2.spectralData['Te'])

energyModes = EnergyModes(N2)
#asdf = energyModes.vibrational(0)

#Q = energyModes.getEnergy(300)

n = 0
v = 10

eVib = energyModes.Evib(n, v)
print(energyModes.matchVrmPair(n, v, eVib))
#print(Q)