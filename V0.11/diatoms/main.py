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
        
        for n in range(0, 7):#len(self.species.spectralData['Te'])):
            Eelec, ge = self.Eelec(n)
            vlim = self.vibrational(n)
            QelecSum += ge * exp(-Eelec * self.constants.Cm_1ToJoules / (self.constants.kB * T))

            for v in range(0, vlim):
                Evib = self.Evib(n, v)
                QvibSum += exp(-Evib * self.constants.Cm_1ToJoules / (self.constants.kB * T))
                
                jlim= self.matchVrmPair(n, v, Evib)
                #print(n, v, max(jlim))
                for J in jlim:
                    Erot = self.Erot(n, v, J)

                    Esum = Erot + Evib + Eelec
                    #print(n, v, J, Erot)
                    QrotSum += (2*J + 1) * exp(-Erot / self.constants.Cm_1ToJoules / (self.constants.kB * T))     
                #print(QelecSum, QvibSum, QrotSum)
        #QvibSum = 1.
        #print(QelecSum, QvibSum, QrotSum)
        
        #QelecSum = 1
        #QvibSum = 1
        QrotSum = 1

        Q = 1./sigma * QelecSum * QvibSum * QrotSum
        print(Q)
        return Q

    def Eelec(self, n):
        energy = self.species.spectralData['Te'][n]
        ge = self.species.spectralData['gi'][n]
        return energy, ge
    
    def Evib(self, n, v):
        spectralData = self.species.spectralData
        #Ediss = spectralData['Ediss'][n]
        omegaE = [spectralData['we'][n], spectralData['wexe'][n], spectralData['weye'][n], spectralData['weze'][n], spectralData['weke'][n]]

        omega0 = [omegaE[0] - omegaE[1] + 3./4.*omegaE[2] + 1./8.*omegaE[3] + 3./16.*omegaE[4],
                  omegaE[1] - 3./2.*omegaE[2] - 3./2.*omegaE[3] - 5./4.*omegaE[4],
                  omegaE[2] + 2.*omegaE[3] + 5./2.*omegaE[4],
                  omegaE[3] + 5./2.*omegaE[4],
                  omegaE[4]
                 ]

        vibrationalEnergy_0 = 1./2.*omegaE[0] + 1./4.*omegaE[1] + 1./8.*omegaE[2] + 1./16.*omegaE[3] + 1./32.*omegaE[4]
        vibrationalEnergy = vibrationalEnergy_0 + omega0[0]*v + omega0[1]*v**2. + omega0[2]*v**3. + omega0[3]*v**4. + omega0[4]*v**5.
        return vibrationalEnergy

    def Erot(self, n, v, J):
        spectralData = self.species.spectralData  

        Be = spectralData['Be'][n]                  #cm^-1
        alphaE = spectralData['alphaE'][n]          #cm^-1
        gammaE = spectralData['gammaE'][n]          #cm^-1
        deltaE = spectralData['deltaE'][n]          #cm^-1
        etaE = spectralData['etaE'][n]              #cm^-1
    
        De = spectralData['De'][n]                  #cm^-1
        betaE = spectralData['betaE'][n]            #cm^-1
        gE = spectralData['gE'][n]                  #cm^-1

        Hv0 = spectralData['Hv0'][n]                #cm^-1
        Hv1 = spectralData['Hv1'][n]                #cm^-1

        #Dunham series expansion
        Bv = Be - alphaE*(v + 1./2.) + gammaE*(v + 1./2.)**2. + deltaE*(v + 1./2.)**3. + etaE*(v + 1./2.)**4.
        Dv = De - betaE*(v + 1./2.) + gE*(v + 1./2.)**2.
        Hv = Hv0 - Hv1*(v + 1./2.)

        energyRot = Bv*J*(J + 1) - Dv*J**2.*(J + 1)**2. + Hv*J**3.*(J + 1)**3.

        return energyRot


    def vibrational(self, n):
        spectralData = self.species.spectralData
        Ediss = spectralData['Ediss'][n]# * self.constants.Cm_1ToJoules
        omegaE = [spectralData['we'][n], spectralData['wexe'][n], spectralData['weye'][n], spectralData['weze'][n], spectralData['weke'][n]]

        omega0 = [omegaE[0] - omegaE[1] + 3./4.*omegaE[2] + 1./8.*omegaE[3] + 3./16.*omegaE[4],
                  omegaE[1] - 3./2.*omegaE[2] - 3./2.*omegaE[3] - 5./4.*omegaE[4],
                  omegaE[2] + 2.*omegaE[3] + 5./2.*omegaE[4],
                  omegaE[3] + 5./2.*omegaE[4],
                  omegaE[4]
                 ]

        #This is weird, could well be wrong to use Capitelli's values (they're negative)
        Ediss = abs(Ediss)

        #negative order as root function assigns largest power first
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
        Eelec = spectralData['Te'][n]
        re = spectralData['re'][n] * 1E-8
        mu = self.species.reducedMass

        Jmax = 1
        while self.rotational(n, Jmax, omegaE) != False:
            Jmax += 1
        if Jmax == 0:
            return 0

        
        JenergyMax = 1
        run = True
        while run:
            #rm, Urm = self.rotational(n, JenergyMax, omegaE)
            #print(JenergyMax)
            #minEnergy = (self.constants.h*100)*(JenergyMax*(JenergyMax + 1.))/(8*pi**2.*self.constants.c*mu*re**2.)
            #print(minEnergy)
            Erot = self.Erot(n, v, JenergyMax)
            #print(JenergyMax, Erot + Evib, Ediss)
            if Erot + Evib < Ediss:
                JenergyMax += 1
            else:
                run = False
        
        #print('Jmax =', Jmax, 'with Urm =', self.rotational(n, Jmax-1, omegaE)[1])
        #print('JmaxEnergy =', JenergyMax, 'with Urm =', self.rotational(n, JenergyMax, omegaE)[1], '<', Ediss)

        #print(JenergyMax)
        Jlist = range(1, Jmax)
        UrmList = []
        vibRotCoupleList = []
        diffList = []
        for J in Jlist:
            rm, Urm = self.rotational(n, J, omegaE)
            vibRotCouple = self.Evib(n, v) + self.Erot(n, v, J)

            diff = vibRotCouple - Urm

            UrmList.append(Urm)
            vibRotCoupleList.append(vibRotCouple)
            diffList.append(diff)

        #print(Jlist[UrmList.index(min(Urm))])
        #print(UrmList[320], UrmList[321], UrmList[322])
        #print(len(UrmList))
        '''
        UrmList[294] = 6.14E4
        UrmList[295] = 6.16E4
        plt.plot(Jlist, UrmList, label='$U^{rot}_m(J)$')
        plt.plot(Jlist, vibRotCoupleList, label='$\\varepsilon_m(J)$')
        plt.grid()
        plt.title('$N_2$ Vib-Rotational Energy and Potential')
        plt.xlabel('Rotational Quantum Number J')
        plt.ylabel('Energy $cm^{-1}$')
        plt.xlim([0,320])
        plt.ylim([0, 1E5])
        plt.ticklabel_format(style='sci', axis='y', scilimits=(1,4))
        plt.legend()
        '''

        minSep = min([abs(energy) for energy in diffList])

        try:
            Jquasi = diffList.index(minSep)
            #print('v=', v, 'Jlim=', Jlim)
            plt.show()

        except ValueError:
            Jquasi = diffList.index(-minSep)
            #print('v=', v, 'Jlim=', Jlim)
            plt.show()

        Jlim = []
        Jlim.extend(range(0, JenergyMax))

        if Jquasi > max(Jlim):
            Jlim.append(Jquasi)
            
        return Jlim     

    def rotational(self, n, J, omegaE):
        mu = self.species.reducedMass               #kg
        spectralData = self.species.spectralData    
        re = spectralData['re'][n] * 1E-8           #cm
        Ediss = abs(spectralData['Ediss'][n])       #cm^-1

        beta = omegaE*sqrt(2*pi**2.*self.constants.c*mu / (Ediss * 100.* self.constants.h))

        def potential(rm, J):
            U_0 = Ediss * (1. - exp(-beta*(rm - re)))**2.
            U = U_0 + (self.constants.h*100)*(J*(J + 1.))/(8*pi**2.*self.constants.c*mu*rm**2.)
            return U
        
        def potentialDerivative(rm, J):
            first = 2.*Ediss*beta*exp(-beta*(rm - re))*(1. - exp(-beta*(rm- re)))
            second = (2./(rm**3.)) * ( self.constants.h*100.*J*(J + 1.) ) / (8*pi**2.*mu*self.constants.c)
            dU_dr = first - second
            return dU_dr
        
        if 0:
            JList = [0, 100, 200, 300, 400]
            #JList = [0, 100, 350, 400]
            #JList = [323]

            rmX = np.linspace(0.6E-8, 3E-8, num=1000)
            xPlot = []
            UPlot = [[] for J in JList]
            dU_drPlot = [[] for J in JList]

            for rm in rmX:
                for i, Jplt in enumerate(JList):
                    UPlot[i].append(potential(rm, Jplt))
                    dU_drPlot[i].append(potentialDerivative(rm, Jplt))
                xPlot.append(rm)                 
            
            plt.figure()
            for i, Jplt in enumerate(JList):
                label = 'J = ' + str(Jplt)
                plt.plot(xPlot, UPlot[i], label=label)
            plt.title('Molecular Potential Well of N2')
            plt.ylim([0, 3.5E5])
            plt.xlim([0.6E-8, 3.0E-8])
            plt.ylabel('Molecular Potential (cm$^{-1}$)')
            plt.xlabel('Internuclear Distance (cm)')
            plt.ticklabel_format(style='sci', axis='y', scilimits=(1,4))
            plt.legend()      
            plt.grid()     
                 
            plt.figure()
            for i, Jplt in enumerate(JList):
                label = 'J = ' + str(Jplt)
                plt.plot(xPlot, dU_drPlot[i], label=label)
            plt.title('Derivative of the Molecular Potential of N2')
            plt.xlabel('Internuclear Distance (cm)')
            plt.ylabel('Molecular Potential Derivative (cm$^{-2}$)')
            #plt.ylim([-1.2E10, 0.4E10])
            plt.ylim([-1E13, 1E13])
            plt.xlim([0.6E-8, 3.0E-8])
            plt.legend()
            plt.grid()
            plt.show()

        '''
        rm = fsolve(potentialDerivative, re, args=J)
        Urm = potential(rm, J)
        print(J, rm, Urm)
        '''


        #rm = fsolve(potentialDerivative, re*1.7, args=J)
        #Urm = potential(rm, J)
        #print(J, rm, Urm)


        
        warnings.simplefilter('error')
        try:
            rm = fsolve(potentialDerivative, re*1.7, args=J)
            Urm = potential(rm, J)        
            #print(J, rm, Urm)
            return rm, Urm
        except RuntimeWarning:
            return False


N2 = species.Species('N2')
O2 = species.Species('O2')
#print(N2.spectralData['Te'])

energyModes = EnergyModes(O2)
#asdf = energyModes.vibrational(0)

Q = energyModes.getEnergy(300)

print(Q)

'''
n = 0
v = 3
J = 1

omegaE = 1580.19

eVib = energyModes.Evib(n, v)
energyModes.vibrational(n)

energyModes.rotational(n, J, omegaE)
'''

#print(energyModes.matchVrmPair(n, v, eVib))
#print(energyModes.rotational(n, 10, omegaE))

#print(Q)