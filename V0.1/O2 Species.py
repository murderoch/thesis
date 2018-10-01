import matplotlib.pyplot as plt

R = 8.3144598

def calcMcBride2002(temp):
    if temp < 1000:
        Cp = -3.4255634E4*temp**-2 + 4.84700097E2*temp**-1 + 1.11901096 + 4.29388924E-03*temp \
             - 6.83630052E-7*temp**2 + -2.0233727E-9*temp**3 + 1.039040018E-12*temp**4
    elif temp < 6000:
        Cp = -1.037939022E6*temp**-2 + 2.344830282E3*temp**-1 + 1.819732036 + 1.267847582E-3*temp \
             -2.188067988E-7*temp**2 + 2.0537195272E-11*temp**3 - 8.19346705E-16*temp**4
    else:
        Cp = 4.9752943E8*temp**-2 - 2.866106874E5*temp**-1 + 6.69035225E1 - 6.169959020E-3*temp \
             + 3.016396027E-7*temp**2 - 7.4214166E-12*temp**3 + 7.27817577E-17*temp**4   
    return Cp * R


def calcGupta1990(temp):
    if temp < 1000:
        Cp = 0.36146E1 - 0.18598E-2*temp + 0.70814E-5*temp**2 - 0.68080E-8*temp**3 \
             + 0.21628E-11*temp**4
    elif temp < 6000:
        Cp = 0.35949E1 + 0.75213E-3*temp - 0.18732E-6*temp**2 + 0.27913E-10*temp**3 \
             - 0.15774E-14*temp**4
    elif temp < 15000:
        Cp = 0.38599E1 + 0.32510E-3*temp - 0.92131E-8*temp**2 - 0.78684E-12*temp**3 \
             + 0.29426E-16*temp**4
    elif temp < 25000:
        Cp = 0.34867E1 + 0.52384E-3*temp - 0.39123E-7*temp**2 + 0.10094E-11*temp**3 \
             -0.88718E-17*temp**4
    else:
        Cp = 0.39620E1 + 0.39446E-3*temp - 0.29506E-7*temp**2 + 0.73975E-12*temp**3 \
             -0.64209E-17*temp**4
    return Cp * R



Jaffe1987 = [[300, 600, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000,
              20000, 25000, 30000, 35000, 40000 ,45000, 50000],
             [1.0344, 1.3597, 1.6949, 2.0457, 2.3127, 2.5181, 2.6741, 2.7955, 2.8571, 2.8329, 2.7247,
              2.5565, 1.6084, 1.0531, 0.7630, 0.5914, 0.4774, 0.3954, 0.3335, 0.2853]]

Capitelli2005 = [[300, 600, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000,
              20000, 25000, 30000, 34000, 40000 ,44000, 50000],
             [29.39, 32.09, 34.88, 37.80, 40.02, 41.72, 43.02, 44.00, 44.45, 44.18, 43.20,
              41.74, 33.71, 29.07, 26.69, 25.31, 24.56, 23.78, 23.39, 22.94]]

for i, Cp in enumerate(Jaffe1987[1]):
    Jaffe1987[1][i] = Jaffe1987[1][i] * R + 5/2*R

McBride2002 = [range(300,20000,50),[]]
Gupta1990 = [range(300,30000,50),[]]


for temp in McBride2002[0]:
    McBride2002[1].append(calcMcBride2002(temp))

for temp in Gupta1990[0]:
    Gupta1990[1].append(calcGupta1990(temp))

plt.figure()
plt.title('Comparison of Thermodynamic Sources for Molecular Oxygen')
plt.plot(Jaffe1987[0], Jaffe1987[1], 'b')
plt.plot(Gupta1990[0], Gupta1990[1], 'g')
plt.plot(McBride2002[0], McBride2002[1], 'r')
plt.plot(Capitelli2005[0], Capitelli2005[1], 'violet')
plt.ylim([0, 50])
plt.legend(['Jaffe 1987', 'Gupta 1990',  'McBride 2002', 'Capitelli 2005'], loc=4)
plt.xlabel('Temperature (K)')
plt.ylabel('Cp (cal / g mol K)')
plt.grid()
plt.show()
