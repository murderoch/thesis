import matplotlib.pyplot as plt

R = 8.3144598

def calcMcBride2002(temp):
    if temp < 1000:
        Cp = -7.9536113E3*temp**-2 + 1.607177787E2*temp**-1 + 1.966226438 + 1.01367031E-03*temp \
             - 1.110415423E-6*temp**2 + 6.5175065E-10*temp**3 - 1.58477925E-13*temp**4
    elif temp < 6000:
        Cp = 2.61902026E5*temp**-2 + -7.29872203E2*temp**-1 + 3.31717727 + -4.28133436E-4*temp \
             + 1.036104594E-7*temp**2 - 9.43830433E-12*temp**3 + 2.725038297E-16*temp**4
    else:
        Cp = 1.77900426E8*temp**-2 - 1.082328257E5*temp**-1 + 2.81077836E1 - 2.97523226E-3*temp \
             + 1.853997534E-7*temp**2 - 5.7623154E-12*temp**3 + 7.191720164E-17*temp**4   
    return Cp * R


def calcGupta1990(temp):
    if temp < 1000:
        Cp = 0.28236E1 -0.89478E-3*temp + 0.83080E-6*temp**2 - 0.16837E-9*temp**3 \
             - 0.73205E-13*temp**4
    elif temp < 6000:
        Cp = 0.25431E1 -0.27551E-4*temp - 0.31028E-8*temp**2 + 0.45511E-11*temp**3 \
             - 0.43681E-15*temp**4
    elif temp < 15000:
        Cp = 0.25460E1 -0.59520E-4*temp + 0.27010E-7*temp**2 - 0.27980E-11*temp**3 \
             + 0.93800E-16*temp**4
    elif temp < 25000:
        Cp = -0.97871E-2 + 0.12450E-2*temp - 0.16154E-6*temp**2 + 0.80380E-11*temp**3 \
             - 0.12624E-15*temp**4
    else:
        Cp = 0.16428E2 - 0.39313E-2*temp + 0.2984E-6*temp**2 - 0.81613E-11*temp**3 \
             + 0.75004E-16*temp**4
    return Cp * R


Gordon1999 = [[100, 200, 500, 700, 1000, 2000, 3000, 5000, 10000, 15000, 20000],
               [23.703, 22.734, 21.257, 21.040, 20.915, 20.827, 20.937, 21.799, 23.150, 23.868, 24.722]]

McBride2002 = [range(300,20000,50),[]]
Gupta1990 = [range(300,20000,50),[]]

for temp in McBride2002[0]:
    McBride2002[1].append(calcMcBride2002(temp))

for temp in Gupta1990[0]:
    Gupta1990[1].append(calcGupta1990(temp))

plt.figure()
plt.title('Comparison of Thermodynamic Sources for Atomic Oxygen')
plt.plot(Gordon1999[0], Gordon1999[1], 'b')
plt.plot(McBride2002[0], McBride2002[1], 'r')
plt.plot(Gupta1990[0], Gupta1990[1], 'g')
plt.ylim([0, 38])
plt.legend(['Gordon 1999', 'McBride 2002', 'Gupta 1990'], loc=4)
plt.xlabel('Temperature (K)')
plt.ylabel('Cp (cal / g mol K)')
plt.grid()
plt.show()
