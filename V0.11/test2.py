from scipy import optimize
import matplotlib.pyplot as plt

Llist = [0, 0, 1, 1, 2, 2, 3, 3, 4, 4]
energyList = [105019.307, 105165.232, 105788.62733333332, 105912.031, 106751.47200000002, 106765.803, 106785.16, 106785.201, 106787.903, 106787.903]

IH = 109678.77177820474 
z = 8 
I = 109837.02186501028

def func(l, A):
    energy = I - IH*(z + 1.)**2. / (l + A)**2.
    return energy

p, e = optimize.curve_fit(func, Llist, energyList) 

print(p, e)
for l in range(0,6):
    plt.plot(l, func(l, *p), 'k.')

plt.plot(Llist, energyList, 'rx')
plt.show()