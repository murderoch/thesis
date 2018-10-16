from scipy import optimize
import matplotlib.pyplot as plt
import numpy as np


def func(x, n, y):
    energy = I - (IH*(z + 1.)**2.) / (n + x[0] + x[1] / (n**2.))**2. - y
    return energy

def func2(x, n, y):
    energy = x[0]*n**-1. + x[1] + x[2]*n + x[3]*n**2. - y   
    return energy


#I = 136657.9932077545
IH = 109678.8162193522
z = 0

x = [3, 4, 5, 6]
y = [86628.32966666667, 99093.81533333333, 103626.15866666667, 105788.62733333332]

x = np.array(x)
y = np.array(y)
I = 136657.9932077545

x = np.array(x)
y = np.array(y)
nCalc = 5

res_lsq = optimize.least_squares(func, x0=[1, 1], args=(x, y))
#res_lsq2 = optimize.least_squares(func2, x0=[1, 1, 1, 1], args=(x, y))

print(res_lsq.x)
xVals = range(2, 8)

A, B = res_lsq.x

yVals = [func([A, B], n, 0) for n in xVals]
#yVals2 = [func2(res_lsq2.x, n, 0) for n in xVals]


plt.figure()
plt.plot(x, y, 'rx')
plt.plot(xVals, yVals, 'b--')
#plt.plot(xVals, yVals2, 'g--')
plt.show()