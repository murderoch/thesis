import matplotlib.pyplot as plt

cells = [30, 60, 120]

CEA = {}
with open('CEA.dat', 'r') as CEAFile:
    for row in CEAFile:
        data = row.split(' ')

        ID = (data[0], int(data[2]))
        #print(data[8])
        CEA[ID] = float(data[8])

Cap = {}
with open('Cap.dat', 'r') as CapFile:
    for row in CapFile:
        data = row.split(' ')

        ID = (data[0], int(data[2]))
        #print(data[8])
        Cap[ID] = float(data[8])

plt.figure(1)
plt.title('pR2.5e-5')
plt.figure(2)
plt.title('pR1.0e-4')

for keys, results in CEA.items():
    uIdx = keys[0].index('u')
    vel = float(keys[0][uIdx+1:])
    if keys[0][:8] == 'pR2.5e-5':
        plt.figure(1)
    else:
        plt.figure(2)

    if keys[1] == 120:
        plt.plot(vel, results, 'k+')
    elif keys[1] == 60:
        plt.plot(vel, results, 'g+')
    elif keys[1] == 30:
        plt.plot(vel, results, 'r+')

for keys, results in Cap.items():
    uIdx = keys[0].index('u')
    vel = float(keys[0][uIdx+1:])
    if keys[0][:8] == 'pR2.5e-5':
        plt.figure(1)
    else:
        plt.figure(2)

    if keys[1] == 120:
        plt.plot(vel, results, 'ko')
    elif keys[1] == 60:
        plt.plot(vel, results, 'go')
    elif keys[1] == 30:
        plt.plot(vel, results, 'ro')

plt.show()