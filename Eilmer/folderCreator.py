import os
import shutil

velocityList = range(2000, 11200, 200)
velocityList = [2000]
rhoRList = ['2.5e-5', '5.0e-5', '1.0e-4', '2.0e-4', '4.0e-4', '1.7e-3']

cwd = os.getcwd()
childDir = cwd + '\\CEA\\'
masterDir = cwd + '\\master\\'

for rhoR in rhoRList:
    for velocity in velocityList:
        folderName = 'pR' + rhoR + '-' + 'u' + str(velocity/1000)
        folderDir = childDir + folderName
        #os.mkdir(folderDir)
        shutil.copytree(masterDir, folderDir)

        with open(folderDir + '\\infProp.lua', 'w') as propertyFile:
            propertyFile.write('u_inf = ' + str(velocity) + ' -- m/s\n')
            propertyFile.write('rR = ' + rhoR + ' -- kg/m^2')