import os
import shutil

velocityList = xrange(2000, 11500, 500)
#velocityList = [2000]

#rhoRList = ['2.5e-5', '5.0e-5', '1.0e-4', '2.0e-4', '4.0e-4', '1.7e-3']
rhoRList = ['2.5e-5', '1.0e-4']# '4.0e-4', '1.7e-3']

cwd = os.getcwd()
childDir = cwd + '/test/'
masterDir = cwd + '/master/'

for i in xrange(len(rhoRList)):
	rhoR = rhoRList[i]
	for j in xrange(len(velocityList)):
		velocity = velocityList[j]
		
		velocityKm = float(velocity)/1000.
		velocityStr = "%.1f" % velocityKm
		
		folderName = 'pR' + rhoR + '-' + 'u' + velocityStr

		print(folderName)

		folderDir = childDir + folderName

		shutil.copytree(masterDir, folderDir)
		
		#print(str(velocity))

		with open(folderDir + '/infProp.lua', 'w') as propertyFile:
			propertyFile.write('u_inf = ' + str(velocity) + ' -- m/s\n')
			propertyFile.write('rR = ' + rhoR + ' -- kg/m^2')



