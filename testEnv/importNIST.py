import sys

def importFromFile(species):
    NistFilename = species + ' NIST.dat'

    badSymbols = [' ', '|', '[', ']', '-', 'i', 't', '+', 'x', '?', '*']

    with open(NistFilename) as NistFile:
        with open(species + '.dat', 'w') as writeFile:
            writeFile.write('J, epsion\n')

            next(NistFile)
            row = NistFile.readline()
            jIdx = row.index('J')

            next(NistFile)
            next(NistFile)

            for rowIdx, row in enumerate(NistFile):
                

                    
                j = row[jIdx-4:jIdx+5]
                level = row[38:60]

                for symbol in badSymbols:
                    j = j.replace(symbol, '')
                    level = level.replace(symbol, '')

                if '/' in j:
                    j = float(j.split('/')[0])/float(j.split('/')[1])
                    j = str(j)

                if j.strip() and level != '':
                    writeFile.write(j + ', ')
                    writeFile.write(level + '\n')
                else:
                    writeFile.write('\n')

importFromFile(sys.argv[1])