class Constants:
    def __init__(self):
        self.kB = 1.38064852E-23                        #
        self.R = 8.3144598                              #
        self.Cm_1ToJoules = 1.60217E-19 * 1.23986E-4    #
        #self.IH = 13.59844 * self.Cm_1ToJoules          #CRC Handbook of Chemistry and Physics (2003)
        self.eVToCm_1 = 8065.544005                     #
        self.IH = 13.59844 * self.eVToCm_1              #


subShellMap = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', \
               'l', 'm', 'n', 'o', 'q', 'r', 't', 'u', 'v', \
               'w', 'x', 'y', 'z']
               
termMap = ['S', 'P', 'D', 'F', 'G', 'H', 'I', 'J', 'K', \
           'L', 'M', 'N', 'O', 'Q', 'R', 'T', 'U', 'V', \
           'W', 'X', 'Y', 'Z']

def getSpectralDataDir():
    return 'spectraSource/'