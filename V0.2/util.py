def getSpectralDataDir():
    return 'spectraSource/'

class Constants:
    def __init__(self):
        self.kB = 1.38064852E-23                        #
        self.R = 8.3144598                              #
        self.Cm_1ToJoules = 1.60217E-19 * 1.23986E-4    #
        self.IH = 13.59844 * self.Cm_1ToJoules          #CRC Handbook of Chemistry and Physics (2003)
