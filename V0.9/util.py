class Constants:
    def __init__(self):
        self.kB = 1.38064852E-23    
        self.electronCharge = 1.60217662E-19                    #
        self.R = 8.3144598                              #
        self.Cm_1ToJoules = self.electronCharge * 1.23986E-4    #
        #self.IH = 13.59844 * self.Cm_1ToJoules          #CRC Handbook of Chemistry and Physics (2003)
        self.eVToCm_1 = 8065.544005                     #
        self.IH = 13.59844 * self.eVToCm_1              #
        self.Rydberg = 13.605693009 * self.eVToCm_1
        self.vacPermiativity = 8.854187817E-12

subShellMap = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', \
               'l', 'm', 'n', 'o', 'q', 'r', 't', 'u', 'v', \
               'w', 'x', 'y', 'z', 'ss', 'pp', 'dd', 'ff', \
               'gg', 'hh',' ii', 'jj', 'kk', 'll', 'mm', 'nn',\
               'oo', 'qq', 'rr', 'tt', 'uu', 'vv', 'ww', 'xx',\
               'yy', 'zz', 'sss', 'ppp', 'ddd', 'fff', \
               'ggg', 'hhh',' iii', 'jjj', 'kkk', 'lll', 'mmm', 'nnn',\
               'ooo', 'qqq', 'rrr', 'ttt', 'uuu', 'vvv', 'www', 'xxx',\
               'yyy', 'zzz']
               
termMap = ['S', 'P', 'D', 'F', 'G', 'H', 'I', 'J', 'K', \
           'L', 'M', 'N', 'O', 'Q', 'R', 'T', 'U', 'V', \
           'W', 'X', 'Y', 'Z' 'SS', 'PP', 'DD', 'FF', 'GG',\
           'HH', 'II', 'JJ', 'KK', 'LL', 'MM', 'NN', 'OO',\
           'QQ', 'RR', 'TT', 'UU', 'VV', 'WW', 'XX', 'YY',\
           'ZZ', 'SSS', 'PPP', 'DDD', 'FFF', 'GGG',\
           'HHH', 'III', 'JJJ', 'KKK', 'LLL', 'MMM', 'NNN', 'OOO',\
           'QQQ', 'RRR', 'TTT', 'UUU', 'VVV', 'WWW', 'XXX', 'YYY',\
           'ZZZ']

def getSpectralDataDir():
    return 'spectraSource/'