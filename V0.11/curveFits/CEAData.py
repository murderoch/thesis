
def getCEA():
    
    '''
    'N+'
    CEA = {(200., 1000.): [
       5.237079210e+03,
       2.299958315e+00,
       2.487488821e+00,
       2.737490756e-05,
      -3.134447576e-08,
       1.850111332e-11,
      -4.447350984e-15,
       2.256284738e+05,
       5.076830786e+00,],

    (1000., 6000.): [
       2.904970374e+05,
      -8.557908610e+02,
       3.477389290e+00,
      -5.288267190e-04,
       1.352350307e-07,
      -1.389834122e-11,
       5.046166279e-16,
       2.310809984e+05,
      -1.994146545e+00,],

    (6000., 50000.): [
       1.646092148e+07,
      -1.113165218e+04,
       4.976986640e+00,
      -2.005393583e-04,
       1.022481356e-08,
      -2.691430863e-13,
       3.539931593e-18,
       3.136284696e+05,
      -1.706646380e+01,]}

    '''
    
    'O'
    CEA = {(200., 1000.): [
      -7.953611300e+03,
       1.607177787e+02,
       1.966226438e+00,
       1.013670310e-03,
      -1.110415423e-06,
       6.517507500e-10,
      -1.584779251e-13,
       2.840362437e+04,
       8.404241820e+00],

    (1000., 6000.): [
       2.619020262e+05,
      -7.298722030e+02,
       3.317177270e+00,
      -4.281334360e-04,
       1.036104594e-07,
      -9.438304330e-12,
       2.725038297e-16,
       3.392428060e+04,
      -6.679585350e-01],

    (6000., 50000.): [
       1.779004264e+08,
      -1.082328257e+05,
       2.810778365e+01,
      -2.975232262e-03,
       1.854997534e-07,
      -5.796231540e-12,
       7.191720164e-17,
       8.890942630e+05,
      -2.181728151e+02]}
    

    '''
    CEA = {(200, 1000): [2.210371497e+04,
                    -3.818461820e+02,
                    6.082738360e+00,
                    -8.530914410e-03,
                    1.384646189e-05,
                    -9.625793620e-09,
                    2.519705809e-12,
                    7.108460860e+02,
                    -1.076003744e+01],
        
        (1000, 6000): [5.877124060e+05,
                        -2.239249073e+03,
                        6.066949220e+00,
                        -6.139685500e-04,
                        1.491806679e-07,
                        -1.923105485e-11,
                        1.061954386e-15,
                        1.283210415e+04,
                        -1.586640027e+01,],
            (6000, 50000): [8.310139160e+08,
                        -6.420733540e+05,
                        2.020264635e+02,
                        -3.065092046e-02,
                        2.486903333e-06,
                        -9.705954110e-11,
                        1.437538881e-15,
                        4.938707040e+06,
                        -1.672099740e+03,]}
        '''

    return CEA