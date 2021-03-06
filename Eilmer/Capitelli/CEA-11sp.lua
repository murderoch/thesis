-- Auto-generated by prep-gas on: 04-Oct-2018 20:23:14

model = 'ThermallyPerfectGas'
species = {'N2', 'N2+', 'O2', 'O2+', 'N', 'N+', 'O', 'O+', 'NO', 'NO+', 'e-', }

db = {}
db['N2'] = {}
db['N2'].M = 0.02801340
db['N2'].sigma = 3.62100000
db['N2'].epsilon = 97.53000000
db['N2'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  segment0 = {
    T_lower = 200.0,
    T_upper = 1000.0,
    coeffs = {
       2.210371497e+04,
      -3.818461820e+02,
       6.082738360e+00,
      -8.530914410e-03,
       1.384646189e-05,
      -9.625793620e-09,
       2.519705809e-12,
       7.108460860e+02,
      -1.076003744e+01,
    }
  },
  segment1 = {
    T_lower = 1000.0,
    T_upper = 6000.0,
    coeffs = {
       5.877124060e+05,
      -2.239249073e+03,
       6.066949220e+00,
      -6.139685500e-04,
       1.491806679e-07,
      -1.923105485e-11,
       1.061954386e-15,
       1.283210415e+04,
      -1.586640027e+01,
    }
  },
  segment2 = {
    T_lower = 6000.0,
    T_upper = 20000.0,
    coeffs = {
       8.310139160e+08,
      -6.420733540e+05,
       2.020264635e+02,
      -3.065092046e-02,
       2.486903333e-06,
      -9.705954110e-11,
       1.437538881e-15,
       4.938707040e+06,
      -1.672099740e+03,
    }
  },
}
db['N2'].viscosity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.2526577e-01,
      B = -3.1779652e+01,
      C = -1.6407983e+03,
      D =  1.7454992e+00,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.7395209e-01,
      B =  5.6152222e+02,
      C = -1.7394809e+05,
      D = -3.9335958e-01,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.8503551e-01,
      B =  9.0902171e+02,
      C = -7.3129061e+05,
      D = -5.3503838e-01,
   },
}
db['N2'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  8.5439436e-01,
      B =  1.0573224e+02,
      C = -1.2347848e+04,
      D =  4.7793128e-01,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.8407146e-01,
      B =  1.3357293e+02,
      C = -1.1429640e+04,
      D =  2.4417019e-01,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  2.4176185e+00,
      B =  8.0477749e+03,
      C =  3.1055802e+06,
      D = -1.4517761e+01,
   },
}
db['N2+'] = {}
db['N2+'].M = 0.02801285
db['N2+'].sigma = 3.62100000
db['N2+'].epsilon = 97.53000000
db['N2+'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  segment0 = {
    T_lower = 298.1,
    T_upper = 1000.0,
    coeffs = {
      -3.474047470e+04,
       2.696222703e+02,
       3.164916370e+00,
      -2.132239781e-03,
       6.730476400e-06,
      -5.637304970e-09,
       1.621756000e-12,
       1.790004424e+05,
       6.832974166e+00,
    }
  },
  segment1 = {
    T_lower = 1000.0,
    T_upper = 6000.0,
    coeffs = {
      -2.845599002e+06,
       7.058893030e+03,
      -2.884886385e+00,
       3.068677059e-03,
      -4.361652310e-07,
       2.102514545e-11,
       5.411996470e-16,
       1.340388483e+05,
       5.090897022e+01,
    }
  },
  segment2 = {
    T_lower = 6000.0,
    T_upper = 20000.0,
    coeffs = {
      -3.712829770e+08,
       3.139287234e+05,
      -9.603518050e+01,
       1.571193286e-02,
      -1.175065525e-06,
       4.144441230e-11,
      -5.621893090e-16,
      -2.217361867e+06,
       8.436270947e+02,
    }
  },
}
db['N2+'].viscosity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.2526577e-01,
      B = -3.1779652e+01,
      C = -1.6407983e+03,
      D =  1.7454992e+00,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.7395209e-01,
      B =  5.6152222e+02,
      C = -1.7394809e+05,
      D = -3.9335958e-01,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.8503551e-01,
      B =  9.0902171e+02,
      C = -7.3129061e+05,
      D = -5.3503838e-01,
   },
}
db['N2+'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  8.5439436e-01,
      B =  1.0573224e+02,
      C = -1.2347848e+04,
      D =  4.7793128e-01,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.8407146e-01,
      B =  1.3357293e+02,
      C = -1.1429640e+04,
      D =  2.4417019e-01,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  2.4176185e+00,
      B =  8.0477749e+03,
      C =  3.1055802e+06,
      D = -1.4517761e+01,
   },
}
db['O2'] = {}
db['O2'].M = 0.03199880
db['O2'].sigma = 3.45800000
db['O2'].epsilon = 107.40000000
db['O2'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  segment0 = {
    T_lower = 200.0,
    T_upper = 1000.0,
    coeffs = {
      -3.425563420e+04,
       4.847000970e+02,
       1.119010961e+00,
       4.293889240e-03,
      -6.836300520e-07,
      -2.023372700e-09,
       1.039040018e-12,
      -3.391454870e+03,
       1.849699470e+01,
    }
  },
  segment1 = {
    T_lower = 1000.0,
    T_upper = 6000.0,
    coeffs = {
      -1.037939022e+06,
       2.344830282e+03,
       1.819732036e+00,
       1.267847582e-03,
      -2.188067988e-07,
       2.053719572e-11,
      -8.193467050e-16,
      -1.689010929e+04,
       1.738716506e+01,
    }
  },
  segment2 = {
    T_lower = 6000.0,
    T_upper = 20000.0,
    coeffs = {
       4.975294300e+08,
      -2.866106874e+05,
       6.690352250e+01,
      -6.169959020e-03,
       3.016396027e-07,
      -7.421416600e-12,
       7.278175770e-17,
       2.293554027e+06,
      -5.530621610e+02,
    }
  },
}
db['O2'].viscosity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.0916180e-01,
      B = -5.2244847e+01,
      C = -5.9974009e+02,
      D =  2.0410801e+00,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.2216486e-01,
      B =  1.7550839e+02,
      C = -5.7974816e+04,
      D =  1.0901044e+00,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  7.3981127e-01,
      B =  3.9194906e+02,
      C = -3.7833168e+05,
      D =  9.0931780e-01,
   },
}
db['O2'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  7.7229167e-01,
      B =  6.8463210e+00,
      C = -5.8933377e+03,
      D =  1.2210365e+00,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  9.0917351e-01,
      B =  2.9124182e+02,
      C = -7.9650171e+04,
      D =  6.4851631e-02,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -1.1218262e+00,
      B = -1.9286378e+04,
      C =  2.3295011e+07,
      D =  2.0342043e+01,
   },
}
db['O2+'] = {}
db['O2+'].M = 0.03199825
db['O2+'].sigma = 3.62100000
db['O2+'].epsilon = 97.53000000
db['O2+'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  segment0 = {
    T_lower = 298.1,
    T_upper = 1000.0,
    coeffs = {
      -8.607205450e+04,
       1.051875934e+03,
      -5.432380470e-01,
       6.571166540e-03,
      -3.274263750e-06,
       5.940645340e-11,
       3.238784790e-13,
       1.345544668e+05,
       2.902709750e+01,
    }
  },
  segment1 = {
    T_lower = 1000.0,
    T_upper = 6000.0,
    coeffs = {
       7.384654880e+04,
      -8.459559540e+02,
       4.985164160e+00,
      -1.611010890e-04,
       6.427083990e-08,
      -1.504939874e-11,
       1.578465409e-15,
       1.446321044e+05,
      -5.811230650e+00,
    }
  },
  segment2 = {
    T_lower = 6000.0,
    T_upper = 20000.0,
    coeffs = {
      -1.562125524e+09,
       1.161406778e+06,
      -3.302504720e+02,
       4.710937520e-02,
      -3.354461380e-06,
       1.167968599e-10,
      -1.589754791e-15,
      -8.857866270e+06,
       2.852035602e+03,
    }
  },
}
db['O2+'].viscosity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.0916180e-01,
      B = -5.2244847e+01,
      C = -5.9974009e+02,
      D =  2.0410801e+00,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.2216486e-01,
      B =  1.7550839e+02,
      C = -5.7974816e+04,
      D =  1.0901044e+00,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  7.3981127e-01,
      B =  3.9194906e+02,
      C = -3.7833168e+05,
      D =  9.0931780e-01,
   },
}
db['O2+'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  7.7229167e-01,
      B =  6.8463210e+00,
      C = -5.8933377e+03,
      D =  1.2210365e+00,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  9.0917351e-01,
      B =  2.9124182e+02,
      C = -7.9650171e+04,
      D =  6.4851631e-02,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -1.1218262e+00,
      B = -1.9286378e+04,
      C =  2.3295011e+07,
      D =  2.0342043e+01,
   },
}
db['N'] = {}
db['N'].M = 0.01400670
db['N'].sigma = 3.29800000
db['N'].epsilon = 71.40000000
db['N'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  segment0 = {
    T_lower = 200.0,
    T_upper = 1000.0,
    coeffs = {
       0.000000000e+00,
       0.000000000e+00,
       2.500000000e+00,
       0.000000000e+00,
       0.000000000e+00,
       0.000000000e+00,
       0.000000000e+00,
       5.610463780e+04,
       4.193905036e+00,
    }
  },
  segment1 = {
    T_lower = 1000.0,
    T_upper = 6000.0,
    coeffs = {
       8.876501380e+04,
      -1.071231500e+02,
       2.362188287e+00,
       2.916720081e-04,
      -1.729515100e-07,
       4.012657880e-11,
      -2.677227571e-15,
       5.697351330e+04,
       4.865231506e+00,
    }
  },
  segment2 = {
    T_lower = 6000.0,
    T_upper = 20000.0,
    coeffs = {
       5.475181050e+08,
      -3.107574980e+05,
       6.916782740e+01,
      -6.847988130e-03,
       3.827572400e-07,
      -1.098367709e-11,
       1.277986024e-16,
       2.550585618e+06,
      -5.848769753e+02,
    }
  },
}
db['N'].viscosity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.3724737e-01,
      B =  4.3997150e+02,
      C = -1.7450753e+05,
      D =  1.0365689e-01,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.9986588e-01,
      B =  1.4112801e+03,
      C = -1.8200478e+06,
      D = -5.5811716e-01,
   },
}
db['N'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.3771661e-01,
      B =  4.4243270e+02,
      C = -1.7578446e+05,
      D =  8.9942915e-01,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  9.0001710e-01,
      B =  1.4141175e+03,
      C = -1.8262403e+06,
      D =  2.4048513e-01,
   },
}
db['N+'] = {}
db['N+'].M = 0.01400615
db['N+'].sigma = 3.62100000
db['N+'].epsilon = 97.53000000
db['N+'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  segment0 = {
    T_lower = 298.1,
    T_upper = 1000.0,
    coeffs = {
       5.237079210e+03,
       2.299958315e+00,
       2.487488821e+00,
       2.737490756e-05,
      -3.134447576e-08,
       1.850111332e-11,
      -4.447350984e-15,
       2.256284738e+05,
       5.076830786e+00,
    }
  },
  segment1 = {
    T_lower = 1000.0,
    T_upper = 6000.0,
    coeffs = {
       2.904970374e+05,
      -8.557908610e+02,
       3.477389290e+00,
      -5.288267190e-04,
       1.352350307e-07,
      -1.389834122e-11,
       5.046166279e-16,
       2.310809984e+05,
      -1.994146545e+00,
    }
  },
  segment2 = {
    T_lower = 6000.0,
    T_upper = 20000.0,
    coeffs = {
       1.646092148e+07,
      -1.113165218e+04,
       4.976986640e+00,
      -2.005393583e-04,
       1.022481356e-08,
      -2.691430863e-13,
       3.539931593e-18,
       3.136284696e+05,
      -1.706646380e+01,
    }
  },
}
db['N+'].viscosity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.3724737e-01,
      B =  4.3997150e+02,
      C = -1.7450753e+05,
      D =  1.0365689e-01,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.9986588e-01,
      B =  1.4112801e+03,
      C = -1.8200478e+06,
      D = -5.5811716e-01,
   },
}
db['N+'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.3771661e-01,
      B =  4.4243270e+02,
      C = -1.7578446e+05,
      D =  8.9942915e-01,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  9.0001710e-01,
      B =  1.4141175e+03,
      C = -1.8262403e+06,
      D =  2.4048513e-01,
   },
}
db['O'] = {}
db['O'].M = 0.01599940
db['O'].sigma = 2.75000000
db['O'].epsilon = 80.00000000
db['O'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  segment0 = {
    T_lower = 200.0,
    T_upper = 1000.0,
    coeffs = {
      -7.953611300e+03,
       1.607177787e+02,
       1.966226438e+00,
       1.013670310e-03,
      -1.110415423e-06,
       6.517507500e-10,
      -1.584779251e-13,
       2.840362437e+04,
       8.404241820e+00,
    }
  },
  segment1 = {
    T_lower = 1000.0,
    T_upper = 6000.0,
    coeffs = {
       2.619020262e+05,
      -7.298722030e+02,
       3.317177270e+00,
      -4.281334360e-04,
       1.036104594e-07,
      -9.438304330e-12,
       2.725038297e-16,
       3.392428060e+04,
      -6.679585350e-01,
    }
  },
  segment2 = {
    T_lower = 6000.0,
    T_upper = 20000.0,
    coeffs = {
       1.779004264e+08,
      -1.082328257e+05,
       2.810778365e+01,
      -2.975232262e-03,
       1.854997534e-07,
      -5.796231540e-12,
       7.191720164e-17,
       8.890942630e+05,
      -2.181728151e+02,
    }
  },
}
db['O'].viscosity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.7269241e-01,
      B =  8.3842977e+01,
      C = -5.8502098e+04,
      D =  8.5100827e-01,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.7669586e-01,
      B =  1.0158420e+03,
      C = -1.0884566e+06,
      D = -1.8001077e-01,
   },
}
db['O'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.7271664e-01,
      B =  8.3989100e+01,
      C = -5.8580966e+04,
      D =  1.5179900e+00,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.7676666e-01,
      B =  1.0170744e+03,
      C = -1.0906690e+06,
      D =  4.8644232e-01,
   },
}
db['O+'] = {}
db['O+'].M = 0.01599885
db['O+'].sigma = 3.62100000
db['O+'].epsilon = 97.53000000
db['O+'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  segment0 = {
    T_lower = 298.1,
    T_upper = 1000.0,
    coeffs = {
       0.000000000e+00,
       0.000000000e+00,
       2.500000000e+00,
       0.000000000e+00,
       0.000000000e+00,
       0.000000000e+00,
       0.000000000e+00,
       1.879352842e+05,
       4.393376760e+00,
    }
  },
  segment1 = {
    T_lower = 1000.0,
    T_upper = 6000.0,
    coeffs = {
      -2.166513208e+05,
       6.665456150e+02,
       1.702064364e+00,
       4.714992810e-04,
      -1.427131823e-07,
       2.016595903e-11,
      -9.107157762e-16,
       1.837191966e+05,
       1.005690382e+01,
    }
  },
  segment2 = {
    T_lower = 6000.0,
    T_upper = 20000.0,
    coeffs = {
      -2.143835383e+08,
       1.469518523e+05,
      -3.680864540e+01,
       5.036164540e-03,
      -3.087873854e-07,
       9.186834870e-12,
      -1.074163268e-16,
      -9.614208960e+05,
       3.426193080e+02,
    }
  },
}
db['O+'].viscosity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.7269241e-01,
      B =  8.3842977e+01,
      C = -5.8502098e+04,
      D =  8.5100827e-01,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.7669586e-01,
      B =  1.0158420e+03,
      C = -1.0884566e+06,
      D = -1.8001077e-01,
   },
}
db['O+'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.7271664e-01,
      B =  8.3989100e+01,
      C = -5.8580966e+04,
      D =  1.5179900e+00,
   },
   segment1 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.7676666e-01,
      B =  1.0170744e+03,
      C = -1.0906690e+06,
      D =  4.8644232e-01,
   },
}
db['NO'] = {}
db['NO'].M = 0.03000610
db['NO'].sigma = 3.62100000
db['NO'].epsilon = 97.53000000
db['NO'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  segment0 = {
    T_lower = 200.0,
    T_upper = 1000.0,
    coeffs = {
      -1.143916503e+04,
       1.536467592e+02,
       3.431468730e+00,
      -2.668592368e-03,
       8.481399120e-06,
      -7.685111050e-09,
       2.386797655e-12,
       9.098214410e+03,
       6.728725490e+00,
    }
  },
  segment1 = {
    T_lower = 1000.0,
    T_upper = 6000.0,
    coeffs = {
       2.239018716e+05,
      -1.289651623e+03,
       5.433936030e+00,
      -3.656034900e-04,
       9.880966450e-08,
      -1.416076856e-11,
       9.380184620e-16,
       1.750317656e+04,
      -8.501669090e+00,
    }
  },
  segment2 = {
    T_lower = 6000.0,
    T_upper = 20000.0,
    coeffs = {
      -9.575303540e+08,
       5.912434480e+05,
      -1.384566826e+02,
       1.694339403e-02,
      -1.007351096e-06,
       2.912584076e-11,
      -3.295109350e-16,
      -4.677501240e+06,
       1.242081216e+03,
    }
  },
}
db['NO'].viscosity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.0262029e-01,
      B = -6.2017783e+01,
      C = -1.3954524e+02,
      D =  2.0268332e+00,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.8009050e-01,
      B =  3.0486891e+02,
      C = -9.4847722e+04,
      D =  5.2873381e-01,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.0580582e-01,
      B =  6.2427878e+02,
      C = -5.7879210e+05,
      D =  2.6516450e-01,
   },
}
db['NO'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  9.5028758e-01,
      B =  7.6667058e+01,
      C = -9.9894764e+03,
      D = -6.2776717e-03,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.6215238e-01,
      B =  4.4568223e+02,
      C = -2.3856466e+05,
      D =  4.6209876e-01,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -1.0377865e+00,
      B = -3.4486864e+04,
      C =  6.7451187e+07,
      D =  2.0928749e+01,
   },
}
db['NO+'] = {}
db['NO+'].M = 0.03000555
db['NO+'].sigma = 3.62100000
db['NO+'].epsilon = 97.53000000
db['NO+'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  segment0 = {
    T_lower = 298.1,
    T_upper = 1000.0,
    coeffs = {
       1.398106635e+03,
      -1.590446941e+02,
       5.122895400e+00,
      -6.394388620e-03,
       1.123918342e-05,
      -7.988581260e-09,
       2.107383677e-12,
       1.187495132e+05,
      -4.398433810e+00,
    }
  },
  segment1 = {
    T_lower = 1000.0,
    T_upper = 6000.0,
    coeffs = {
       6.069876900e+05,
      -2.278395427e+03,
       6.080324670e+00,
      -6.066847580e-04,
       1.432002611e-07,
      -1.747990522e-11,
       8.935014060e-16,
       1.322709615e+05,
      -1.519880037e+01,
    }
  },
  segment2 = {
    T_lower = 6000.0,
    T_upper = 20000.0,
    coeffs = {
       2.676400347e+09,
      -1.832948690e+06,
       5.099249390e+02,
      -7.113819280e-02,
       5.317659880e-06,
      -1.963208212e-10,
       2.805268230e-15,
       1.443308939e+07,
      -4.324044462e+03,
    }
  },
}
db['NO+'].viscosity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.0262029e-01,
      B = -6.2017783e+01,
      C = -1.3954524e+02,
      D =  2.0268332e+00,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.8009050e-01,
      B =  3.0486891e+02,
      C = -9.4847722e+04,
      D =  5.2873381e-01,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.0580582e-01,
      B =  6.2427878e+02,
      C = -5.7879210e+05,
      D =  2.6516450e-01,
   },
}
db['NO+'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  9.5028758e-01,
      B =  7.6667058e+01,
      C = -9.9894764e+03,
      D = -6.2776717e-03,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.6215238e-01,
      B =  4.4568223e+02,
      C = -2.3856466e+05,
      D =  4.6209876e-01,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -1.0377865e+00,
      B = -3.4486864e+04,
      C =  6.7451187e+07,
      D =  2.0928749e+01,
   },
}
db['e-'] = {}
db['e-'].M = 0.00000055
db['e-'].sigma = 3.62100000
db['e-'].epsilon = 97.53000000
db['e-'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 1, 
  segment0 = {
    T_lower = 298.1,
    T_upper = 20000.0,
    coeffs = {
       0.000000000e+00,
       0.000000000e+00,
       2.500000000e+00,
       0.000000000e+00,
       0.000000000e+00,
       0.000000000e+00,
       0.000000000e+00,
      -7.453750000e+02,
      -1.172081224e+01,
    }
  },
}
db['e-'].viscosity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.2526577e-01,
      B = -3.1779652e+01,
      C = -1.6407983e+03,
      D =  1.7454992e+00,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.7395209e-01,
      B =  5.6152222e+02,
      C = -1.7394809e+05,
      D = -3.9335958e-01,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.8503551e-01,
      B =  9.0902171e+02,
      C = -7.3129061e+05,
      D = -5.3503838e-01,
   },
}
db['e-'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  8.5439436e-01,
      B =  1.0573224e+02,
      C = -1.2347848e+04,
      D =  4.7793128e-01,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.8407146e-01,
      B =  1.3357293e+02,
      C = -1.1429640e+04,
      D =  2.4417019e-01,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  2.4176185e+00,
      B =  8.0477749e+03,
      C =  3.1055802e+06,
      D = -1.4517761e+01,
   },
}