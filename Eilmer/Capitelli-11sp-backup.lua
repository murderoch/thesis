-- Auto-generated by prep-gas on: 04-Oct-2018 20:23:14

model = 'ThermallyPerfectGas'
species = {'N2', 'N2+', 'O2', 'O2+', 'N', 'N+', 'O', 'O+', 'NO', 'NO+', 'e-', }

db = {}
db['N2'] = {}
db['N2'].M = 0.02801340
db['N2'].sigma = 3.62100000
db['N2'].epsilon = 97.53000000
db['N2'].thermoCoeffs = {
  origin = 'Capitelli',
  nsegments = 4, 
  segment0 = {
    T_lower = 100.0,
    T_upper = 600.0,
    coeffs = {
      -1.129106765e+04,
       3.387587534e+02,
       2.501903115e+01,
       2.524957008e-02,
      -8.264339696e-05,
       1.310252018e-07,
      -7.061592582e-11,
       7.108460860e+02,
      -1.076003744e+01,
    }
  },
  segment1 = {
    T_lower = 600.0,
    T_upper = 10000.0,
    coeffs = {
       2.636648705e+05,
      -3.907918154e+03,
       3.461048957e+01,
       2.384106439e-03,
      -4.760812086e-07,
       3.040729356e-11,
       4.238089730e-16, 
       1.283210415e+04,
      -1.586640027e+01,
    }
  },
  segment2 = {
    T_lower = 10000.0,
    T_upper = 20000.0,
    coeffs = {
       1.000000000e+00,
       2.799520946e+05,
       3.726947237e+02,
      -1.174683142e-01,
       1.325013565e-05,
      -6.021237717e-10,
       9.630241714e-15, 
       4.938707040e+06,
      -1.672099740e+03,
    }
  },
  segment3 = {
    T_lower = 20000.0,
    T_upper = 50000.0,
    coeffs = {
       1.000000000e+00,
      -4.986436792e+06,
       1.014857931e+03,
      -6.265312980e-02,
       1.829018680e-06,
      -2.565469699e-11,
       1.402009591e-16, 
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
  nsegments = 4, 
  segment0 = {
    T_lower = 100.0,
    T_upper = 600.0,
    coeffs = {
      -2.085032553e+04,
       6.175387972e+02,
       2.191834689e+01,
       4.216661389e-02,
      -1.318557560e-04,
       2.043714058e-07,
      -1.117730024e-10, 
       1.790004424e+05,
       6.832974166e+00,
    }
  },
  segment1 = {
    T_lower = 600.0,
    T_upper = 5000.0,
    coeffs = {
       6.937309068e+06,
      -3.314867804e+04,
       8.428701501e+01,
      -3.885612695e-02,
       1.670968854e-05,
      -3.026589783e-09,
       1.975317090e-13, 
       1.340388483e+05,
       5.090897022e+01,
    }
  },
  segment2 = {
    T_lower = 5000.0,
    T_upper = 12000.0,
    coeffs = {
      -4.520071673e+00,
      -2.228295563e+05,
       1.826586087e+02,
      -3.277685441e-02,
       3.876754437e-06,
      -2.202140022e-10,
       4.687005180e-15, 
      -2.217361867e+06,
       8.436270947e+02,
    }
  },
  segment3 = {
    T_lower = 12000.0,
    T_upper = 50000.0,
    coeffs = {
       1.000000000e+00,
      -6.510671094e+05,
       2.046882470e+02,
      -1.271434605e-02,
       3.993703215e-07,
      -6.005904694e-12,
       3.508932960e-17, 
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
  nsegments = 4, 
  segment0 = {
    T_lower = 100.0,
    T_upper = 600.0,
    coeffs = {
       5.744384178e+04,
      -1.637741843e+03,
       4.690296471e+01,
      -9.237551182e-02,
       2.310710672e-04,
      -2.426680520e-07,
       9.577979365e-11, 
      -3.391454870e+03,
       1.849699470e+01,
    }
  },
  segment1 = {
    T_lower = 600.0,
    T_upper = 4000.0,
    coeffs = {
       3.046735999e+06,
      -1.709054929e+04,
       6.007365685e+01,
      -1.715866110e-02,
       7.240091563e-06,
      -1.322677088e-09,
       9.051981576e-14,  
      -1.689010929e+04,
       1.738716506e+01,
    }
  },
  segment2 = {
    T_lower = 4000.0,
    T_upper = 10000.0,
    coeffs = {
       1.000000000e+00,
      -7.302484847e+04, 
       1.034764471e+02,
      -2.335795812e-02,
       4.500216979e-06,
      -3.942032374e-10,
       1.233272713e-14, 
       2.293554027e+06,
      -5.530621610e+02,
    }
  },
  segment3 = {
    T_lower = 10000.0,
    T_upper = 50000.0,
    coeffs = {
       1.000000000e+00,
      -1.219251804e+03,
       7.572344370e+01,
      -4.937086952e-03,
       1.869972023e-07,
      -3.253756015e-12,
       2.134575541e-17, 
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
  nsegments = 4, 
  segment0 = {
    T_lower = 100,
    T_upper = 600.0,
    coeffs = {
       5.017108592e+02,
       1.733163473e+01,
       2.830578603e+01,
       1.013734967e-02,
      -5.525799734e-05,
       1.281584333e-07,
      -8.586280594e-11, 
       1.345544668e+05,
       2.902709750e+01,
    }
  },
  segment1 = {
    T_lower = 600.0,
    T_upper = 4000.0,
    coeffs = {
       2.246339449e+06,
      -1.072752442e+04,
       4.278342779e+01,
       4.472812495e-04,
      -1.130911106e-06,
       3.679949819e-10,
      -3.600847363e-14, 
       1.446321044e+05,
      -5.811230650e+00,
    }
  },
  segment2 = {
    T_lower = 4000.0,
    T_upper = 14000.0,
    coeffs = {
       1.000000000e+00,
      -1.982234170e+05,
       2.019384523e+02,
      -5.174422583e-02,
       7.560562328e-06,
      -4.834630433e-10,
       1.110114565e-14,
      -8.857866270e+06,
       2.852035602e+03,
    }
  },
  segment3 = {
    T_lower = 14000.0,
    T_upper = 50000.0,
    coeffs = {
       1.000000000e+00,
      -3.541917709e+05,
       1.704387219e+02,
      -1.135183659e-02,
       3.756837351e-07,
      -5.846473929e-12,
       3.497837650e-17,  
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
  nsegments = 4, 
  segment0 = {
    T_lower = 100.0,
    T_upper = 6000.0,
    coeffs = {
      -5.590000755e+03,
       9.659666997e+01,
       2.028729896e+01,
       9.572521531e-04,
      -7.102182217e-07,
       1.930899393e-10,
      -1.285558277e-14, 
       5.610463780e+04,
       4.193905036e+00,
    }
  },
  segment1 = {
    T_lower = 6000.0,
    T_upper = 14000.0,
    coeffs = {
       1.000000000e+00,
       1.786326718e+06,
      -9.178987451e+02,
       1.795820294e-01,
      -1.435530622e-05,
       3.680822777e-10,
       4.249124861e-15, 
       5.697351330e+04,
       4.865231506e+00,
    }
  },
  segment2 = {
    T_lower = 14000.0,
    T_upper = 22000.0,
    coeffs = {
       1.000000000e+00,
       7.158015501e+08,
      -1.917832823e+05,
       2.005332998e+01,
      -1.021609579e-03,
       2.541348504e-08,
      -2.474915338e-13, 
       2.550585618e+06,
      -5.848769753e+02,
    }
  },
  segment3 = {
    T_lower = 22000.0,
    T_upper = 50000.0,
    coeffs = {
       1.000000000e+00,
       6.754335060e+07,
      -8.347472926e+03,
       4.178817705e-01,
      -1.045944997e-05,
       1.308061834e-10,
      -6.525470987e-16, 
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
    T_lower = 100.0,
    T_upper = 14000.0,
    coeffs = {
       3.865724693e+04,
      -1.482683058e+01,
       2.109662652e+01,
      -5.758445296e-04,
       2.510658467e-07,
      -2.471227560e-11,
       7.646392900e-16, 
       2.256284738e+05,
       5.076830786e+00,
    }
  },
  segment1 = {
    T_lower = 14000.0,
    T_upper = 30000.0,
    coeffs = {
       1.000000000e+00,
       1.065739478e+08,
      -2.759773013e+04,
       2.809070599e+00,
      -1.396969008e-04,
       3.380476647e-09,
      -3.156583987e-14, 
       2.310809984e+05,
      -1.994146545e+00,
    }
  },
  segment2 = {
    T_lower = 30000.0,
    T_upper = 50000.0,
    coeffs = {
       1.000000000e+00,
      -1.472353695e+09,
       1.623665113e+05,
      -6.943885753e+00,
       1.446403525e-04,
      -1.468153781e-09,
       5.796401890e-15, 
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
  nsegments = 4, 
  segment0 = {
    T_lower = 100.0,
    T_upper = 6000.0,
    coeffs = {
      -5.245989807e+04,
       1.022292838e+03,
       1.851927477e+01,
       2.265275401e-03,
      -1.052097686e-06,
       2.231996977e-10,
      -1.570453037e-14, 
       2.840362437e+04,
       8.404241820e+00,
    }
  },
  segment1 = {
    T_lower = 6000.0,
    T_upper = 14000.0,
    coeffs = {
       8.077273429e-01,
       2.667359043e+05,
      -8.262243182e+01,
       7.021227396e-03,
       1.894997890e-06,
      -3.077756113e-10,
       1.278844755e-14, 
       3.392428060e+04,
      -6.679585350e-01,
    }
  },
  segment2 = {
    T_lower = 14000.0,
    T_upper = 22000.0,
    coeffs = {
       1.000000000e+00,
      -1.761582783e+08,
       6.102013693e+04,
      -8.190825779e+00,
       5.325267720e-04,
      -1.673312132e-08,
       2.037565818e-13, 
       8.890942630e+05,
      -2.181728151e+02,
    }
  },
  segment3 = {
    T_lower = 22000.0,
    T_upper = 50000.0,
    coeffs = {
       1.000000000e+00,
       1.967424433e+07,
      -1.111494939e+03,
      -1.883089050e-03,
       1.367393682e-06,
      -3.204795298e-11,
       2.274579777e-16, 
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
    T_lower = 100.0,
    T_upper = 14000.0,
    coeffs = {
       8.101330223e+02,
       1.372643755e-01,
       2.068773813e+01,
       3.500623539e-04,
      -2.629218475e-07,
       5.173605857e-11,
      -2.236857505e-15, 
       1.879352842e+05,
       4.393376760e+00,
    }
  },
  segment1 = {
    T_lower = 14000.0,
    T_upper = 32000.0,
    coeffs = {
       1.000000000e+00,
       1.973225528e+06,
       -3.354957071e+02,
       1.767521682e-02,
       4.313176661e-07,
       -5.413291259e-11,
       1.083334764e-15, 
       1.837191966e+05,
       1.005690382e+01,
    }
  },
  segment2 = {
    T_lower = 32000.0,
    T_upper = 50000.0,
    coeffs = {
       1.000000000e+00,
       1.000000000e+00,
       1.480726772e+04,
      -1.673498335e+00,
       6.875747693e-05,
      -1.212151580e-09,
       7.773305141e-15, 
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
  nsegments = 4, 
  segment0 = {
    T_lower = 100.0,
    T_upper = 600.0,
    coeffs = {
       1.044566445e+04,
      -2.763600573e+02,
       3.170615882e+01,
      -1.000730380e-02,
       9.788202706e-06,
       1.574749417e-08,
      -2.472300092e-17,  
       9.098214410e+03,
       6.728725490e+00,
    }
  },
  segment1 = {
    T_lower = 600.0,
    T_upper = 4000.0,
    coeffs = {
       2.373252190e+06,
      -1.129079898e+04,
       4.376721986e+01,
      -5.578895730e-04,
      -5.565517909e-07,
       1.547229603e-10,
      -4.829677995e-22, 
       1.750317656e+04,
      -8.501669090e+00,
    }
  },
  segment2 = {
    T_lower = 4000.0,
    T_upper = 15000.0,
    coeffs = {
       1.000000000e+00,
      -9.146969123e+04,
       1.042888146e+02,
      -1.797499614e-02,
       2.110097746e-06,
      -8.455209317e-11,
       2.295046110e-24,  
      -4.677501240e+06,
       1.242081216e+03,
    }
  },
  segment3 = {
    T_lower = 15000.0,
    T_upper = 50000.0,
    coeffs = {
       1.000000000e+00,
      -8.181616280e+05,
       2.178428749e+02,
      -1.173147234e-02,
       2.954603736e-07,
      -2.859699048e-12,
       1.743856119e-27,   
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
  nsegments = 4, 
  segment0 = {
    T_lower = 100.0,
    T_upper = 600.0,
    coeffs = {
      -1.917720143e+04,
       5.345975362e+02,
       2.334518906e+01,
       3.034515510e-02,
      -7.920663723e-05,
       8.414912057e-08,
      -4.197872027e-17,  
       1.187495132e+05,
      -4.398433810e+00,
    }
  },
  segment1 = {
    T_lower = 600.0,
    T_upper = 12000.0,
    coeffs = {
       4.790599864e+06,
      -1.778128874e+04,
       4.742748091e+01,
      -1.650854289e-03,
       8.186429456e-11,
       1.596397161e-11,
       1.338670909e-25,
       1.443308939e+07,
      -4.324044462e+03,
    }
  },
  segment2 = {
    T_lower = 12000.0,
    T_upper = 20000.0,
    coeffs = {
       1.000000000e+00,
       9.947184157e+06,
      -2.886816035e+03,
       3.143926162e-01,
      -1.458841634e-05,
       2.560817443e-10,
      -8.517006662e-25, 
       1.443308939e+07,
      -4.324044462e+03,
    }
  },
  segment3 = {
    T_lower = 20000.0,
    T_upper = 50000.0,
    coeffs = {
       1.000000000e+00,
      -3.560023444e+06,
       7.641397536e+02,
      -4.516404876e-02,
       1.165050584e-06,
      -1.151626986e-11, 
       7.467918393e-27, 
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
    T_upper = 50000.0,
    coeffs = {
       0.000000000e+00,
       0.000000000e+00,
       2.079000000e+00,
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