# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 15:49:21 2023

@author: Main
"""

import numpy as np

V = np.array([0.8, 1.0, 1.2, 1.5, 2.0])

A = np.array([
    [0.26,  0.32, 0.385, 0.48, 0.64],
    [0.255, 0.32, 0.38,  0.47, 0.63],
    [0.23,  0.29, 0.345, 0.42, 0.57],
    [0.24,  0.30, 0.36,  0.44, 0.60],
    [0.23,  0.29, 0.35,  0.43, 0.58],
    [0.23,  0.29, 0.35,  0.43, 0.58]
    ])

R = np.zeros((6,5))
Ravg = np.zeros((6,1))

for coil in range(6):
    for i in range(5):
        R[coil,i] = V[i]/A[coil,i]
    
for coil in range(6):
    Ravg[coil,0] = sum(R[coil,:])/len(R[coil,:])

print(R.round(1))
print(Ravg.round(2))