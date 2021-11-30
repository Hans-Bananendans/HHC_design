#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 12:45:21 2021

@author: johan
"""

import numpy as np

def wirestack(i,j):
#              i=6 
#          <--------->
#       ^  o o o o o o
#       |   o o o o o 
#   j=5 |  o o o o o o
#       |   o o o o o
#       v  o o o o o o 
    return i*j-np.floor(j/2)



imax = 18
jmax = 16
results = np.zeros((imax,jmax))

for i in range(0,imax):
    for j in range(0,jmax):
        results[i][j] = wirestack(i,j)
        

# Solutions for 83 at (6,15), (8,11) and (17,5)
# Given an estimated coil width of ~22 mm and the diameter of AWG13 is 1.83 mm
#   it is most likely that the Helmholtz coil has a (8,11) setup:
#        
#                i=8
#         _______________________
#        |_______________________|
#        | |o o o o o o o o
#        | | o o o o o o o 
#        | |o o o o o o o o
#        | | o o o o o o o 
#        | |o o o o o o o o
#  j=11  | | o o o o o o o 
#        | |o o o o o o o o
#        | | o o o o o o o 
#        | |o o o o o o o o
#        | | o o o o o o o 
#        |_|o o o o o o o o______
#        |_______________________|