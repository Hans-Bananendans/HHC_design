#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 15:26:12 2021

@author: Johan Monster
"""

import time
import numpy as np
import matplotlib.pyplot as plt

from hhc_classes import PowerSupply, HHCoil, HHCage


time1 = time.perf_counter()

def scale_vector(vB: np.ndarray, B: float):
#    assert (len(vB)==3), "Error: Argument 'vB' must be a numpy.ndarray of length 3."
    assert (isinstance(vB, np.ndarray)), "Error: Argument 'vB' must be a numpy.ndarray."
    B = float(B)
    vB = vB.astype('float64')
    
    norm = np.linalg.norm(vB)
    
    for i in range(len(vB)):
#        print("B / norm * vB[i] = ", B, "/", norm, "*", vB[i], "=", B/norm*vB[i]) # DEBUG
        vB[i] = B/norm*vB[i]
        
    return vB


#%% ====== Transform EMF vector from compass frame to the body frame of the cage ======

# From [poppenk2007]
# X = +East, Y= +North, Z = +Azimuth
vEMF_compass = np.array([-0.1769, 19.0764, -44.8871])

# Orientation of TU Delft Helmholtz cage with respect to compass directions.
# Angle between North heading and axis of X coil pair (see Google maps)
a_LR = 23*np.pi/180

# Assume HHC setup at TU Delft is aligned with the Aerospace building walls
#  and that it is level with the local horizon.
Rz = np.array([[np.cos(-a_LR), -np.sin(-a_LR), 0],
               [np.sin(-a_LR),  np.cos(-a_LR), 0],
               [            0,              0, 1]])
vEMF = np.dot(Rz,vEMF_compass)

# Sanity check:
v1 = vEMF_compass[0:2]/np.linalg.norm(vEMF_compass[0:2])
v2 = vEMF[0:2]/np.linalg.norm(vEMF[0:2])
a_check = 180/np.pi*np.arccos(np.dot(v1,v2)) # if a_check == a_LR -> CHECK!


#%% ====== Define components ========

Vmax_supply = 150           # [V] Maximum supply voltage of maximum current
Imax_supply = 15            # [A] Maximum supply current
eff_supply = 0.99           # [-] Supply efficiency
tt_supply = 0.050           # [s] Transient time for 10%-90% power (=slew rate)
                            #       typical figure for lab supplies O(-2) s

supplyX = PowerSupply(Vmax_supply,Imax_supply,eff_supply, tt_supply)
supplyY = PowerSupply(Vmax_supply,Imax_supply,eff_supply, tt_supply)
supplyZ = PowerSupply(Vmax_supply,Imax_supply,eff_supply, tt_supply)

# Assuming 30 meters of AWG13 wires, plus four connectors per coil, with each
#   connector having a resistance of 20 mOhm. Each coil pair will get a set.
Rwiring = 30*0.00328 + 4*0.020
# Rwiring = 1 # Dummy testing for seeing effect

N_windings = 83
sidelength_coil = [1.85, 1.95, 2.05] # [m] X,Y,Z
width_coil = 0.05 # [m] Coil thickness estimated 
Rdl = 0.0062 # [Ohm/m] Impedance of coil wire (per meter)

coilX = HHCoil(sidelength_coil[0],width_coil,N_windings,Rdl,Rwiring,supplyX)
coilY = HHCoil(sidelength_coil[1],width_coil,N_windings,Rdl,Rwiring,supplyY)
coilZ = HHCoil(sidelength_coil[2],width_coil,N_windings,Rdl,Rwiring,supplyZ)

# Define cage object
cage = HHCage(coilX, coilY, coilZ, vEMF=vEMF)

#%% ====== Calculations =======

print("For a supply slew rate of", supplyX.t_transient*1000, "ms:")
print("VLt =", [round(coil.VLt_max(),3) for coil in cage.coils()], " [V]")


# Case: Cancel the EMF and create a net zero magnetic field vector
print("\n============ CASE: Cancel vEMF =============")
print(" -> generate a field of vB =", -1*cage.vEMF.round(3), "[uT]")
Ireq_EMF = cage.Ireq_breakdown(np.array([0,0,0]))
cage.properties_vB_req(np.array([0,0,0]), cancelEMF=True)


# Case: Cancel the EMF and generate 250 uT in +Z direction (Z-coil has highest 
#   impedance)
print("\n============ CASE: 250 uT in +Z ============")
vB_desired = scale_vector(np.array([0,0,1]),250)
print(" -> generate a field of vB =", vB_desired.round(3), "[uT]")
cage.properties_vB_req(vB_desired, cancelEMF=True)


# Case: Cancel the EMF and generate 250 uT diagonally (X, Y, Z coils all 
#   generate the same magnetic field strength)
print("\n======= CASE: 250 uT in XYZ diagonal =======")
vB_desired = scale_vector(np.array([1,1,1]),250)
print(" -> generate a field of vB =", vB_desired.round(3), "[uT]")
cage.properties_vB_req(vB_desired, cancelEMF=True)


# Case: Cancel the EMF and generate 250 uT exactly opposite to the EMF vector
print("\n======= CASE: 250 uT opposite to EMF =======")
vB_desired = scale_vector(-1*vEMF,250)
print(" -> generate a field of vB =", vB_desired.round(3), "[uT]")
cage.properties_vB_req(vB_desired, cancelEMF=True)


# Plot the characteristics of the cage
cage.plot_current_performance()

#%% Timing stuff
time2 = time.perf_counter()
print("\n  Time:", round(time2-time1,6), "[s]")