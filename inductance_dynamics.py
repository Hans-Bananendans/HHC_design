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

# From IGRF2020 and WMM-2020
# X = +East, Y= +North, Z = +Azimuth
vEMF_compass = np.array([0.72645, 19.17495, -45.45745])

# Orientation of TU Delft Helmholtz cage with respect to compass directions.
# Angle between North heading and axis of X coil pair (see Google maps)
a_LR = 23*np.pi/180

# Assume HHC setup at TU Delft is aligned with the Aerospace building walls
#  and that it is level with the local horizon.
Rz = np.array([[np.cos(-a_LR), -np.sin(-a_LR), 0],
               [np.sin(-a_LR),  np.cos(-a_LR), 0],
               [            0,              0, 1]])
vEMF = np.dot(Rz,vEMF_compass)


#%% ====== Define components ========

Vmax_supply = 30            # [V] Maximum supply voltage of maximum current
Imax_supply = 10            # [A] Maximum supply current
eff_supply = 0.865          # [-] Supply efficiency
tt_supply = 0.050           # [s] Transient time for 10%-90% power (=slew rate)
                            #       typical figure for lab supplies O(-2) s

supplyX = PowerSupply(Vmax_supply,Imax_supply,eff_supply, tt_supply)
supplyY = PowerSupply(Vmax_supply,Imax_supply,eff_supply, tt_supply)
supplyZ = PowerSupply(Vmax_supply,Imax_supply,eff_supply, tt_supply)

# Resistance has been measured:
Rcoils = [6.27, 6.85, 6.92] # Ohm

N_windings = 83
sidelength_coil = [1.85, 1.95, 2.05] # [m] X,Y,Z
# width_coil = 0.05 # [m] Coil thickness estimated
width_coil = 0.019 # [m] Coil thickness (from Poppenk2009)
Rdl = 0.00657 # [Ohm/m] Impedance of coil wire (per meter)

coilX = HHCoil(sidelength_coil[0],width_coil,N_windings,0,Rdl,supplyX)
coilY = HHCoil(sidelength_coil[1],width_coil,N_windings,0,Rdl,supplyY)
coilZ = HHCoil(sidelength_coil[2],width_coil,N_windings,0,Rdl,supplyZ)
coilX.set_R(Rcoils[0])
coilY.set_R(Rcoils[1])
coilZ.set_R(Rcoils[2])

# Define cage object
cage = HHCage(coilX, coilY, coilZ, vEMF=vEMF)

#%% ====== Calculations =======

# Time steps
bmax = 241
rpm = 300
period = 60/rpm
t = np.linspace(0, period, 100)

def bxt(t_array, bmax, period):
    bx = np.zeros(len(t_array))
    ix = np.zeros(len(t_array))
    dix = np.zeros(len(t_array))
    
    imax = cage.Xcoil.calc_Ireq(bmax)
    for i in range(len(t_array)):
        bx[i] = bmax * np.sin(2*np.pi / period * t_array[i])
        ix[i] = imax * np.sin(2*np.pi / period * t_array[i])
        dix[i] = imax * 2*np.pi / period * np.cos(2*np.pi / period * t_array[i])
    return bx, ix, dix


def bzt(t_array, bmax, period):
    bz = np.zeros(len(t_array))
    iz = np.zeros(len(t_array))
    diz = np.zeros(len(t_array))
    
    imax = cage.Zcoil.calc_Ireq(bmax)
    for i in range(len(t_array)):
        bz[i] = bmax * np.cos(2*np.pi / period * t_array[i])
        iz[i] = imax * np.cos(2*np.pi / period * t_array[i])
        diz[i] = -imax * 2*np.pi / period * np.sin(2*np.pi / period * t_array[i])
    return bz, iz, diz

# -bmax for x-axis to flip sinusoid (which changes direction of rotation)
bx, ix, dix = bxt(t, -bmax, period)
bz, iz, diz = bzt(t, bmax, period)

# Voltage due to supply of current for magnetic field
vx_field = np.array([cage.Xcoil.Vneeded(i) for i in ix])
# Voltage to overcome inductive voltage when field changes
vx_induction = np.array([cage.Xcoil.L*di for di in dix])

vz_field = np.array([cage.Zcoil.Vneeded(i) for i in iz])
vz_induction = np.array([cage.Zcoil.L*di for di in diz])

vx_total = vx_field + vx_induction
vz_total = vz_field + vz_induction


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13,7))
ax1.axhline(0, color='k')
ax2.axhline(0, color='k')

# X-coil plots
ax1.plot(t, vx_field, c="r", label="Generative component")
ax1.plot(t, vx_induction, "#FA0", label="Inductive component")
ax1.plot(t, vx_total, "k", label="Total")
# Dotted lines
vxmf = max(vx_field)
vxmt = max(vx_total)
ax1.plot([t[0], (t[-1]-t[0])/2, t[-1]], [vxmf, vxmf, vxmf],
          c="r", linestyle="dotted")
ax1.plot([t[0], (t[-1]-t[0])/2, t[-1]], [vxmt, vxmt, vxmt], 
          c="k",  linestyle="dotted")

# Z-coil plots
ax2.plot(t, vz_field, c="b", label="Generative component")
ax2.plot(t, vz_induction, "#FA0", label="Inductive component")
ax2.plot(t, vz_total, "k", label="Total")
# Dotted lines
vzmf = max(vz_field)
vzmt = max(vz_total)

print("Max baseline:",vzmf,"V")
print("Max total:",vzmt,"V")

ax2.plot([t[0], (t[-1]-t[0])/2, t[-1]], [vzmf, vzmf, vzmf], 
          "b", linestyle="dotted")
ax2.plot([t[0], (t[-1]-t[0])/2, t[-1]], [vzmt, vzmt, vzmt], 
          "k",  linestyle="dotted")

# Labels
ax1.set_xlabel("Time [s]")
ax2.set_xlabel("Time [s]")
ax1.set_ylabel("Supplied voltage [V]")

ax1.text((t[-1]-t[0])/2, 1.03*vxmt, "peak +{}%".format(round(100*(vxmt/vxmf-1)-1, 2)),
         fontsize=13, c="k")
ax2.text((t[-1]-t[0])/2, 1.03*vzmt, "peak +{}%".format(round(100*(vzmt/vzmf-1)-1, 2)),
         fontsize=13, c="k")

ax1.legend(loc="lower right")
ax2.legend(loc="lower right")

plot_title="""
Supplied voltage envelope of X and Z during the rotation ({} rpm, {} $\mu T$).
""".format(rpm, bmax)

fig.suptitle(plot_title)

filename = "./figures/inductance_dynamics_{}rpm.png".format(rpm)
fig.savefig(filename, dpi=150)
#print("R", cage.Xcoil.R, cage.Ycoil.R, cage.Zcoil.R)
#print("L", round(1000*cage.Xcoil.L, 1), round(1000*cage.Ycoil.L, 1), round(1000*cage.Zcoil.L, 1))


#%% Timing stuff
time2 = time.perf_counter()
print("\n  Time:", round(time2-time1,6), "[s]")