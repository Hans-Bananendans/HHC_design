from hhc_classes import PowerSupply, HHCoil
import numpy as np
import matplotlib.pyplot as plt

global DEBUG
DEBUG = 1

#%% Coil parameters
N_windings = 83
sidelength_coil = [1.85, 1.95, 2.05] # [m] X,Y,Z
width_coil = 0.019 # [m] Coil thickness (from Poppenk2009)
Rcoils = [6.27, 6.85, 6.92]  # Ohm

supply = PowerSupply(60, 5, 0.86, 1)

coilX = HHCoil(sidelength_coil[0],
               width_coil,
               N_windings,
               0,
               0,
               supply
               )
coilX.set_R(Rcoils[0])

# def get_Vx_ref(Bref, Q, n=83, slope=1, offset=0):
#     return (2/5E-7) * Bref / (n * Q)

def B2V(Bref, Q, n=83, slope=1, offset=0):
    """
    Translates magnetic flux density in T to supply control voltage range
    """    
    return (2/5E-7) * Bref / (n * Q)


def PID(dt, K, ut_1, et, et_1, et_2):
    """
    Computes PID control signal
    """
    [KP, KI, KD] = K
    ut = ut_1
    + (KP + KI*dt + KD/dt) * et
    + (-1*KP - 2*KD/dt) * et_1
    + KD / dt * et_2
    return ut


def compute_Ic(Vc):
    # Returns a positive, linear supply current according to:
    # I_supply = 5/2 * V_in  | Also automatically caps at 5A
    Ic_abs = min(5/2 * abs(Vc), 5)
    
    # Returns the sign of Ic separately
    if Vc < 0:
        sign_I = -1
    else:
        sign_I = 1
    
    return Ic_abs, sign_I
    
# def supply(V_in):
#     """
#     Returns a positive, linear supply current according to:
#     I_supply = 5/2 * V_in
#     Automatically caps at 5A
#     """
#     return min(5/2 * abs(V_in), 5)

# def H_Bridge(V_in):
#     """
#     Extracts the sign of the input voltage (-1 or 1)
#     """
#     return V_in/abs(V_in)

# def fI_coil(I_in, sign):
#     """
#     Signs the supply current.
#     """
#     return I_in * sign

def B_coil(I_in, Q, n=83):
    return 1E-7 * n * I_in / Q

def V_sensor(B):
    """
    Outputs a voltage proportional to the magnetic field.
    If the value is out of bound of the sensors +/- 2000 mG, or 200 uT,
    it will return the saturation voltage.
    """
    if B > 200E-6:
        return 2.0
    if B < -200-6:
        return -2.0
    else:
        return B/100*1E6

def B_sensor(V_sensor):
    """
    Converts readings of the magnetometer to flux density values [T]
    """
    return V_sensor * 1E-6 * 100

#%% Simulation configuration

t0 = 0
tmax = 1
dt = 0.1
t=0
i=0

Bx = 0
Bx_meas = 0
Bx_ref = 175E-6
Vc = 0

# Gains
[KP, KI, KD] = [10, 10, 1]
K = [KP, KI, KD]

# Pre-allocate data vectors
ldata = int((tmax-t0)/dt+1)

data_i = np.zeros(ldata)
data_t = np.zeros(ldata)
data_Bx = np.zeros(ldata)
data_Bx_meas = np.zeros(ldata)
data_Bx_ref = np.zeros(ldata)
data_Vc = np.zeros(ldata)
data_e = np.zeros(ldata)


if DEBUG == 1:
    print("|   i   |  Bx_ref  |   Bx   |    e    |    Vc   |")
    print("|-------|----------|--------|---------|---------|")

#%% Loop
while t <= tmax:        

    # Convert Bref to Vref
    Vref = B2V(Bx_ref, coilX.Q)
    
    # Startup conditionals
    if i == 0:
        ut_1 = 0
        et_1 = 0
        et_2 = 0
        Bmeas = 0
        Vmeas = B2V(Bmeas, coilX.Q)
        et = Vref - Vmeas
    else:
        # Convert Bmeas to Vmeas
        Vmeas = B2V(Bmeas, coilX.Q)
        
        ut_1 = data_Vc[i-1]
        
        et = Vref - Vmeas
        et_1 = data_e[i-1]
        if i == 1:
            et_2 = 0
        else:
            et_2 = data_e[i-2]
            
    
    # PID control voltage:
    Vc = PID(dt, K, ut_1, et, et_1, et_2)
    
    # Calculate the power coming out of the power supplies at Vc
    Ic_abs, sign_I = compute_Ic(Vc)
    
    Bx = B_coil(Ic_abs*sign_I, coilX.Q)
    
    Bmeas = B_sensor(V_sensor(Bx))
    
    # Store data
    data_t[i] = t
    data_t[i] = i
    data_Bx[i] = Bx_meas
    data_Bx_meas[i] = Bx_meas
    data_Bx_ref[i] = Bx_ref
    data_Vc[i] = Vc
    data_e[i] = et
    
    if DEBUG == 1:
        print("|", str(i).rjust(5), 
              "|", str(round(Bx_ref, 6)).rjust(6),
              "|", str(round(Bx, 6)).rjust(6), 
              "|", str(round(et, 3)).rjust(7), 
              "|", str(round(Vc, 3)).rjust(7), "|")
    # Increment
    i += 1
    t += dt
    
    
#%% Plotting

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13,7))
ax1.axhline(0, color='k')
ax2.axhline(0, color='k')

# B plots
ax1.plot(data_t[0:-1], data_Bx_ref[0:-1], c="cyan", label="B_ref")
ax1.plot(data_t[0:-1], data_Bx[0:-1], "#C00", label="B_coil")

# V plots
ax2.plot(data_t[0:-1], data_Vc[0:-1], "#00F", label="V_control")

# Labels
ax1.set_xlabel("Time [s]")
ax1.set_ylabel("Magnetic field [T]")
ax1.legend(loc="lower right")
plot_title="Control envelope"

fig.suptitle(plot_title)