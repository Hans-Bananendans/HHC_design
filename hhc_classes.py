#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 12:37:03 2021

@author: Johan Monster
"""

import numpy as np
import matplotlib.pyplot as plt

class PowerSupply:
    def __init__(self, Vmax, Amax, efficiency, t_transient):
        self.Vmax = Vmax
        self.Amax = Amax
        self.eff = efficiency
        self.Vdotmax = 0.8*Vmax/t_transient # Slew rate for voltage
        self.Adotmax = 0.8*Amax/t_transient # Slew rate for current
        self.t_transient = t_transient
    
    def Pmax(self):
        return self.Amax*self.Vmax/self.eff
    
    def P(self, V, A):
        return self.A*self.V/self.eff

    def listmethods(self):
        """Function that lists all non-reserved methods in the class"""
        object_methods = [method_name for method_name in dir(self) if callable(getattr(self, method_name))]        
        for method in object_methods:
            if method[0:2] != "__":
                print(str(method+"()"))
    

class HHCoil:
    """
    Class for a pair of SQUARE Helmholtz coils
    """
    #                                lz
    #     __________               <---->                    
    #    |          |           ^  oooooo               oooooo
    #    |          |           |  ------               ------
    #    |          |lx      lx |            spacing            
    #    |          |           |        <------------->
    #    |__________|           |  ------               ------
    #        ly=lx              v  oooooo               oooooo <- N_windings 
    #                                    
    
    def __init__(self, 
                 lx,                    # [m] Side length of coil in X-direction
                 lz,                    # [m] Width of coil perpendicular to the windings
                 N_windings,            # [-] Number of windings of coil
                 Rdl,                   # [Ohm/m] Resistance of coil wire per meter
                 Rconnectors,           # [Ohm] Total resistance of connectors
                 supply: PowerSupply,   # PowerSupply object powering the coil
                 spacing_ratio=0.5445   # Ratio between coil spacing and lx
                 ):
        self.lx = lx            
        self.lz = lz            
        self.N = N_windings     
        self.Rdl = Rdl          
        self.Rconnectors = Rconnectors    
        
        # Crude way of validating that supply is the correct class:
        assert ("PowerSupply" in str(type(supply))), \
        "Error: argument 'supply' is not a 'PowerSupply' class object!"
        
        self.supply = supply
        
        self.mu0 = 1.256637062E-6   # [H/m] Vacuum permeability
        
        self.spacing = spacing_ratio*lx # [m] Distance between coil pair
        
        self.l_wire = self.l_wire()     # [m] Total length of BOTH coil wire
        self.R = self.l_wire*self.Rdl   # [Ohm] Total resistance of coil wire
        
        self.L = 2*self.L_1coil()       # [H] Inductance of BOTH coils

    def l_wire(self):
        """
        Calculates minimum wire length needed to make BOTH coils, assuming the 
            coils are wired in series. Includes the length of wire between 
            the coils twice, as it is assumed that the connection to the supply
            is located at one of the coils.
        """
        l_wire_1coil = 4*self.lx*self.N
        
        return 2*l_wire_1coil + 2*self.spacing # [m] 

    def L_1coil(self): 
        """
         Function from [wheeler1982] for computing the short-coil inductance 
          of a square-shaped coil.
           - lx is the side-length of the square coil in [m]
           - lz is the width of the coil in [m]
           - mu0 is the vacuum permeability of 1.256637062E-6 [H/m]
           - N is the number of windings in the coil [-]
           - L is the coil inductance in [H]
        """   
    
        a = self.lx/2
        b = self.lz
    
        # Calculate the short-coil, square coil inductance:
        L = self.mu0*self.N**2*a*4/np.pi * \
            (np.log(1 + np.pi * a/b) + 1/(3.644 + 2.0*b/a + 0.5098*(b/a)**2))
        
        return L # [H]
    
    def Rtot(self):
        return self.R + self.Rconnectors
    
    def VLt_max(self):
        """
        Calculates the voltage needed to overcome the inductance L of the coil
         pair at maximum slew rate dI/dt, according to Lenz's Law.
        This voltage is negative by default.
        """
        return -self.L * 0.8*self.supply.Adotmax
    
    def Vneeded(self, I):
        """
        Calculate voltage needed to maintain a certain current through the 
         coils.
        """
        return I*(self.R + self.Rconnectors)
    
    def Pneeded(self, I):
        """
        Calculate power needed to maintain a certain current through the 
         coils. This ignores the power that the coil needs to push against
         the Earth Magnetic Field (EMF)
        """
        return self.Vneeded(I)*I
        
    
    def print_summary(self):
        """
        Prints a summary of the coil characteristics
        """
        for key, val in self.__dict__.items():
            if val is not None:
                print("%s: %s" % (key, val))
    
    def compute_Q(self, coor):
        """
        Compute "Q-term", a mathematical combination of several coil parameters
         that is need to apply the Biot-Savart law for square Helmholtz coils.
         For details, see [frix1994].
        """
        [x,y,z] = coor
        
        a = self.lx/2
        xn = [a, -a]
        yn = [a, -a]
        zn = [self.spacing/2, -self.spacing/2] # Vertical positions of the coils
        
        Q = 0
        for i in range(0,2):
            for j in range(0,2):
                for k in range(0,2):
                    # Split up long equation into several parts for readability:
                    q1 = (-1)**((i+1)+(j+1))
                    q2 = (x-xn[i])*(y-yn[j]) / ((x-xn[i])**2+(y-yn[j])**2+(z-zn[k])**2)**0.5
                    q3 = ( 1/((x-xn[i])**2+(z-zn[k])**2) + 1/((y-yn[j])**2+(z-zn[k])**2) )
                    Q += q1*q2*q3
        return Q
        
        
    def calc_B_HHC(self, coor, I): 
        """
         Function from [frix1994] for computing the field strength at a point 
           [x,y,z] located in between a SQUARE Helmholtz coil pair.
           - coor is a 3D Cartesian coordinate in [m]
           - coil_side is the side-length of the square coil in [m]
           - coil_spacing is the distance between the coils
           - mu0 is the vacuum permeability of 1.256637062E-6 [H/m]
           - n is the number of windings in the coil [-]
           - I is the coil current in [A]
           - B is the field strength at the coordinate in [uT] (default) or [T]
        """
        # Calculate the field strength at point (Biot-Savart)
        Bz = self.mu0 * self.N * I/(4*np.pi) * self.compute_Q(coor)
        
        # Convert from [T] to [uT]
        Bz *= 1E6
            
        return Bz # [uT]
    
    def calc_Bmid(self, I):
        """
        Calculates Bmid, the field strength in [uT] (default) or [T] in the
         exact geometric middle of the Helmholtz coil pair. This point lies
         on the commom axis of the two coils, and at 1/2 the coil spacing.
        """
        return self.calc_B_HHC([0,0,0.5*self.spacing], I)
        
    def calc_Bmid_max(self):
        """
        Calculates Bmid, the field strength in [uT] (default) or [T] in the
         exact geometric middle of the Helmholtz coil pair. This point lies
         on the commom axis of the two coils, and at 1/2 the coil spacing.
        """
        return self.calc_B_HHC([0,0,0.5*self.spacing], self.supply.Amax)
    
    def calc_Ireq_coor(self, coor, Breq): 
        """
         Rewritten equation from [frix1994] for computing the field strength at 
           a point [x,y,z] located in between a SQUARE Helmholtz coil pair.
           - coor is a 3D Cartesian coordinate in [m]
           - coil_side is the side-length of the square coil in [m]
           - coil_spacing is the distance between the coils
           - mu0 is the vacuum permeability of 1.256637062E-6 [H/m]
           - n is the number of windings in the coil [-]
           - Breq is the required field strength in [T] or [uT] at coor
           - Ireq is the coil current in [A] required for a certain Breq
        """
        # Convert from [uT] to [T]
        Breq *= 1E-6
    
        # Calculate the field strength at point (Biot-Savart)
        Ireq = 1/self.N * (4*np.pi)/self.mu0 * Breq/self.compute_Q(coor)
        
        return Ireq # [A]

    def calc_Ireq(self, Breq): 
        """
         Rewritten equation from [frix1994] for computing the field strength 
          in the middle of a SQUARE Helmholtz coil pair.
           - coor is a 3D Cartesian coordinate in [m]
           - coil_side is the side-length of the square coil in [m]
           - coil_spacing is the distance between the coils
           - mu0 is the vacuum permeability of 1.256637062E-6 [H/m]
           - n is the number of windings in the coil [-]
           - Breq is the required field strength in [T] or [uT] at coor
           - Ireq is the coil current in [A] required for a certain Breq
        """
        # Convert from [uT] to [T]
        Breq *= 1E-6
    
        # Define coordinate of middle point on the coil axis
        coor = [0,0,0.5*self.spacing]
        
        # Calculate the field strength at point (Biot-Savart)
        Ireq = 1/self.N * (4*np.pi)/self.mu0 * Breq/self.compute_Q(coor)
        
        return Ireq # [A]
    
    def plot_current_performance(self, Arange='oops'):
        """
        Plot the Vneeded, Pneeded, and Bmid as function of a range of values
         for the current (Arange), and plot these.
        """
        if isinstance(Arange, str):
            np.linspace(0,self.supply.Amax,101)
        
        Vrange = []
        Prange = []
        Brange = []
        
        for i in range(len(Arange)):
            Vrange.append(self.Vneeded(Arange[i]))
            Prange.append(self.Pneeded(Arange[i]))
            Brange.append(self.calc_Bmid(Arange[i]))
        
        fig, axs = plt.subplots(3,1)
        
        varinst = [Vrange, Prange, Brange]
        varnames = ["Voltage [V]", "Power [W]","Field strength [uT]"]
        for i in range(len(varinst)):
            axs[i].plot(Arange,varinst[i],'r')
            axs[i].set_xlim(min(Arange),max(Arange))
            axs[i].set_xlabel("Current [A]")
            axs[i].set_ylabel(varnames[i])
            axs[i].grid(True)
            
        fig.suptitle("Voltage, Consumed power, and magnetic fieldstrength in \n the middle of the pair, all as function of current")
        fig.tight_layout()
        plt.show
        
    def T_equilibrium(self, A, l1, l2, t1, t2, dwire):
        """
        Estimated from pixel counting pictures in paper:
            l1 ~= 69 mm
            l2 ~= 28 mm
        Estimated from memory:
            t1 = t2 = 3 mm
        
        Calculate the equilibrium temperature of the coil for a given current. 
         Coil cross section is defined as follows:
        """
        #          t2
        #         <->
        #          ___________________
        #         |___________________| t1          
        #      ^  | |oooooooooooooo                 
        #      |  | |oooooooooooooo                 
        #      |  | |oooooooooooooo <--- dwire        
        #   l2 |  | |oooooooooooooo                 
        #      |  | |oooooooooooooo
        #      |  | |oooooooooooooo
        #      v  |_|oooooooooooooo___
        #         |___________________|
        #         
        #         <------------------->
        #                  l1
        pass
        # TODO: Implement this
    
    def listmethods(self):
        """Function that lists all non-reserved methods in the class"""
        object_methods = [method_name for method_name in dir(self) if callable(getattr(self, method_name))]        
        for method in object_methods:
            if method[0:2] != "__":
                print(str(method+"()"))


class HHCage: 
    def __init__(self, 
                 Xcoil: HHCoil,
                 Ycoil: HHCoil,
                 Zcoil: HHCoil,
                 vEMF=np.array([0,0,0])
                 ):
        
        self.Xcoil = Xcoil
        self.Ycoil = Ycoil
        self.Zcoil = Zcoil
        self.vB = np.array([0,0,0])
        self.vEMF = vEMF
        
        for obj in self.coils():
            assert ("HHCoil" in str(type(obj))), \
            "Error: argument 'supply' is not a 'HHCoil' class object!"
        
    def coils(self):
        return [self.Xcoil, self.Ycoil, self.Zcoil]
    
    def Bmax(self):
        """
        Calculate the maximum field strength that can be generated in any 
         direction, which corresponds to the lowest maximum field strength of 
         the three coils.
        """
        B = []
        for coil in self.coils():
            B.append(coil.calc_Bmid_max())
        return min(B)
    
    def vB_breakdown(self):
        """
        Breaks down a desired magnetic field vector vB into magnetic field
         vectors for each individual coil.
        """
        return [np.array([self.vB[0],          0,          0]), 
                np.array([         0, self.vB[0],          0]), 
                np.array([         0,          0, self.vB[0]])]
        
    def vB_len(self):
        """Returns the length of vB"""
        return np.linalg.norm(self.vB)
        
    def Ireq_breakdown(self, vB, cancelEMF=False):
        """
        Breaks down a desired magnetic field vector vB into required coil
        currents for each individual coil.
        """
        Ireq = []
        for i, coil in enumerate(self.coils()):
            Ireq.append(coil.calc_Ireq(vB[i]))
        return np.array(Ireq)

    def Ireq_cancelEMF(self):
        """
        Returns required coil currents for each individual coil to cancel out
         the Earth Magnetic Field (EMF) at the Helmholtz Coil Site
        """
        Ireq = []
        for i, coil in enumerate(self.coils()):
            Ireq.append(coil.calc_Ireq(-self.vEMF[i]))
        return np.array(Ireq)
            
    
#    def vEMF(self, vEMF: np.ndarray):
#        assert (len(vEMF)==3), "Error: Argument 'vEMF' must be a numpy.ndarray of length 3."
#        assert (isinstance(vEMF, np.ndarray)), "Error: Argument 'vEMF' must be a numpy.ndarray."
#        self.vEMF = vEMF
    
    def properties_vB_req(self, vB: np.ndarray, cancelEMF=False):
        """
        Compute Vneeded, Pneeded, and Bmid as function of the current needed
         for a desired magnetic field vector.
        """
        
        assert (len(vB)==3), "Error: Argument 'vB' must be a numpy.ndarray of length 3."
        assert (isinstance(vB, np.ndarray)), "Error: Argument 'vB' must be a numpy.ndarray."

        if cancelEMF:
            Ireq = self.Ireq_breakdown(vB)+self.Ireq_cancelEMF()
        else:
            Ireq = self.Ireq_breakdown(vB)
        
        Vreq = []
        Preq = []
        Bresult = []
        
        for i, coil in enumerate(self.coils()):
            Vreq.append(coil.Vneeded(Ireq[i]))
            Preq.append(coil.Pneeded(Ireq[i]))
            Bresult.append(coil.calc_Bmid(Ireq[i]))
        
        print("Ireq =", np.array(Ireq).round(3), "[A]")
        print("Vreq =", np.array(Vreq).round(3), "[V]")
        print("Preq =", np.array(Preq).round(3), "[W]")
        print("Ptotal =", round(sum(Preq),3), "[W]")
        print("Bresult =", np.array(Bresult).round(3), "[uT]")
        print("vEMF =", self.vEMF.round(3), "[uT]")
        print("Resulting field: ", (np.array(Bresult)+self.vEMF).round(3), "[uT]")

    def properties_cancelEMF(self):
        """
        Compute Vneeded, Pneeded, and Bmid as function of the current needed
         to cancel the local Earth Magnetic Field at the cage.
        """

        Ireq = self.Ireq_cancelEMF()
        Vreq = []
        Preq = []
        Bresult = []
        
        for i, coil in enumerate(self.coils()):
            Vreq.append(coil.Vneeded(Ireq[i]))
            Preq.append(coil.Pneeded(Ireq[i]))
            Bresult.append(coil.calc_Bmid(Ireq[i]))
        
        print("Ireq =", np.array(Ireq).round(3), "[A]")
        print("Vreq =", np.array(Vreq).round(3), "[V]")
        print("Preq =", np.array(Preq).round(3), "[W]") 
        print("Bresult =", np.array(Bresult).round(3), "[uT]")
        print("vEMF =", self.vEMF.round(3), "[uT]")
        print("Resulting field: ", (np.array(Bresult)+self.vEMF).round(3), "[uT]")
        
    def listmethods(self):
        """Function that lists all non-reserved methods in the class"""
        object_methods = [method_name for method_name in dir(self) if callable(getattr(self, method_name))]        
        for method in object_methods:
            if method[0:2] != "__":
                print(str(method+"()"))
                

class ThermalModel:
    def __init__(self):  
        pass

    
    def listmethods(self):
        """Function that lists all non-reserved methods in the class"""
        object_methods = [method_name for method_name in dir(self) if callable(getattr(self, method_name))]        
        for method in object_methods:
            if method[0:2] != "__":
                print(str(method+"()"))
    
    
    def T_equilibrium(self, A, l1, l2, t1, t2, dwire):
        """
        Estimated from pixel counting pictures in paper:
            l1 ~= 69 mm
            l2 ~= 28 mm
        Estimated from memory:
            t1 = t2 = 3 mm
        
        Calculate the equilibrium temperature of the coil for a given current. 
         Coil cross section is defined as follows:
        """
        #          t2
        #         <->
        #          ___________________
        #         |___________________| t1          
        #      ^  | |oooooooooooooo                 
        #      |  | |oooooooooooooo                 
        #      |  | |oooooooooooooo <--- dwire        
        #   l2 |  | |oooooooooooooo                 
        #      |  | |oooooooooooooo
        #      |  | |oooooooooooooo
        #      v  |_|oooooooooooooo___
        #         |___________________|
        #         
        #         <------------------->
        #                  l1
        pass
        # TODO: Implement this