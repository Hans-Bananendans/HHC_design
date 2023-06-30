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
    
    def Rtot(self, Teq=20):
        # Use linearized relationship (constant dR/dT) to find extra resistance
        # caused by elevated temperatures of copper (see [ladino2015])
        # Assumption: Every conductor (including connectors) are purely copper
        
        alpha_copper = 0.00427 # [Ohm/K] see [ladino2015]
        
        return (self.R + self.Rconnectors)*(1+alpha_copper*(Teq-20))
    
    def VLt_max(self):
        """
        Calculates the voltage needed to overcome the inductance L of the coil
         pair at maximum slew rate dI/dt, according to Lenz's Law.
        This voltage is negative by default.
        """
        return -self.L * 0.8*self.supply.Adotmax
    
    def Vneeded(self, I, Teq=20):
        """
        Calculate voltage needed to maintain a certain current through the 
         coils.
        """
        return I*self.Rtot(Teq)
    
    def Pneeded(self, I, Teq=20):
        """
        Calculate power needed to maintain a certain current through the 
         coils. This ignores the power that the coil needs to push against
         the Earth Magnetic Field (EMF)
        """
        return self.Vneeded(I, Teq)*I
        
    
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
         Function (7) from [frix1994] for computing the field strength at a
           point [x,y,z] located in between a SQUARE Helmholtz coil pair.
           - coor is a 3D Cartesian coordinate in [m]
           - coil_side is the side-length of the square coil in [m]
           - coil_spacing is the distance between the coils
           - mu0 is the vacuum permeability of 1.256637062E-6 [H/m]
           - n is the number of windings in the coil [-]
           - I is the coil current in [A]
           - B is the field strength at the coordinate in [uT] (default) or [T]
        
        Note: the Cartesian coordinate system is centered on the AXIS LINE of
            the Helmholtz coil pairs at the exact MIDDLE OF THE COILS. In other
            words, the coordinate system is located exactly in the middle.
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
        # return self.calc_B_HHC([0,0,0.5*self.spacing], I)
        return self.calc_B_HHC([0,0,0], I)
        
    def calc_Bmid_max(self):
        """
        Calculates Bmid, the field strength in [uT] (default) or [T] in the
         exact geometric middle of the Helmholtz coil pair. This point lies
         on the commom axis of the two coils, and at 1/2 the coil spacing.
        """
        # return self.calc_B_HHC([0,0,0.5*self.spacing], self.supply.Amax)
        return self.calc_B_HHC([0,0,0], self.supply.Amax)
    
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
        # coor = [0,0,0.5*self.spacing]
        coor = [0,0,0]
        
        # Calculate the field strength at point (Biot-Savart)
        Ireq = 1/self.N * (4*np.pi)/self.mu0 * Breq/self.compute_Q(coor)
        
        return Ireq # [A]
 
    """Unused method - Due removal """
#    def plot_current_performance(self, Arange='nothing', Teq=20):
#        """
#        Plot the Vneeded, Pneeded, and Bmid as function of a range of values
#         for the current (Arange), and plot these.
#        """
#        if isinstance(Arange, str):
#            np.linspace(0,self.supply.Amax,101)
#        
#        Vrange = []
#        Prange = []
#        Brange = []
#        
#        for i in range(len(Arange)):
#            Vrange.append(self.Vneeded(Arange[i]))
#            Prange.append(self.Pneeded(Arange[i]))
#            Brange.append(self.calc_Bmid(Arange[i]))
#        
#        fig, axs = plt.subplots(3,1)
#        
#        varinst = [Vrange, Prange, Brange]
#        varnames = ["Voltage [V]", "Power [W]","Field strength [uT]"]
#        for i in range(len(varinst)):
#            axs[i].plot(Arange,varinst[i],'r')
#            axs[i].set_xlim(min(Arange),max(Arange))
#            axs[i].set_xlabel("Current [A]")
#            axs[i].set_ylabel(varnames[i])
#            axs[i].grid(True)
#            
#        fig.suptitle("Voltage, Consumed power, and magnetic fieldstrength in \n the middle of the pair, all as function of current")
#        fig.tight_layout()
#        plt.show

    
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
            
    """Unused method - Due removal """
#    def vEMF(self, vEMF: np.ndarray):
#        assert (len(vEMF)==3), "Error: Argument 'vEMF' must be a numpy.ndarray of length 3."
#        assert (isinstance(vEMF, np.ndarray)), "Error: Argument 'vEMF' must be a numpy.ndarray."
#        self.vEMF = vEMF
    
    def properties_vB_req(self, vB: np.ndarray, cancelEMF=False, Teq=20):
        """
        Compute Vneeded, Pneeded, and Bmid as function of the current needed
         for a desired magnetic field vector.
        Use Teq to set the expected equilibrium temperature of the coils (default: 20 C)
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
            Vreq.append(coil.Vneeded(Ireq[i], Teq=Teq))
            Preq.append(coil.Pneeded(Ireq[i], Teq=Teq))
            Bresult.append(coil.calc_Bmid(Ireq[i]))
        
        print("Ireq =", np.array(Ireq).round(3), "[A]")
        print("Vreq =", np.array(Vreq).round(3), "[V]")
        print("Preq =", np.array(Preq).round(3), "[W]")
        print("Ptotal =", round(sum(Preq),3), "[W]")
        print("Bresult =", np.array(Bresult).round(3), "[uT]")
        print("vEMF =", self.vEMF.round(3), "[uT]")
        print("Resulting field: ", (np.array(Bresult)+self.vEMF).round(3), "[uT]")

    """Unused method - Due removal """
#    def properties_cancelEMF(self, Teq):
#        """
#        Compute Vneeded, Pneeded, and Bmid as function of the current needed
#         to cancel the local Earth Magnetic Field at the cage.
#        """
#
#        Ireq = self.Ireq_cancelEMF()
#        Vreq = []
#        Preq = []
#        Bresult = []
#        
#        for i, coil in enumerate(self.coils()):
#            Vreq.append(coil.Vneeded(Ireq[i]), Teq=Teq)
#            Preq.append(coil.Pneeded(Ireq[i]), Teq=Teq)
#            Bresult.append(coil.calc_Bmid(Ireq[i]))
#        
#        print("Ireq =", np.array(Ireq).round(3), "[A]")
#        print("Vreq =", np.array(Vreq).round(3), "[V]")
#        print("Preq =", np.array(Preq).round(3), "[W]") 
#        print("Bresult =", np.array(Bresult).round(3), "[uT]")
#        print("vEMF =", self.vEMF.round(3), "[uT]")
#        print("Resulting field: ", (np.array(Bresult)+self.vEMF).round(3), "[uT]")
        
    def listmethods(self):
        """Function that lists all non-reserved methods in the class"""
        object_methods = [method_name for method_name in dir(self) if callable(getattr(self, method_name))]        
        for method in object_methods:
            if method[0:2] != "__":
                print(str(method+"()"))
                
    def plot_current_performance(self, Arange='nothing', Teq=20):
        """
        Plot the Vneeded, Pneeded, and Bmid as function of a range of values
          for the current (Arange), and plot these.
        Use Teq to set the expected equilibrium temperature of the coils (default: 20 C)
        """
        # If no current range is given, use a default instead
        if isinstance(Arange, str):
            defaultlen = 101
            Amap = [np.linspace(0,self.coils()[0].supply.Amax,defaultlen),
                    np.linspace(0,self.coils()[1].supply.Amax,defaultlen),
                    np.linspace(0,self.coils()[2].supply.Amax,defaultlen)]
        # Else, apply the given current range to all three coils
        else:
            Amap = [np.linspace(0,max(Arange),len(Arange)),
                    np.linspace(0,max(Arange),len(Arange)),
                    np.linspace(0,max(Arange),len(Arange))]     
        
        
        # Initialize dataset containers
        Vrange = [[],[],[]]
        Prange = [[],[],[]]
        Brange = [[],[],[]]
        
        for i, coil in enumerate(self.coils()):
            for j in range(len(Amap[i])):
                Vrange[i].append(coil.Vneeded(Amap[i][j], Teq=Teq))
                Prange[i].append(coil.Pneeded(Amap[i][j], Teq=Teq))
                Brange[i].append(coil.calc_Bmid(Amap[i][j]))
        
        fig, axs = plt.subplots(3,1)
        
        varinst = [Vrange, Prange, Brange]

        varnames = ["Voltage [V]", "Power [W]", "Field strength [uT]"]
        # For every dataset in varinst
        for i in range(len(varinst)):
            # Plot the dataset for the Xcoil in red, Y in green, Z in blue
            axs[i].plot(Amap[0],varinst[i][0],'r', label = "Xcoil")
            axs[i].plot(Amap[1],varinst[i][1],'g', label = "Ycoil")
            axs[i].plot(Amap[2],varinst[i][2],'b', label = "Zcoil")

            # Find the smallest and largest currents in the whole set Amap, 
            #   and set it as the plot boundaries
            axs[i].set_xlim(
                    min([min(A) for A in Amap]),
                    max([max(A) for A in Amap]))
            
            # Label the axes
            axs[i].set_xlabel("Current [A]")
            axs[i].set_ylabel(varnames[i])
            
            # Put a grid in
            axs[i].grid(True)
            
            # Put a legend in
            axs[i].legend(loc="upper left")
            
        fig.suptitle("Voltage, Consumed power, and magnetic fieldstrength in \n the middle of the pair, all as function of current")
        fig.tight_layout()
        plt.show
        
        return varinst, Amap[0]
        
    def plot_dRdT(self, Arange='nothing'):
        """
        Plot the R-T relationship as function of the equilibrium temperature
          Teq for each coil.
        """
        # If no current range is given, use a default instead (0 C - 120 C)
        if isinstance(Arange, str):
            defaultlen = 101
            Amap = np.linspace(0, 120 ,defaultlen)
        # Else, apply the given temperature range to all three coils
        else:
            Amap = np.linspace(min(Arange),max(Arange),len(Arange))
        
        Rrange = [[],[],[]]
        
        for i, coil in enumerate(self.coils()):
            for j in range(len(Amap)):
                Rrange[i].append(coil.Rtot(Teq=Amap[j]))
        
        fig, ax = plt.subplots()
        
        ax.plot(Amap, Rrange[0], 'r', label = "Xcoil")
        ax.plot(Amap, Rrange[1], 'g', label = "Ycoil")
        ax.plot(Amap, Rrange[2], 'b', label = "Zcoil")
        
        ax.set_xlim(min(Amap), max(Amap))
                    
        ax.set_title("Total coil resistance as function of temperature")
        ax.set_xlabel("Conductor temperature [deg C]")
        ax.set_ylabel("Total resistance [Ohm]")
        
        # Put a grid in
        ax.grid(True)
        
        # Put a legend in
        ax.legend(loc="upper left")
        
        plt.show

class ThermalModel:
    """
    Estimated from pixel counting pictures in paper:
        l1 ~= 35mm
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
    
    """ Currently implemented model: """
    
    #       t2 = 3mm
    #         <->
    #          __________________
    #         |__________________| t1 = 3mm          
    #      ^  | |              | |
    #      |  | |          T1->| |<-T2      T_amb   
    #      |  | |              | |
    # l2 = |  | |  q_current   | |---> q_radiation
    # 22mm |  | |              | |---> q_convection
    #      |  | |             -|-|-> q_conduction
    #      v  |_|______________|_|
    #         |__________________|
    #         
    #         <------------------->
        #               l1 = 28mm      

    """ Fair warning: results as of now (01-12-2021), the methods in this class
        are NOT validated with reality, and may give poor results.
    """
    # TODO: Generalize constructor
    def __init__(self):  
        self.N = 83
        self.k = 237 # [W/mK] k_aluminium
        self.h = 10 #[W/m2k] Wild guess for free convection (typically 2-25)
        self.l1 = 28E-3 # Sidelength block
        self.t1 = 3E-3 # Thickness of Al layer
        self.t2 = self.t1
        self.l2 = self.l1-2*self.t1
        
        self.A1 = (2*(self.l1-2*self.t2) + 2*self.l2) * 1 # [m2] per meter length
        self.A2 = (2*self.l1 + 2*self.l2) * 1 # [m2] per meter length
        
        self.emiss = 0.18 # Emissivity roughly polished aluminium -> https://www.engineeringtoolbox.com/radiation-heat-emissivity-aluminum-d_433.html
        self.sigma = 5.67E-8 # [W/m2K4] Stefan-Boltzmann constant
    
    def q_current(self, I, R):
        return I**2 * R * self.N
    
    def dT_cond(self, I, R):
        dT_cond = - self.t1 * (self.q_current(I, R))/(self.k*self.A1)
        return dT_cond
    
    def T2(self, I, R, Tamb=20):
        dT_cond = self.dT_cond(I, R)
        poly4 = self.emiss*self.sigma/self.k
        poly1 = self.h/self.k
        poly0 = - (self.h/self.k*Tamb+self.emiss*self.sigma/self.k*Tamb - dT_cond*self.A1/self.A2)
        sols = np.roots([poly4, 0, 0, poly1, poly0])
        # Return only positive, non-imaginary roots:
        return np.real([i for i in list(sols) if np.imag(i)==0 and np.real(i)>=0][0])
    
    def T1(self, I, R, Tamb=20):
        return self.T2(I, R, Tamb) - self.dT_cond(I, R)
    
    def listmethods(self):
        """Function that lists all non-reserved methods in the class"""
        object_methods = [method_name for method_name in dir(self) if callable(getattr(self, method_name))]        
        for method in object_methods:
            if method[0:2] != "__":
                print(str(method+"()"))