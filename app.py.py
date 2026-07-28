# while loop (total number of muons)
import numpy as np

import streamlit as stl

import matplotlib as plt

# Create a narrow column (ratio 1) and an empty wide space (ratio 3)
col1, _ = stl.columns([1, 3])

with col1: # stremlit layout
    N0_input = stl.number_input("Initial Number of Muons (N0):", value=563)
    initial_v_ratio = stl.number_input("Initial Speed (as fraction of c):", value=0.9952, format="%.4f")
    altitude_m = stl.number_input("Travel Distance / Altitude (metres):", value=1917.0)

    if stl.button("Calculate Surviving Muons"): #if that button is pressed
        N0 = N0_input 
        alt= altitude_m 
        c = 2.99792*10**8 #speed of light in a vacuum (because it has to be before v0
        v0 = initial_v_ratio * c #variables on website

        
        mrest = 1.88353*10**(-28) #rest mass of a muon
        dx = 0.01 # each step of distance
        Lorentz_factor0 = 1/(np.sqrt(1-(v0)**2/c**2))
        mulifetime = 2.19698*10**(-6) #average muon lifetime
        E0 = Lorentz_factor0*mrest*c**2 # original total energy (relativistic)
        pi = np.pi #the number pi
        me = 9.10938*10**(-31) # electron mass
        e = 1.60218*10**(-19) #electron charge
        e0 = 8.85419*10**(-12) #vacuum pemittivity
        Na = 6.022141*10**(23) #avogadro's number
        Z = 7.23 #average atomic number of dry air
        z = -1 #charge of muon in multiples of electron charge
        I = 10*Z*e #Mean excitation energy, I, for air in the atmosphere
        p = 1.225 #air density
        A = 28.96 #average atomic mass of the air
        Mu = 1*10**(-3)#molar mass constant
        
        x = 0 # intial value of x
        mux = 0
        
        MUX = [mux] #list for the distance the Earth moves toward the muon
        n0 = [N0] #list for number of muons left
        X = [x]#list for distance
      
        while x < alt:
      
          n = (Na*Z*p)/(A*Mu) #
          B = v0/c

          #Bethe bloch energy change
          BB1 = 4*pi/(me*c**2)
          BB2 = n*z**2/(B)**2
          BB3 = (e**2/(4*pi*e0))**2
          BB4 = np.log((2*me*c**2*B**2)/(I*(1-(v0/c)**2))) - (v0/c)**2  
          de = BB1*BB2*BB3*BB4*dx
          E = E0 - de
      
    #Lorentz factor calculations
          Lorentz_factor1 = E/(mrest*c**2)
          Lorentz_factormean = (Lorentz_factor0 + Lorentz_factor1)/2

          #Velocity calculations  
          v1 = c*np.sqrt(1 - (mrest*c**2/E)**2)        
          vmean= (v0 + v1)/2

          #updating distance and time  
          dt = dx/vmean
          x += dx
          mu_step = dx/Lorentz_factormean
          mux += mu_step

          #decay of muons in given time  
          Nt = N0*np.exp(-dt/(Lorentz_factormean*mulifetime))

          #resetting cycle  
          E0 = E
          v0 = v1
          Lorentz_factor0 = Lorentz_factor1
          N0 = Nt

          #recording values for graphing  
          X.append(x)
          n0.append(N0)
          MUX.append(mux)

        stl.success(f"Muons remaining at target distance: {n0[-1]:.1f}")
        fig, ax = plt.subplots(figsize=(8, 4))

# Plotting Earth-frame distance (X) vs. Surviving Muons (n0)
        ax.plot(X, n0, color="blue", linewidth=2, label="Surviving Muons")

        ax.set_xlabel("Earth Frame Distance (m)")
        ax.set_ylabel("Number of Muons")
        ax.set_title("Muon Decay Over Traveled Distance")
        ax.grid(True, linestyle="--", alpha=0.6)
        ax.legend()

# Display the figure in Streamlit
        stl.pyplot(fig)
