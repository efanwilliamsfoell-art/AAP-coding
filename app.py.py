# while loop (total number of muons)
import numpy as np

import streamlit as stl

import matplotlib.pyplot as plt

stl.set_page_config(page_title="Cosmic ray muon decay",layout="wide")  # <-- This spreads the app to fill the entire screen width

# 2. Create 2 side-by-side columns 
# [1, 2] means col1 takes 1/3 of screen width, col2 takes 2/3
col1, col2 = stl.columns([1, 2])

# 3. YOUR SNIPPET (Left Side: Inputs & Controls)
with col1:
    stl.subheader("Model Inputs")
    N0_input = stl.number_input("Initial Number of Muons (N0):", value=563)
    initial_v_ratio = stl.number_input("Initial Speed (as fraction of c):", value=0.9952, format="%.4f")
    altitude_m = stl.number_input("Travel Distance / Altitude (metres):", value=1907.0)

# 4. Right Side: Results & Graphs (Stretches across remaining space)
with col2:
    stl.subheader("Results of the model")
    

    if stl.button("Calculate Surviving Muons"): #if that button is pressed
        N0 = N0_input 
        N0clas = N0_input
        alt= altitude_m 
        c = 2.99792*10**8 #speed of light in a vacuum (because it has to be before v0
        v0 = initial_v_ratio * c #variables on website

        
        mrest = 105.65837 #rest mass of a muon in MeV
        dx = 0.01 # each step of distance
        gamma0 = 1/(np.sqrt(1-(v0)**2/c**2))
        mulifetime = 2.19698*10**(-6) #average muon lifetime
        E0 = gamma0*mrest # original total energy (relativistic)
        pi = np.pi #the number pi
        me = 0.51099895 # electron mass in MeV 
        e = 1.60218*10**(-19) #electron charge
        e0 = 8.85419*10**(-12) #vacuum pemittivity
        Na = 6.022141*10**(23) #avogadro's number
        Z = 7.23 #average atomic number of dry air
        z = -1 #charge of muon in multiples of electron charge
        I = 85.7*10**(-6) #Mean excitation energy, I, for air in the atmosphere
        p = 0.001225 #air density at sea level in g/cm^3
        A = 28.966 #average atomic mass of the air
        Mu = 1*10**(-3)#molar mass constant
        K = 0.307075 # units are MeV /mol cm^2

        p0 = 1.225 # air density at sea level in kg/m^3
        L = 0.0065 #temperature lapse rate K/m
        T0 = 288.15 # temperature at sea level
        R = 8.31446 # ideal gas constant J/(mol·K)
        g = 9.80665 # gravitational acceleration
        m = A/1000 # average atomic mass of dry air in kg/mol
        
        x = 0 # intial value of x
        mux = 0
        
        MUX = [mux] #list for the distance the Earth moves toward the muon
        n0 = [N0] #list for number of muons left (relativity)
        n0clas = [N0clas]
        X = [x]#list for distance
      
        while x < alt:
          M = gamma0*mrest #relativistic mass
          h = alt - x
  
          #figuring out maximum energy lost in a single collision
          Bsquared = 1 - 1/(gamma0**2)
          r = me/mrest # ratio of rest mass to electron mass
          numerator = 2*me*Bsquared*gamma0**2
          denominator = 1 + 2*gamma0*r + r**2
          Wmax = numerator/denominator
          SP1 = K*z**2*Z/(A*Bsquared)
          SP2 = 0.5*np.log(2*me*Bsquared*gamma0**2*Wmax/I**2) - Bsquared

          p = p0*(1- (L*h)/T0)**(g*m/(R*L) - 1) #pressure at each altitude

          de = SP1*SP2*(dx*100)*p/1000
          E = E0 - de
      
    #Lorentz factor calculations
          gamma1 = E/(mrest)
          gammamean = (gamma0 + gamma1)/2

          #Velocity calculations  
          v1 = c*np.sqrt(1 - (mrest/E)**2)        
          vmean= (v0 + v1)/2

          #updating distance and time  
          dt = dx/vmean
          x += dx
          mu_step = dx/gammamean
          mux += mu_step

          #decay of muons in given time  
          Ntrel = N0*np.exp(-dt/(gammamean*mulifetime))
          Nclas = N0clas*np.exp(-dt/mulifetime)

          #resetting cycle  
          E0 = E
          v0 = v1
          gamma0 = gamma1
          N0 = Ntrel
          N0clas = Nclas

          #recording values for graphing  
          X.append(x)
          n0.append(N0)
          MUX.append(mux)
          n0clas.append(N0clas)

        stl.success(f"Muons remaining at target distance: {n0[-1]:.0f}")
        stl.success(f"Distance the muons feel earth travelled towards them: {MUX[-1]:.1f} m")
        fig, ax = plt.subplots(figsize=(8, 4))

        exp_x = [0, alt]  # Mt. Washington (0m) and Sea Level (1907m)
        exp_y = [563, 408]  # Measured counts
        exp_yerr = [10, 9]  # Reported uncertainties

        ax.errorbar(exp_x, exp_y, yerr=exp_yerr, fmt="ko", capsize=5, capthick=1.0, markersize=4,label="Real data",)

# Plotting Earth-frame distance (X) vs. Surviving Muons (n0)
        ax.plot(X, n0, color="blue", linewidth=2, label="Relativistic model")
        ax.plot(X, n0clas, color ="red", linewidth=2, label="Classical model")

        ax.set_xlabel("Earth frame distance (m)")
        ax.set_ylabel("Number of muons")
        ax.set_title("Muon decay over Earth-frame distance")
        ax.grid(True, linestyle="--", alpha=0.6)
        ax.legend()

# Display the figure in Streamlit
        stl.pyplot(fig)

        
