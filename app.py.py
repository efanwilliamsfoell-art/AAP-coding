# while loop (total number of muons)
import numpy as np

import streamlit as stl

N0 = stl.number_input("Initial Number of Muons (N0):", value=563)
initial_v_ratio = stl.number_input("Initial Speed (as fraction of c):", value=0.9952, format="%.4f")
altitude_m = stl.number_input("Travel Distance / Altitude (meters):", value=1917.0)

if stl.button("Calculate Surviving Muons"):
  N0 = 563
  alt= 1907 #altitude at Mt. Washington
  mrest = 1.88353*10**(-28) #rest mass of a muon
  c = 2.99792*10**8 #speed of light in a vacuum
  v0 = 0.9952*c
  dx = 0.001
  Lorentz_factor0 = 1/(np.sqrt(1-(v0)**2/c**2))
  mulifetime = 2.19698*10**(-6)
  E0 = Lorentz_factor0*mrest*c**2
  pi = np.pi
  me = 9.10938*10**(-31) # electron mass
  e = 1.60218*10**(-19) #electron charge
  e0 = 8.85419*10**(-12) #vacuum pemittivity
  Na = 6.022141*10**(23) #avogadro's number
  Z = 3.63591 #average atomic number
  z = -1 #charge of muon in multiples of electron charge
  I = 10*Z*e #Mean excitation energy, I, for air in the atmosphere
  p = 1.225 #
  A = 28.96
  Mu = 1*10**(-3)
  x = 0

  n0 = [N0] #list for number of muons left
  X = [x]

  while x < alt:

    n = (Na*Z*p)/(A*Mu)
    B = v0/c
  

    BB1 = 4*pi/(me*c**2)
    BB2 = n*z**2/(B)**2
    BB3 = (e**2/(4*pi*e0))**2
    BB4 = np.log((2*me*c**2*B**2)/(I*(1-(v0/c)**2))) - (v0/c)**2  

    de = BB1*BB2*BB3*BB4*dx
    E = E0 - de

    v1 = c*np.sqrt(1 - (mrest*c**2/E)**2)

    Lorentz_factor1 = E/(mrest*c**2)
    Lorentz_factormean = (Lorentz_factor0 + Lorentz_factor1)/2
    vmean= (v0 + v1)/2
    dt = dx/vmean
    x += dx

    Nt = N0*np.exp(-dt/(Lorentz_factormean*mulifetime))

    E0 = E
    v0 = v1
    Lorentz_factor0 = Lorentz_factor1
    N0 = Nt

    X.append(x)
    n0.append(N0)

stl.success(f"Muons remaining at target distance: 'n0[-1]')

