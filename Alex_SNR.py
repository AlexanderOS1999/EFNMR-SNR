# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 14:13:12 2020

@author: AlexO'S
"""

import numpy as np

# Set constants

perm = (4*np.pi)*(10**(-7)) # permeability
w0 = (2*np.pi)*2000         # larmor frequency in rad/s
g = 2.675E8                 # gyromagnetic ratio in rad/s /Tesla
hbar = 1.0546E-34           # hbar
rhow = 997                  # density of water in kg/m^3
mw = 3E-26                  # water molecule mass in kg
KB = 1.381E-23              # Boltzmann constant in J/K
T = 288                     # temperature in K
b = (2*np.pi)*2000          # bandwidth in rad/s        
rhoe = 1.68E-8              # copper resistivity in Ohm m 

# Sample dimensions
rs = (float(input('Active sample dimensions: Enter radius (cm) - ')))/100
ls = (float(input('Enter length (cm) - ')))/100 
N = (2*(np.pi)*(rs**2)*ls*rhow)/mw 
#enter radius as 4cm and length as 10cm
#results in sample volume of ~ 0.5L

# Coil parameters
Hc = (float(input('Enter coil height (cm) - ')))/100
Dc = (float(input('Enter coil diameter (cm) - ')))/100
Nl = int(input('Enter number of layers - '))

# list of numbers with given interval and range 
def createList(r1, r2): 
    return np.arange(r1, r2+0.01, 0.01) 
      
# Driver Code 
r1, r2 = 0.1, 0.5

Dw = createList(r1,r2) #List of diameter of wire in mm, 
                       #range: 0.1-0.5mm,interval: 0.01mm

Nc = np.reciprocal(Dw)*1000*Hc #Hence list of number of turns per layer

# Other parameters
Bp = (float(input('Enter polarising field strength (mT) - ')))/1000

# Calculations
n1 = perm*w0*(np.pi)*N*(g**2)*(hbar**2)*Hc
d1 = 8*(np.sqrt(rhoe*(KB**3)*(T**3)*b*Dc))
C = n1/d1 #constant
n2 = (np.sqrt(Nl))*Bp
d2 = np.sqrt(Nc)
p = n2/d2 #parameter
SNR = C*p
SNRc = SNR/(np.sqrt((2*np.pi)*((Hc**2)+(Dc**2)))) #List of corrected
                                                  #SNR values to be plotted
print(SNRc)                                                  
n3 = perm*w0*N*(g**2)*(hbar**2)*Nc*Nl*Bp
d3 = 4*KB*T*np.sqrt(((Hc**2)+(Dc**2)))
V = (n3/d3)*(10**(6)) #Voltage induced in microvolts

# Plots against Dw or Nc
import matplotlib.pyplot as plt

plt.figure()
plt.subplot(211)
plt.plot(Dw, SNRc)
plt.xlabel('Diameter of wire [mm]')
plt.ylabel('Finite multilayer SNR')

plt.subplot(212)
plt.plot(Dw, V)
plt.xlabel('Diameter of wire [mm]')
plt.ylabel('EMF [μV]')
plt.show()