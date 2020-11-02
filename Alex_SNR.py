# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 14:13:12 2020

@author: AlexO'S
"""

import numpy as np
import matplotlib.pyplot as plt


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
#enter radius as 4cm and length as 10cm
#results in sample volume of ~ 0.5L
rs = (float(input('Enter sample radius (cm) - ')))/100
ls = (float(input('Enter sample length (cm) - ')))/100 
N = (2*(np.pi)*(rs**2)*ls*rhow)/mw 


# Coil parameters
Dc = (float(input('Enter coil diameter (cm) - ')))/100
Hc = (float(input('Enter coil height (cm) - ')))/100

# list of numbers with given interval and range 
def createList(r1, r2): 
    return np.arange(r1, r2+0.01, 0.01)     
#Functon generates list of diameter of wire in mm, 
#range: 0.1-0.5mm,interval: 0.01mm 
r1, r2 = 0.1, 0.5
Dw = createList(r1,r2) 

#Hence list of number of turns per layer
Nc = np.reciprocal(Dw)*1000*Hc 


# Other parameters
Bp = (float(input('Enter polarising field strength (mT) - ')))/1000


# Calculations
#constant
n1 = perm*w0*(np.pi)*N*(g**2)*(hbar**2)*Hc
d1 = 8*(np.sqrt(rhoe*(KB**3)*(T**3)*b*Dc))
C = n1/d1

#list is a protected keyword by python, cannot use it
#without causing errors
#This will hold multiple SNRc and V lists
collective_SNRcs = []
collective_Vs = []

#for loop for each desirable value  
#of the number of layers Nl from 1 to 10
#11 numbers i.e. 0-10, but don't count zero!
for Nl in range(1,11):
    n2 = (np.sqrt(Nl))*Bp
    d2 = np.sqrt(Nc)
    p = n2/d2 #parameter
    SNR = C*p
    SNRc = SNR/(np.sqrt((2*np.pi)*((Hc**2)+(Dc**2)))) #corrected SNR
    #append the vector of SNRc values for a given nl 
    #to the list of lists 
    collective_SNRcs.append(SNRc)
    
    #repeat process for EMF
    n3 = perm*w0*N*(g**2)*(hbar**2)*Nc*Nl*Bp
    d3 = 4*KB*T*np.sqrt(((Hc**2)+(Dc**2)))
    V = (n3/d3)*(10**(6)) #Voltage induced in microvolts
    collective_Vs.append(V) 
#having constructed the lists of SNRc EMF values 
#for each given Nl, plot them
  
    
# Plots against Dw
#first SNR    
fig, ax1 = plt.subplots(1,1)
ax1.set_ylabel('Finite multilayer SNR',fontsize='x-large')
ax1.set_xlabel('Diameter of wire [mm]',fontsize='x-large')
#call plot for each one list of SNRc values in the list of lists.
for i in range(np.size(collective_SNRcs,axis=0)):
    ax1.plot(Dw, collective_SNRcs[i])
#label each line with its respective Nl value    
ax1.legend(('Nl = 1', 'Nl = 2', 'Nl = 3', 'Nl = 4', 'Nl = 5', 'Nl = 6', 'Nl = 7', 'Nl = 8', 'Nl = 9', 'Nl = 10'))    
 
#repeat for EMF   
fig, ax2 = plt.subplots(1,1)
ax2.set_ylabel('EMF [μV]',fontsize='x-large')
ax2.set_xlabel('Diameter of wire [mm]',fontsize='x-large')
for i in range(np.size(collective_Vs,axis=0)):
    ax2.plot(Dw, collective_Vs[i])
ax2.legend(('Nl = 1', 'Nl = 2', 'Nl = 3', 'Nl = 4', 'Nl = 5', 'Nl = 6', 'Nl = 7', 'Nl = 8', 'Nl = 9', 'Nl = 10'))    
plt.show() 