# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 14:13:12 2020

@author: AlexO'S
"""

import numpy as np
import matplotlib.pyplot as plt


# Set constants
perm = (4*np.pi)*(10**(-7)) # permeability of free space
w0 = (2*np.pi)*2000         # larmor frequency in rad/s
g = 2.675E8                 # gyromagnetic ratio in rad/s /Tesla
hbar = 1.0546E-34           # hbar
rhow = 997                  # density of water in kg/m^3
mw = 3E-26                  # water molecule mass in kg
KB = 1.381E-23              # Boltzmann constant in J/K
T = 298                     # temperature in K
b = (2*np.pi)*2000          # bandwidth in rad/s        
rhoe = 1.68E-8              # copper resistivity in Ohm m 
f = 2000                    # operating frequency in Hz 
eta = 0.9                   # ratio of wire diameter to avg. distance 
                            # between two adjacent windings

# Sample dimensions
#enter diameter as 8cm and length as 10cm
#results in sample volume of ~ 0.5L
Ds = (float(input('Enter sample diameter (cm) - ')))/100
ls = (float(input('Enter sample length (cm) - ')))/100 
N = (2*(np.pi)*(Ds**2)*ls*rhow)/(4*mw) 


# Coil parameters
#note this is the diameter of a 1 layer coil
#convert to metres
Dc = (float(input('Enter coil inner diameter (cm) - ')))/100 
Hc = (float(input('Enter coil height (cm) - ')))/100

# list of numbers with given interval and range 
def createList(r1, r2): 
    return np.arange(r1, r2+0.005, 0.005)     
#Function generates list of diameter of wire in mm, 
#range: 0.1-1.15mm,interval: 0.005mm 
r1, r2 = 0.1, 1.15
Dw = createList(r1,r2) 
 
#Hence list of number of turns per layer
#convert Dw to m for correct calculation
Nc = Hc/(Dw*(10**-3)) 


# Other parameters
Bp = (float(input('Enter polarising field strength (mT) - ')))/1000
print()

# Calculations
#constants
n1 = perm*w0*(np.pi)*N*(g**2)*(hbar**2)*Hc
d1 = 8*(np.sqrt(rhoe*(KB**3)*(T**3)*b*Dc))
C = n1/d1
Ac = ((np.pi)*(Dc**2))/4 #area of coil in m^2
delta = np.sqrt((rhoe/(np.pi*f*perm*1))) #skin depth of copper in m
                                        #note function of f 
                                        #& relative perm approx.= 1
                                        
#list is a protected keyword by python, cannot use it
#without causing errors
#This will hold multiple e.g. SNRc lists
collective_SNRc = []
collective_V = []
collective_R = []
collective_L_inf = []
collective_L_fin = []
collective_Q = []
collective_Va = []
collective_SNRtot = []
collective_Amp = []
collective_Cap = []

#for loop for each desirable value  
#of the number of layers Nl from 1 to 14
#15 numbers i.e. 0-14, but don't count zero!
for Nl in range(1,15):
    n2 = (np.sqrt(Nl))*Bp
    d2 = np.sqrt(Nc)
    p = n2/d2 #parameter
    SNR = C*p
    SNRc = SNR/(np.sqrt((2*np.pi)*((Hc**2)+(Dc**2)))) #corrected SNR
    #append the vector of SNRc values for a given Nl 
    #to the list of lists 
    collective_SNRc.append(SNRc) #OVERESTIMATED?
    
    #repeat process for EMF
    n3 = perm*w0*N*(g**2)*(hbar**2)*Nc*Nl*Bp
    d3 = 4*KB*T*np.sqrt(((Hc**2)+(Dc**2)))
    V = (n3/d3)*(10**(6)) #Voltage induced in microvolts
    collective_V.append(V) #OVERESTIMATED?
    
    #repeat process for DC resistance
    R = (4*rhoe*Dc*(Nc**3)*Nl)/((Hc**2))
    collective_R.append(R)
    
    #repeat process for inductance of an infinte multilayer solenoid
    L_inf = (perm*(Nc**2)*(Nl**2)*Ac)/Hc
    collective_L_inf.append(L_inf)
    
    #repeat process for inductance of a finite multilayer solenoid
    #first calculate the average diameter of the coil
    #for given number of layers Nl
    #convert Dw to m in calculation!
    Dm = Dc+(Nl*Dw*(10**-3))
    L_fin = ((31.6*(10**-6))*((Dm/2)**2)*(Nc**2)*(Nl**2))/((3*Dm)+(9*Hc)+(10*Nl*(Dw*(10**-3))))
    collective_L_fin.append(L_fin)
    
    #repeat process for quality factor Q
    #take into account AC resistance
    #F is the ratio of total resistance to DC resistance
    #it is a function of parameter A
    A = ((((np.pi)/4)**(3/4))*(Dw*(10**-3))*np.sqrt(eta))/delta
    F = 1+((((5*(Nl**2))-1)*(A**4))/45) #CORRECT EQUATION?
    #from below plots, infinite treatment is an overestimation of L
    #therefore use finite treatment in calculation     
    Q = (w0*L_fin)/(R*F)
    collective_Q.append(Q)
    
    #repeat process for EMF multiplied
    #by resonant circuit
    #to be amplified by amplifier
    Va = Q*V #in microvolts
    collective_Va.append(Va)
    
    #repeat process for total SNR
    #replace constant bandwidth b 
    #with one that depends on Q
    #convert to dB
    SNRt = SNRc*np.sqrt(((b*Q)/w0))
    SNRtot = 20*np.log10(SNRt)
    collective_SNRtot.append(SNRtot)
    
    #repeat for amplication factor required for 1V in dB
    #convert Va to volts
    Amp = 20*np.log10((1/(Va*(10**-6))))
    collective_Amp.append(Amp)
    
    #repeat process for capacitance Cap
    Cap = L_fin/(((w0**2)*(L_fin**2))+((R*F)**2))
    collective_Cap.append(Cap)
#having constructed the list of lists for e.g. SNRtot values 
#for each given Nl, plot them
  
    
# Plots against Dw or CSA of the wire
CSA = ((np.pi)*(Dw**2))/4 #in mm^2
##uncomment this code below for thermal optimisation 
#of SNR & EMF only:
##first SNR    
#fig, ax1 = plt.subplots(1,1)
#ax1.set_ylabel('Finite multilayer SNR',fontsize='x-large')
#ax1.set_xlabel('Diameter of wire [mm]',fontsize='x-large')
##call plot for each one list of SNRc values in the list of lists.
#for i in range(np.size(collective_SNRcs,axis=0)):
#    ax1.plot(Dw, collective_SNRcs[i])
##label each line with its respective Nl value    
#ax1.legend(('Nl = 1', 'Nl = 2', 'Nl = 3', 'Nl = 4', 'Nl = 5', 'Nl = 6', 'Nl = 7', 'Nl = 8', 'Nl = 9', 'Nl = 10'))    
# 
##repeat for EMF   
#fig, ax2 = plt.subplots(1,1)
#ax2.set_ylabel('EMF [μV]',fontsize='x-large')
#ax2.set_xlabel('Diameter of wire [mm]',fontsize='x-large')
#for i in range(np.size(collective_Vs,axis=0)):
#    ax2.plot(Dw, collective_Vs[i])
#ax2.legend(('Nl = 1', 'Nl = 2', 'Nl = 3', 'Nl = 4', 'Nl = 5', 'Nl = 6', 'Nl = 7', 'Nl = 8', 'Nl = 9', 'Nl = 10'))    

#repeat for DC resistance   
fig, ax3 = plt.subplots(1,1)
ax3.set_ylabel('DC Resistance [Ω]',fontsize='x-large')
ax3.set_xlabel('Wire CSA [mm^2]',fontsize='x-large')
ax3.set_xlim(0.0, 1.0)
ax3.set_ylim(0, 40)
for i in range(np.size(collective_R,axis=0)):
    ax3.plot(CSA, collective_R[i])  
ax3.legend(('Nl = 1', 'Nl = 2', 'Nl = 3', 'Nl = 4', 'Nl = 5', 'Nl = 6', 'Nl = 7', 'Nl = 8', 'Nl = 9', 'Nl = 10', 'Nl = 11', 'Nl = 12', 'Nl = 13', 'Nl = 14')) 

##uncomment this code below for infinite treatment:
###repeat for inductance of infinite solenoid   
##fig, ax4 = plt.subplots(1,1)
##ax4.set_ylabel('Infinite multilayer inductance [H]',fontsize='x-large')
##ax4.set_xlabel('Diameter of wire [mm]',fontsize='x-large')
##for i in range(np.size(collective_L_infs,axis=0)):
##    ax4.plot(Dw, collective_L_infs[i])
##ax4.legend(('Nl = 1', 'Nl = 2', 'Nl = 3', 'Nl = 4', 'Nl = 5', 'Nl = 6', 'Nl = 7', 'Nl = 8', 'Nl = 9', 'Nl = 10'))

#repeat for inductance of finite solenoid   
fig, ax5 = plt.subplots(1,1)
ax5.set_ylabel('Finite multilayer inductance [mH]',fontsize='x-large')
ax5.set_xlabel('Wire CSA [mm^2]',fontsize='x-large')
ax5.set_xlim(0.0, 1.0)
ax5.set_ylim(0, 150)
for i in range(np.size(collective_L_fin,axis=0)):
    ax5.plot(CSA, ((10**3)*collective_L_fin[i]))
ax5.legend(('Nl = 1', 'Nl = 2', 'Nl = 3', 'Nl = 4', 'Nl = 5', 'Nl = 6', 'Nl = 7', 'Nl = 8', 'Nl = 9', 'Nl = 10', 'Nl = 11', 'Nl = 12', 'Nl = 13', 'Nl = 14'))

##uncomment this code below for comparison of inductance treatments:
##fig, ax6 = plt.subplots(1,1)
##ax6.set_ylabel('L [H]',fontsize='x-large')
##ax6.set_xlabel('Diameter of wire [mm]',fontsize='x-large')
##ax6.plot(Dw, collective_L_fins[0])
##ax6.plot(Dw, collective_L_infs[0])
##ax6.legend(('finite Nl = 1', 'infinite Nl = 1'))

#repeat for quality factor Q  
fig, ax7 = plt.subplots(1,1)
ax7.set_ylabel('Q',fontsize='x-large')
ax7.set_xlabel('Wire CSA [mm^2]',fontsize='x-large')
ax7.set_xlim(0.0, 1.0)
for i in range(np.size(collective_Q,axis=0)):
    ax7.plot(CSA, collective_Q[i])
ax7.legend(('Nl = 1', 'Nl = 2', 'Nl = 3', 'Nl = 4', 'Nl = 5', 'Nl = 6', 'Nl = 7', 'Nl = 8', 'Nl = 9', 'Nl = 10', 'Nl = 11', 'Nl = 12', 'Nl = 13', 'Nl = 14'))

#repeat for multiplied EMF Va  
fig, ax8 = plt.subplots(1,1)
ax8.set_ylabel('Multiplied EMF [μV]',fontsize='x-large')
ax8.set_xlabel('Wire CSA [mm^2]',fontsize='x-large')
ax8.set_xlim(0.0, 1.0)
for i in range(np.size(collective_Va,axis=0)):
    ax8.plot(CSA, collective_Va[i])
#ax8.legend(('Nl = 1', 'Nl = 2', 'Nl = 3', 'Nl = 4', 'Nl = 5', 'Nl = 6', 'Nl = 7', 'Nl = 8', 'Nl = 9', 'Nl = 10', 'Nl = 11', 'Nl = 12', 'Nl = 13', 'Nl = 14'))

#repeat for total SNR 
fig, ax9 = plt.subplots(1,1)
ax9.set_ylabel('Total SNR [dB]',fontsize='x-large')
ax9.set_xlabel('Wire CSA [mm^2]',fontsize='x-large')
ax9.set_xlim(0.0, 1.0)
for i in range(np.size(collective_SNRtot,axis=0)):
    ax9.plot(CSA, collective_SNRtot[i])
#ax9.legend(('Nl = 1', 'Nl = 2', 'Nl = 3', 'Nl = 4', 'Nl = 5', 'Nl = 6', 'Nl = 7', 'Nl = 8', 'Nl = 9', 'Nl = 10', 'Nl = 11', 'Nl = 12', 'Nl = 13', 'Nl = 14'))

#repeat for amplification factor required for 1V 
fig, ax10 = plt.subplots(1,1)
ax10.set_ylabel('Amplification factor [dB]',fontsize='x-large')
ax10.set_xlabel('Wire CSA [mm^2]',fontsize='x-large')
ax10.set_xlim(0.0, 1.0)
for i in range(np.size(collective_SNRtot,axis=0)):
    ax10.plot(CSA, collective_Amp[i])
ax10.legend(('Nl = 1', 'Nl = 2', 'Nl = 3', 'Nl = 4', 'Nl = 5', 'Nl = 6', 'Nl = 7', 'Nl = 8', 'Nl = 9', 'Nl = 10', 'Nl = 11', 'Nl = 12', 'Nl = 13', 'Nl = 14'))
plt.show()


# Optimisation
#index of max Q for user-inputted number of layers
optNl = int(input('Enter chosen optimised number of layers - '))
index_of_maxQ = np.where(collective_Q[(optNl-1)] == max(collective_Q[(optNl-1)]))
##uncomment this code below to verify that above function
##indeed finds the index of the maximum in the array:
#print('The index of max Q:', index_of_maxQ)
#print(collective_Q[(optNl-1)][index_of_maxQ])

#index of max total SNR
index_of_maxSNRtot = np.where(collective_SNRtot[(optNl-1)] == max(collective_SNRtot[(optNl-1)]))
#index of max multiplied EMF
index_of_maxVa = np.where(collective_Va[(optNl-1)] == max(collective_Va[(optNl-1)]))
#print maxima of curves corresponding to chosen number of layers
#and the corresponding CSA values at which they occur
print()
print("For Nl =",optNl,":")
print()
print("Maximum Q =",max(collective_Q[(optNl-1)]),"; at wire CSA =",CSA[index_of_maxQ],"mm^2")
print("Maximum total SNR =",max(collective_SNRtot[(optNl-1)]),"dB; at wire CSA =",CSA[index_of_maxSNRtot],"mm^2")
print("Maximum multiplied EMF =",max(collective_Va[(optNl-1)]),"μV; at wire CSA =",CSA[index_of_maxVa],"mm^2")
print()
# print("Note these will correspond to different values of wire diameter:")
# #find the CSA at which Q max occurs 
# print("Wire CSA at which maximum Q occurs =",CSA[index_of_maxQ],"mm^2")
# #find the CSA at which max total SNR occurs 
# print("Wire CSA at which maximum total SNR occurs =",CSA[index_of_maxSNRtot],"mm^2")
# #find the CSA at which max multiplied EMF occurs
# print("Wire CSA at which maximum multiplied EMF occurs =",CSA[index_of_maxVa],"mm^2")
# print()
print("Choose optimised CSA less than the CSA that gives rise to maximum Q,",CSA[index_of_maxQ],"mm^2, as the total SNR will still lie on a rough plateau around its maximum and the multiplied EMF will lie closer to its maximum. Also want to reduce Q somewhat as we do not want the circuit to be more resonant than the sample.")
print()
#Find SNRtot at CSA for which multiplied EMF is max 
print("Total SNR at which multiplied EMF is maximum =",collective_SNRtot[(optNl-1)][index_of_maxVa],"dB - too much of an appreciable drop from maximum.")
print("Therefore, choose optimised CSA by finding the point on the plots where there is a compromise between the drops from the respective maxima.")
print()


# Calculate theoretical predictions for quantities
#based on chosen wire size
# I estimate this to be about 0.24 mm^2 for 10 layers
# for our chosen dimensions earlier inputted
optCSA = float(input('Enter chosen optimised CSA (mm^2) - '))
print()
#find closest index in CSA array to this
#make new array containing difference between input and CSA values
𝛿_CSA = np.absolute(CSA-optCSA)
#index of optimised value for chosen number of layers
index_of_optCSA = np.where(𝛿_CSA == min(𝛿_CSA))

#coil specifications:
#find optimised diameter in Dw array using above index
print("Optimised coil specifications:")
print("Wire diameter =",Dw[index_of_optCSA],"mm")
#repeat for number of turns per layer
print("Number of turns per layer =",Nc[index_of_optCSA])

#calculate total length of wire needed for chosen number of layers
#Create empty array of diameters for given layer Dl
collective_Dl = []
for Nl in range(1,(optNl+1)):
    Dl = Dc+(2*(Nl-1)*(Dw[index_of_optCSA]*(10**-3)))
    collective_Dl.append(Dl)
#calculation accounts for variation in diameter for each layer    
lw = Nc[index_of_optCSA]*np.pi*sum(collective_Dl)     
print("Total length of wire required =",lw,"m")
print("Outer diameter of coil =",((10**(2))*collective_Dl[(optNl-1)]),"cm")    
print()

print("Calculated theoretical predictions for quantities:")
#repeat for DC resistance 
print("DC resistance =",collective_R[(optNl-1)][index_of_optCSA],"Ω")
#repeat for finite inductance 
print("Finite inductance =",((10**3)*collective_L_fin[(optNl-1)][index_of_optCSA]),"mH")
#calculate tuned capacitance based on above quantities 
print("Therefore, tuned capacitance =",((10**(9))*collective_Cap[(optNl-1)][index_of_optCSA]),"nF")
print()
#repeat for quality factor
print("Quality factor =",collective_Q[(optNl-1)][index_of_optCSA])
#repeat for total SNR
print("Total SNR =",collective_SNRtot[(optNl-1)][index_of_optCSA],"dB")
#repeat for multiplied EMF 
print("Multiplied EMF =",collective_Va[(optNl-1)][index_of_optCSA],"μV")
#calculate amp factor based on this EMF 
print("Amplification required for 1V =",collective_Amp[(optNl-1)][index_of_optCSA],"dB")