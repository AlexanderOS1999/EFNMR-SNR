"""
Created on Sun Nov  8 11:10:28 2020
@authors: Luis Mestre & Ewan Kilpatrick

"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.optimize
from scipy.fft import fft, ifft

SAMPLE_PERIOD = 0.00001


#function for signal exponential
def func_sig(t2, A2, K2, C2):
    return A2 * np.exp(-1*t2/K2) + C2
#linear fit function for actual signal
def fit_exp_linear_sig(t, y3, C4=1):
    y3 = y3 - C4
    y3 = np.log(y3)
    K2, A_log = np.polyfit(t, y3, 1)
    A2 = np.exp(A_log)
    return A2, K2


#plotting the graphs on the same window if you want, with the legend of the fits as well
def plot(ax, t, y, noisy_y, fit_y, fit_parms):
    #A0, K0, C0 = orig_parms
    A, K, C = fit_parms

    #ax.plot(t, y, 'k--', 
      #label='Actual Function:\n $y = %0.2f e^{%0.2f t} + %0.2f$' % (A0, K0, C0))
    ax.plot(t, fit_y, 'b-',
      label='Fitted Function:\n $y = %0.2f e^{%0.2f t} + %0.2f$' % (A, K, C))
    ax.plot(t, noisy_y, 'ro')
    ax.legend(fancybox=True, shadow=False)
    

#generating data
def createList(r1, r2): 
    #SAMPLE_PERIOD (Time between samples)
    return np.arange(r1, r2+SAMPLE_PERIOD, SAMPLE_PERIOD)     


# Set constants
perm = (4*np.pi)*(10**(-7)) # permeability
g = 2.675E8                 # gyromagnetic ratio in rad/s /Tesla
hbar = 1.0546E-34           # hbar
rhow = 997                  # density of water in kg/m^3
mw = 3E-26                  # water molecule mass in kg
KB = 1.381E-23              # Boltzmann constant in J/K
T = 298                     # temperature in K

#parameters for curve
V = 332                 # sample volume (ml)
Bp = 5                  # prepolarising field strength (mT)
T1 = 0.35               # sample t1 (s)
W_0 = 2.086             # Larmor frequency (krad per second)
SNR =  100              # SNR ratio

#calculating N           
N = (2*V*1E-6*rhow)/mw
#calculating max magnetisation Curie's law
M_0 = ((perm*hbar**2*g**2*(N))*Bp*1E-3/(4*KB*T))*1E16

# Generate signal points from a pulse sequence T1, T2 based on the paramaters 
tmin, tmax = 0, T1+1
num = 30
t = np.linspace(tmin, tmax, num)
y3 = func_sig(t, M_0, T1, M_0)

# Add noise to signal exponential fit
noisy_y2 = y3 + 5*V*Bp*1E-3* (np.random.random(num) - 0.5)
  

# FID plot data
Dw2 = createList(0, T1+1) 
Data = M_0 * np.sin(W_0 * 1000 * 2 * np.pi * Dw2) * np.exp(-(Dw2)/T1)
size = len(Data)


#powerline noise model

"""function to generate a the noise (when summed over, it will generate a convolution of
a series of lorenztians convoluted with a gaussian to simulate the powerline 50Hz 
harmonics picked up with a tuned LC circuit)"""


def fct_powerline(n):
    A3 = 1 #amplitude at 2.086kHz
    omega0 = 2086
    omega50 = 50
    sigma_envelope = 200
    phase = np.random.random()
    harmonic = n * omega50
    return A3 * (1/np.sqrt(2*np.pi*sigma_envelope**2))*np.sin(2*np.pi*harmonic*Dw2 + phase) * np.exp(-(harmonic-omega0)**2 /(sigma_envelope)**2) 

noise = 0
noise_peaks = 80

for n in range(0 , noise_peaks):
    noise  = noise + fct_powerline(n)


"""
generating data with a normal distribution noise (not uniform distribution like np.random.random)
this is the closed analogue to jonhson noise while random.random is closer to white noise
note increase the interval inside=more random noise with varying ampl. while multiplying
the overall noise function just increases the overall noise ampl. not the "noisiness"""
noise2 = (np.random.normal(0, 10, len(Dw2))-0.5)

# Combining noises with data
Data_noise = Data + 4E3*noise + noise2


# Linear Fit plot for signal
fig, ax3 = plt.subplots(1,1)
A2, K2 = fit_exp_linear_sig(t, y3, M_0)
fit_y2 = func_sig(t, M_0, T1, M_0)
plot(ax3, t, y3, noisy_y2, fit_y2, (A2, K2, A2))
ax3.set_title('Simulated pulse sequence results fit')
ax3.set_ylabel('Signal amplitude (a.u)' ,fontsize='x-large')
ax3.set_xlabel('Time [s]',fontsize='x-large')
ax3.grid(True, which='both')
plt.show()

#FID plot 
fig, ax4 = plt.subplots(1,1)
ax4.plot(Dw2 , Data_noise , 'k')
ax4.set_title('Simulated FID with powerline and thermal noise')
ax4.axhline(y=0, color='black')
ax4.set_ylabel('signal amplitude (a.u)',fontsize='x-large')
ax4.set_xlabel('Time (s)',fontsize='x-large')
ax4.grid(True, which='both')
plt.show()

# Convert FID data to look like Arduino output
x = np.linspace(1,len(Dw2),len(Dw2))
Data_noise_arduino = []
Data_noise_shift = Data_noise - np.amin(Data_noise)
m = np.amax(Data_noise_shift)
for i in range(len(Data_noise_shift)):
    Data_noise_arduino.append((1024/m)*(Data_noise_shift[i]))

# FID data print to file
writeFile = open('FIDdata.csv','w')
for i in range(len(Dw2)):
    writeFile.write("%d , %d\n" % (x[i],np.round(Data_noise_arduino[i])))    
        
    
