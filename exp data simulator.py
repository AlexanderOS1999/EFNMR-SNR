# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 11:10:28 2020

@author: Luis Mestre 
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.optimize
from scipy.fft import fft, ifft
import csv

#Effectively 1/(your sampling frequency) Used in creating Time axis and Frequency axis!
SAMPLE_PERIOD = 0.00001

#function for magnetisation
def func_mag(t, A, K, C):
    return -1*A * np.exp(t/K) + C
#function for signal exponential
def func_sig(t2, A2, K2, C2):
    return A2 * np.exp(-1*t2/K2) + C2
#linear fit function for magnetisation
def fit_exp_linear_mag(t, y, C3=1):
    y = y - C3
    y = np.log(-1*y)
    K, A_log = np.polyfit(t, y, 1)
    A = np.exp(A_log)
    return A, K
#linear fit function for actual signal
def fit_exp_linear_sig(t, y3, C4=1):
    y3 = y3 - C4
    y3 = np.log(y3)
    K2, A_log = np.polyfit(t, y3, 1)
    A2 = np.exp(A_log)
    return A2, K2

#generating data
def createList(r1, r2): return np.arange(r1, r2+0.00001, 0.00001) 
    
#FID eq
def FID_func(t, A, K, W_0):
    y = np.sin (W_0 * 1000 * t)*A * np.exp(-t/K)
    return y
#non-linear fit for exponential
def fit_exp_nonlinear(t, y):
    opt_parms, parm_cov = sp.optimize.curve_fit(func_mag, t, y, maxfev=1000)
    A, K, C = opt_parms
    return A, K, C
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
    
#just to set the tick frequency to make the grid better etc    
def Ticks(mi, ma, I):
    ticks = np.arange(mi, ma , I)
    return ticks

#generating data
def createList(r1, r2): 
    #SAMPLE_PERIOD (Time between samples)
    return np.arange(r1, r2+SAMPLE_PERIOD, SAMPLE_PERIOD)     



def main():
    
    # Set constants
    perm = (4*np.pi)*(10**(-7)) # permeability
    w0 = (2*np.pi)*2000         # larmor frequency in rad/s
    g = 2.675E8                 # gyromagnetic ratio in rad/s /Tesla
    hbar = 1.0546E-34           # hbar
    rhow = 997                  # density of water in kg/m^3
    mw = 3E-26                  # water molecule mass in kg
    KB = 1.381E-23              # Boltzmann constant in J/K
    T = 298                     # temperature in K
    b = (2*np.pi)*2000          # bandwidth in rad/s        
    rhoe = 1.68E-8              # copper resistivity in Ohm m 
    
    #parameters for curve
    V = 500 #float(input('enter volume (ml) - ' ))
    Bp = 10 #float(input('enter the propolarising field (mT) - '))
    T1 = 0.35 #float(input('enter sample T1 (s) - '))
    W_0 = 2*np.pi*2#*float(input('enter larmor frequency (KHz) - '))
    SNR =  100# float(input('enter your expected SNR - '))
    
    #calculating N           
    N = (2*V*1E-6*rhow)/mw
    #print (N)
    #caclulating max magnetisation curie's law
    M_0 = ((perm*hbar**2*g**2*N)*Bp*1E-3/(4*KB*T))*1E16
    #print (M_0)

    # Generate magnetisation curve data T1, T2 based on the paramaters 
    tmin, tmax = 0, T1+1
    num = 30
    t2 = np.linspace(tmin, tmax, num)
    y = func_mag(t2, M_0, -1*T1, M_0)

    # Add noise to exponential fit
    noisy_y = y + 5*V*Bp*1E-3* (np.random.random(num) - 0.5)
    
    #Generate signal points from a pulse sequence T1, T2 based on the paramaters 
    y3 = func_sig(t2, M_0, T1, M_0)

    # Add noise to exponential fit
    noisy_y2 = y3 + 5*V*Bp*1E-3* (np.random.random(num) - 0.5)

    #fig, ax1 = plt.subplots(1,1)
    #ax1 = fig.add_plot(2,1,1)
    
    #ax2 = fig.add_plot(2,1,2)

    # Non-linear Fit don't really need it now but might come in handy later
    #fig, ax1 = plt.subplots(1,1)
    #A, K, C = fit_exp_nonlinear(t, noisy_y)
    #fit_y = model_func(t, A, K, C)
    #plot(ax1, t, y, noisy_y, fit_y, (A0, K0, C0), (A, K, C0))
    #ax1.set_title('Non-linear Fit')
    
    #FID plot data
    Dw2 = createList(0, T1+1) 
    Data = FID_func(Dw2, M_0, T1, W_0)

    
    #Function generates list of data
    Data = FID_func(Dw2, M_0, T1, W_0)
    size = len(Data)

    
    #powerline noise model
    
    #function to generate a the noise (when summed over, it will generate a convolution of a series of 
    #lorenztians convoluted with a gaussian to simulate the powerline 50Hz harmonics picked up with
    # a tuned LC circuit)
    
    A2 = 1 #amplitude at 2kHz
    omega0 = 2000
    omega50 = 50
    sigma_envelope = 200
    
    #function to generate a the noise (when summed over, it will generate a convolution of a series of 
    #lorenztians convoluted with a gaussian to simulate the powerline 50Hz harmonics picked up with
    # a tuned LC circuit)
    def fct_powerline(n):
        global Dw2,sigma_envelope,omega0,A3
        Dw2 = createList(0, T1+1) 
        #time = np.arange(0,T1+4,0.00001)
        A3 = 1 #amplitude at 2kHz
        omega0 = 2000
        omega50 = 50
        sigma_envelope = 200
        phase = np.random.random()
        omega_n = n*omega50
        return A3 * (1/np.sqrt(2*np.pi*sigma_envelope**2))*np.exp(-(omega_n-omega0)**2 /(sigma_envelope)**2) * np.sin(2*np.pi*omega_n*Dw2 + phase)
    
    #time = np.arange(0,7,0.00001)
    noise = 0
    noise_peaks = 80
    
    for n in range(0,noise_peaks):
        noise  = noise + fct_powerline(n)
    
    SAMPLE_PERIOD = 0.00001
    #FFT of noise
    noise_FFT = fft(noise) / len(noise) 
    
    #c center of noise level (average of distribution) and max_noise is the maximum noise amplitude possible
    c = 0
    max_noise = M_0/SNR
    #c center of noise level (average of distribution) and max_noise is the maximum noise amplitude possible
    #########################################
    #generating data with a normal distribution noise (not uniform distribution like np.random.random)
    #this is the closed analogue to jonhson noise while random.random is closer to white noise
    #note increase the interval inside=more random noise with varying ampl. while multiplying
    #the overall noise function just increases the overall noise ampl. not the "noisyness"
    noise2 = (np.random.normal(c, 100, size)-0.5)
    #test noise
    #noise_r = 17*(np.random.random(size) - 0.5)
    
    #total data with both normally distributed noise and powerline noise
    x = np.size(noise_FFT)
    f_noise = np.fft.fftfreq(len(noise_FFT), SAMPLE_PERIOD)[:x//2] 
    y2_noise = np.abs(noise_FFT[:x//2])
    
    Data_noise = Data + 4E3*noise + (np.random.normal(c, 10, size)-0.5)
 
    #computing FFT of data
    #FFT is a cunt when it comes to scaling! Divide fft by the length of data to preserve magnitude
    Data_FFT = fft(Data_noise) / len(Data)
    #print (Data_FFT)
    x = np.size(Data) 
    #Way to generate the frequencies based on SAMPLE_PERIOD (also since you dont like negatives, throw them away)
    f = np.fft.fftfreq(len(Data), SAMPLE_PERIOD)[:x//2] 
    #generating y axis 
    y2 = np.abs(Data_FFT[0:x//2])
    #When you only take the positive frequencies, the magnitudes actually double (since half the power is technically in the negative freq...)
    #But, the DC (index 0) is a bit special and is not effected by this so only multiply the frequecies above 0 by 2
    y2[1:] *= 2
    #print (y.shape)
    #print(t.shape)
  
    # Linear Fit plot for magnetization
    fig, ax2 = plt.subplots(1,1)
    A, K = fit_exp_linear_mag(t2, y, M_0)
    fit_y = func_mag(t2, M_0, -1*T1, M_0)
    A2=-A
    plot(ax2, t2, y, noisy_y, fit_y, (A2, K, A))
    ax2.set_title('Simulated pulse sequence results fit for magnetization')
    ax2.set_ylabel('Magnetization signal (a.u)' ,fontsize='x-large')
    ax2.set_xlabel('Time [s]',fontsize='x-large')
    #controlling frequency of ticks
    #plt.xticks(Ticks(0, T1+6, 0.5))
    #limiting x axis
    #plt.xlim(0, T1+2.5)
    #adding grid
    ax2.grid(True, which='both')
    plt.show()
    
    # Linear Fit plot for signal
    fig, ax3 = plt.subplots(1,1)
    A2, K2 = fit_exp_linear_sig(t2, y3, M_0)
    fit_y2 = func_sig(t2, M_0, T1, M_0)
    plot(ax3, t2, y3, noisy_y2, fit_y2, (A2, K2, A2))
    ax3.set_title('Simulated pulse sequence results fit')
    ax3.set_ylabel('Signal amplitude (a.u)' ,fontsize='x-large')
    ax3.set_xlabel('Time [s]',fontsize='x-large')
    #controlling frequency of ticks
    #plt.xticks(Ticks(0, T1+6, 0.5))
    #limiting x axis
    #plt.xlim(0, T1+2.5)
    #adding grid
    ax3.grid(True, which='both')
    plt.show()
    
    #FID plot 
    fig, ax4 = plt.subplots(1,1)
    ax4.plot(Dw2, Data_noise)
    ax4.set_title('Simulated FID with noise')
    ax4.axhline(y=0, color='black')
    ax4.set_ylabel('signal amplitude (a.u)',fontsize='x-large')
    ax4.set_xlabel('Time (s)',fontsize='x-large')
    #controlling frequency of ticks
    #plt.xticks(Ticks(0, T1+5, 0.5))
    #limiting x axis
    #plt.xlim(0, T1+4)
    #adding grid
    ax4.grid(True, which='both')
    plt.show() 
    
    #FFT plot 
    fig, ax5 = plt.subplots(1,1)
    ax5.set_title('Simulated spectrum')
    ax5.set_ylabel('signal amplitude (a.u)',fontsize='x-large')
    ax5.set_xlabel('Frequency (Hz)',fontsize='x-large')
    #!!!!!!!!!!!!!!!!!!!!!!!!! Plot FFT using log frequency axis
    ax5.plot(f, y2)
    #set limit on x-axis for why just replace x with y
    plt.xlim(1600, 2400)
    #ax4.axhline(y=0, color='black')
    #plt.xticks(Ticks(1950, 2050, 10))
    ax5.grid(True, which='both')
    plt.show() 
    
if __name__ == '__main__':
    main()