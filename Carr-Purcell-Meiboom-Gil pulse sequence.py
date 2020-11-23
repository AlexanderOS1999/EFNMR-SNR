# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 20:46:13 2020

@author: luis1
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.optimize
from scipy.fft import fft, ifft
import csv


def createList(r1, r2): return np.arange(r1, r2+0.00001, 0.00001) 

def main():
    
    T1 = 2
    T_2 = 2
    T_2_star = 0.15
    T_E = 1
    Amp = 1
    t = createList(0, T1+4) 
    
    def fct_pulse_sequence(n, T_2_star, T_2, T_E, Amp):
            global omega0
            #t = createList(0, T1+1) 
            #time = np.arange(0,T1+4,0.00001)
            #A3 = 1 #amplitude at 2kHz
            omega0 = 2000
            omega50 = T_E
            sigma_envelope = T_2_star
            phase = np.random.random()
            omega_n = n*omega50
    
            return Amp * ((1/np.sqrt(2*np.pi*sigma_envelope**2))*
                          np.exp(-(t-omega_n)**2 /(sigma_envelope)**2) * 
                          np.sin(2*np.pi*2000*t + phase))*np.exp(-t/T_2)
            
    pulse = 0
    for n in range(0, T1+10):

            pulse  = pulse + fct_pulse_sequence(n, T_2_star, T_2, T_E, Amp)
   
    pulse_FFT = fft(pulse) / len(pulse) 
    x = np.size(pulse_FFT)
    SAMPLE_PERIOD =0.00001
    f = np.fft.fftfreq(len(pulse_FFT), SAMPLE_PERIOD)[:x//2] 
    y2 = np.abs(pulse_FFT[0:x//2])
   
    
   
    
    plt.plot(t,pulse)
    plt.xlabel('Time (s)',fontsize='x-large')
    plt.ylabel('signal amplitude (a.u)',fontsize='x-large')
    plt.show()

    plt.plot(f,y2)
    plt.xlim(1990,2010)
    plt.show()

if __name__ == '__main__':
    main()