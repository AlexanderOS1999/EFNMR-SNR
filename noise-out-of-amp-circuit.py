# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 19:32:57 2020

@author: luis1
"""

import numpy as np
import matplotlib.pyplot as plt

K_B = 1.381*1E-23
T = 298
R = 20
Delta_F = 2000
noise_preamp = 1.5*1E-9
noise_filter = 3.8*1E-9
preamp_gain = 2000
amp_gain = 500
n = 3
Signal_amp = 1E-6

def Noise_circuit (Delta_F, R, T, n, amp_gain):
    global K_B, preamp_gain
    noise_J = np.sqrt(4*K_B*T*R*Delta_F)
    noise = (noise_J + noise_preamp*np.sqrt(Delta_F))*preamp_gain
    noise2 = n*noise_filter*np.sqrt(Delta_F)
    noise_tot = (noise + noise2)*amp_gain
    return noise_tot

def signal (amp_gain):
    global preamp_gain
    signal = Signal_amp*preamp_gain*amp_gain
    return signal
    


signal_tot = signal(amp_gain)

print (signal_tot, 'V')

noise_out  = Noise_circuit (Delta_F, R, T, n, amp_gain)

print (noise_out, 'V')

Signal_error = (noise_out / signal_tot)*100

print (Signal_error, '%')

Signal_unc = signal_tot*Signal_error/100

print (signal_tot, '+-' ,Signal_unc, 'V')









