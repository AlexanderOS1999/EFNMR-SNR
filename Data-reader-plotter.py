"""
Created on Tue Nov 24 14:10:10 2020

@authors: ewank
"""

import numpy as np
import matplotlib.pyplot as plt

def lim_calc(x , y , width , sensitivity):
    # Function analyses two arrays and returns appropriate x limits based on y variation
    # Width controls how far either side of useful data to plot
    # Sensitivity (< 1) controls what y value relative to y_max is useful data
    # Calc max
    buzz1 , buzz2 , i = 0 , 0 , 0
    ma = np.amax(y)
    while buzz1 == 0:
        if y[i] > sensitivity * ma:
            buzz1 = 1
        else:
            i += 1
    # Calc min
    j = len(y) - 1
    while buzz2 == 0:
        if y[j] > sensitivity * ma:
            buzz2 = 1
        else:
            j -= 1
            
    diff = x[j] - x[i]
    lower = x[i] - (width * diff)
    upper = x[j] + (width * diff)
    return lower , upper

def read_in(file):
    x , y = [] , []
    readFile = open(file,'r')
    readFile.readline()
    for line in readFile:
        splitUp = line.split(',')
        x.append(int(splitUp[0]))
        y.append(int(splitUp[1]))
    readFile.close()
    return x , y

# Fourier Transform function: returns arrays representing Fourier Transform of input arrays
def FT(x , y):
    y_FFT = np.fft.fft(y) / len(y)  # Dividing by length of data to preserve magnitude
    freqs = np.fft.fftfreq(len(y) , (np.amax(x) / len(x)))[:len(y)//2]  # Generate frequencies using sample rate
    y_FFT_y = 2 * np.abs(y_FFT[0:len(y)//2])    # Generating y axis , doubled (modulus)
    return freqs , y_FFT_y

# Reading in FID data from csv file
x , amp = read_in('FIDdata.csv')

# Set sample rate
sample_rate = 0.00001

# Set time scale based on sample rate
t  = np.arange(0 , len(amp) * sample_rate , sample_rate)
amp = np.array(amp)
amp = amp - np.mean(amp)
amp = 150 * amp/np.amax(amp)

# Plot data directly: FID
plt.plot(t , amp , 'k')
plt.title('Simulated FID with noise (imported)' , fontsize='x-large')
plt.ylabel('signal amplitude (a.u)' , fontsize='x-large')
plt.xlabel('Time (s)' , fontsize = 'x-large')
plt.axhline(y = 0 , color = 'black')
plt.grid(True)
plt.show()

# Calculate FT of input data
freqs , amp_FFT_y = FT(t , amp)

# Plot FT
plt.plot(freqs , amp_FFT_y , 'k')
plt.xlabel('Frequency (Hz)' , fontsize='x-large')
plt.ylabel('signal amplitude (a.u)' , fontsize='x-large')
plt.title('Simulated spectrum (imported)')
# Set limits
lower, upper = lim_calc(freqs , amp_FFT_y , 0.1 , 0.02)
plt.xlim(lower , upper)
plt.grid(True)
plt.show()

Larmor_peak = freqs [ np.argmax(amp_FFT_y) ]
print("\nLarmor peak at %5.0f Hz" % Larmor_peak)