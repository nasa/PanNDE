
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy import signal
import numpy as np
import math

def smoothPulse(t,pulseWidth=1.0):
  return 0.5*(1.0-np.cos(2.0*math.pi*t/pulseWidth))*(t<pulseWidth)
  
# Wavespeed = 6313 #m/s Titanium
# cellperwave = 20 # cells per wave length  

# block_dims = [0.75e-3, 0.75e-3, 0.5e-3] #m
# ds = [2.0e-6, 2.0e-6, 2.0e-6] #m
# ToF = 2*block_dims[2]/Wavespeed
# pulseWidth = 0.2 *ToF 

# tStart = 0.0
# dt = 0.1e-9
# tStop = 400.0e-9

#### Aluminum

lam = 51.74993908e9 # lame 1
mu = 26.664117020e9 # lame 2
rho = 2780 # density
CFL_ideal=0.85

C = np.zeros((21))
C[0] = lam+2*mu
C[1] = lam
C[2] = lam
C[6] = lam+2*mu
C[7] = lam
C[11] = lam+2*mu
C[15] = mu
C[18] = mu
C[20] = mu

cellperwave = 10 # cells per wave length  

block_dims = [7.5e-2, 2.0e-2, 2.0e-3] #m
ds = [40.0e-6, 40.0e-6, 40.0e-6] #m

tStart = 0.0
dt = 2e-9
tStop = 60.0e-6

Wavespeed = np.sqrt(C[0]/rho)
ToF = 2*block_dims[2]/Wavespeed
pulseWidth = 0.2 *ToF 

tVec = np.arange(tStart, tStop, dt)
wf_Tend = np.argmin(np.abs(tVec-pulseWidth))
wf_tVec = tVec[:wf_Tend]
wf=smoothPulse(tVec, pulseWidth)

N = wf.size
dt = tVec[1] - tVec[0]
fVec = np.fft.fftfreq(N, dt)

freqAmp = np.fft.fft(wf)

half = int(len(fVec)/2)

y_spl = UnivariateSpline(fVec[:half],abs(freqAmp[:half]),s=0,k=4)
y_spl_1d = y_spl.derivative(n=1)
y_spl_2d = y_spl.derivative(n=2)

xmax = 200.0e6
fVec_dense = np.arange(fVec[:half].min(),xmax , 100)

fmax = fVec_dense[np.argmax(y_spl_2d(fVec_dense))]

wavelength = Wavespeed/fmax

min_ds = wavelength/cellperwave

dt_max = 1/(Wavespeed)/(1/ds[0]+1/ds[1]+1/ds[2]) 
dt_ideal = CFL_ideal/(Wavespeed)/(1/ds[0]+1/ds[1]+1/ds[2]) 
CLF_used = dt*Wavespeed*(1/ds[0]+1/ds[1]+1/ds[2])

print("Stiffness Matrix")
print(C)

print("Vbulk = %f m/s" %(Wavespeed))

print("Max dt = %0.3e s, Ideal dt = %0.3e s, using %0.3e s" %(dt_max, dt_ideal, dt))

print("Courant-Friedrichs-Lewy number = %f, should be <1.0" % CLF_used)

print("Max Frequency is %.2E Hz" %(fmax))

print("Max cell size is %.2E m" %(min_ds))

print("Min sim time is %.2E secs" %(1.2*ToF))

print("Block size %dx%dx%d elements" %(int(block_dims[0]/ds[0]), int(block_dims[1]/ds[1]), int(block_dims[2]/ds[2])))


# Windowed 
domFreq = 250e3   # Hz
numCycles = 3
pulseLen = (numCycles/domFreq) # Pulse length (seconds) given # of cycles
pulseSteps = int(np.ceil(pulseLen/dt)) # # of pulse time steps 

amp = 1
wf_tVec = np.arange(0, (pulseSteps)*dt, dt)
wf = (amp*np.sin(wf_tVec*domFreq*2*np.pi)).astype(np.double) 
window = signal.windows.hann(pulseSteps)
wf = np.multiply(wf, window) 
wf = wf/wf.max()


