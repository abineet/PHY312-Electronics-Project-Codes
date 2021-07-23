import numpy as np
import matplotlib.pyplot as pl
from matplotlib.widgets import Slider
from scipy.fft import fft, fftfreq
from scipy import interpolate
from scipy.special import jv

def input(t, f_l):
    f = np.cos(-2*np.pi*f_l*t)
    return f

def output(x, f_l, f_s, W, A):
    f = np.cos( A*W*2*np.pi*f_l*np.sin(2*np.pi*f_s*x) - (2*np.pi*f_l*x) )
    return f

a = 2e-5
c = 3e8
f_l = 5e14
f_s = 5e9
W = np.linspace(0,4,100)
N = int(1e6)
T = 1.0/(2*(f_l+10*f_s))

x = np.linspace(0.0,T*N,N,endpoint=False)
req_v = np.linspace(f_l, f_l+4*f_s, 5)

y = np.array([[0,0,0,0,0]])
for w in W:
    signal = output(x, f_l, f_s, w*1e-2, a/c)
    ft = 2.0/N * np.abs(fft(signal)[0:N//2])
    v = fftfreq(N,T)[0:N//2]
    f = interpolate.interp1d(v, ft)
    y = np.concatenate((y,f(req_v).reshape((1,5))), axis = 0)

fig, ax = pl.subplots(nrows=5, ncols=1, sharex=True)
color = ['red', 'green', 'blue', 'purple', 'black']
for i in range(5):
    ax[i].plot(W,jv(i,a*W*1e-2*2*np.pi*f_l/c)**2,label = 'Expected (Bessel function)\u00b2 of order '+str(i), color = color[i], linestyle = '--', linewidth = 1)
    ax[i].scatter(W, y[1:,i]**2, label = 'Data points obtained for q = '+str(i), color = color[i], marker = '^', s=20)
    ax[i].set_ylim((0,1))
    ax[i].grid()
    ax[i].legend()
ax[4].set_xlabel('Width (cm)')
ax[2].set_ylabel('Intensity normalised')
pl.show()