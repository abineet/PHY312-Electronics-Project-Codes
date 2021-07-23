#%%
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.widgets import Slider
from scipy.fft import fft, fftfreq

def input(t, f_l):
    f = np.cos(-2*np.pi*f_l*t)
    return f

def output(x, f_l, f_s, W, A):
    f = np.cos( A*W*2*np.pi*f_l*np.sin(2*np.pi*f_s*x) - (2*np.pi*f_l*x) )
    return f

a = 2e-5
c = 3e8
f_l0 = 5e14
f_s0 = 5e9
W0 = 0.01
N = int(1e6)
T = 1.0/(2*(f_l0+10*f_s0))

x = np.linspace(0.0,T*N,N,endpoint=False)
signal = output(x, f_l0, f_s0, W0, a/c)
ft = 2.0/N * np.abs(fft(signal)[0:N//2])
v = fftfreq(N,T)[0:N//2]

fig, ax = pl.subplots(nrows=3, ncols=1, sharex=False, sharey=False)
pl.subplots_adjust( left =0.25, bottom = 0.2, top = 1 )

t=np.linspace(0,2/f_l0,100)
linput, = ax[0].plot(t,input(t,f_l0), label = 'Input Signal', color = 'red')
lsignal, = ax[1].plot(t,output(t, f_l0, f_s0, W0, a/c), label = 'Output Signal', color = 'green')
lftrans, = ax[2].plot(v, ft, label = 'Fourier Transform')

ax[2].set_xlim((f_l0-10*f_s0, f_l0+10*f_s0))
ax[2].set(xlabel = 'Frequency (units of 1e11 Hz, centered at 5e14 Hz)')
ax[1].set(xlabel = 'time (1e-15 sec)', ylabel = 'Electric field normalised')
ax[0].set(ylabel = 'Electric field normalised')
for i in range(3):
    ax[i].grid()
    ax[i].legend(loc = 'upper right')

axcolor = 'lightgoldenrodyellow'
axdfreqL = pl.axes([0.25, 0.02, 0.65, 0.03], facecolor=axcolor)
axfreqS = pl.axes([0.25, 0.06, 0.65, 0.03], facecolor=axcolor)
axWidth = pl.axes([0.25, 0.10, 0.65, 0.03], facecolor=axcolor)

sdfreqL = Slider(axdfreqL, 'Change Freq Light (GHz)', -5, 5, valinit=0)
sfreqS = Slider(axfreqS, 'Freq Sound (GHz)', 0, 10.0, valinit=5)
sWidth = Slider(axWidth, 'Width (cm)', 0, 2.0, valinit=1)

def update(val):
    f_l = sdfreqL.val*1e9+f_l0
    f_s = sfreqS.val*1e9
    W = sWidth.val*1e-2

    T = 1.0/(2*(f_l+10*f_s))

    x = np.linspace(0.0,T*N,N,endpoint=False)
    signal = output(x, f_l, f_s, W, a/c)
    ft = 2.0/N * np.abs(fft(signal)[0:N//2])
    v = fftfreq(N,T)[0:N//2]

    linput.set_ydata(input(t,f_l))
    lsignal.set_ydata(output(t, f_l, f_s, W, a/c))
    lftrans.set_ydata(ft)
    lftrans.set_xdata(v)
    fig.canvas.draw_idle()


sdfreqL.on_changed(update)
sfreqS.on_changed(update)
sWidth.on_changed(update)

pl.show()