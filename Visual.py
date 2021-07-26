import numpy as np
import matplotlib.pyplot as pl
from matplotlib.widgets import Slider
from scipy.fft import fft, fftfreq

def input(t, f_l):
    f = np.cos(-2*np.pi*f_l*t)
    return f

def output(x, f_l, f_s, W, A):
    f = np.cos( A*W*2*np.pi*f_s*np.sin(2*np.pi*f_s*x) - (2*np.pi*f_l*x) + W*2*np.pi*f_l*n/c)
    return f

n = 1
a = 1
c = 1
f_l = 1
f_s = 0.1
W = 0.1
N = int(1e6)

fig, ax = pl.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
pl.subplots_adjust( left =0.25, bottom = 0.25 )

t=np.linspace(0,2/f_l,1000)
linput, = ax[0].plot(t,input(t,f_l))
lsignal, = ax[1].plot(t,output(t, f_l, f_s, W, a/c))

ax[1].set_xlabel('t (sec)')
ax[1].set_ylabel('Electric field (normalised)')
ax[1].grid()
ax[0].grid()
ax[0].set_title('Input Signal')
ax[1].set_title('Output Signal')

axcolor = 'lightgoldenrodyellow'
axdfreqL = pl.axes([0.25, 0.02, 0.65, 0.03], facecolor=axcolor)
axfreqS = pl.axes([0.25, 0.07, 0.65, 0.03], facecolor=axcolor)
axWidth = pl.axes([0.25, 0.12, 0.65, 0.03], facecolor=axcolor)

sdfreqL = Slider(axdfreqL, 'Freq Light (Hz)', 0, 10, valinit=1)
sfreqS = Slider(axfreqS, 'Freq Sound (Hz)', 0, 5.0, valinit=0.1)
sWidth = Slider(axWidth, 'Width (m)', 0, 5.0, valinit=0.1)

def update(val):
    f_l = sdfreqL.val
    f_s = sfreqS.val
    W = sWidth.val

    linput.set_ydata(input(t,f_l))
    lsignal.set_ydata(output(t, f_l, f_s, W, a/c))
    fig.canvas.draw_idle()


sdfreqL.on_changed(update)
sfreqS.on_changed(update)
sWidth.on_changed(update)

pl.show()