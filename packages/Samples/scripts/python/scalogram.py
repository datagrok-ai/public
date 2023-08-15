#name: Scalogram Python
#description: CWT (Morlet wavelet) scalogram plot
#reference: https://en.wikipedia.org/wiki/Scaleogram
#language: python
#environment: scalogram
#tags: demo, viewers, hide-suggestions
#sample: sensors/ecg.csv
#input: dataframe data [Input data table]
#input: column signal {type:numerical; allowNulls:false} [Column with signal]
#input: double sampleRate = 256.0 [Signal sample rate, in Hz]
#input: int w0 = 4 [Tradeoff between time and frequency resolution]
#input: bool removeDc = TRUE [Flag to force DC component removal]
#output: graphics [Scalogram plot]

import numpy as np
import matplotlib.pyplot as plt
from obspy.imaging.cm import obspy_sequential
from obspy.signal.tf_misfit import cwt

npts = len(data.index)
dt = 1.0 / sampleRate
t = np.linspace(0.0, dt * npts, npts)
f_min = 1.0
f_max = sampleRate / 10.0

signal = data[[signal]].values.ravel()

scalogram = cwt(signal, dt, w0, f_min, f_max, wl='morlet')
if removeDc:
    signal = signal - np.mean(signal)

fig = plt.figure()
ax = fig.add_subplot(111)

x, y = np.meshgrid(t, np.logspace(np.log10(f_min), np.log10(f_max), scalogram.shape[0]))

ax.pcolormesh(x, y, np.abs(scalogram), cmap=obspy_sequential)
ax.set_xlabel("Time, s")
ax.set_ylabel("Frequency, Hz")
ax.set_yscale('log')
ax.set_ylim(f_min, f_max)
plt.show()
