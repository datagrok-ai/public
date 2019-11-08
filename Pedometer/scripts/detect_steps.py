#name: Detect Steps
#description: Detects positions of steps based on accelerometer data
#language: python
#tags: template, demo, accelerometer
#sample: accelerometer.csv
#input: dataframe accel [Accelerometry data table]
#input: column x {semType: Accelerometer-X} [X axis]
#input: column y {semType: Accelerometer-Y} [Y axis]
#input: column z {semType: Accelerometer-Z} [Z axis]
#input: double sample_rate = 32 [Sample rate, in Hz]
#output: dataframe steps {action:join(accel)} [Steps positions]

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, savgol_filter, filtfilt, butter

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter(order, cutoff / (fs / 2.0))
    return filtfilt(b, a, data)

num_samples = len(accel.index)
x = accel[x]
y = accel[y]
z = accel[z]

force = x**2 + y**2 + z**2
force[np.isnan(force)] = 0.0
filtered = butter_lowpass_filter(force, 3.0, sample_rate)
filtered = filtered - savgol_filter(filtered, int(sample_rate * 20) + 1, 2)
min_peak_height = np.std(np.abs(filtered)) / 3.0
idx_peaks, _ = find_peaks(filtered, height=min_peak_height, distance=sample_rate/2.0)
peaks = np.zeros(num_samples, dtype=np.int)
peaks[idx_peaks] = 1

steps = pd.DataFrame(peaks, columns=['steps'])
