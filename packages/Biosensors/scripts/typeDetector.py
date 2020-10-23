#name: typeDetector
#language: python
#input: dataframe dat
#input: int fsamp
#output: signalType

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

n_peaks = 10
x = np.array(dat.iloc[:,0]).ravel()
threshold = (x.max() - x.min())/2
peaks, _ = find_peaks(x, height=threshold)

if(peaks >= n_peaks*0.7 and peaks <= n_peaks*1.3):
    signalType = 'ecg'
else:
    signalType = 'eda'
