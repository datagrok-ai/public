#name: typeDetector
#language: python
#input: dataframe dat
#input: int npeaks
#input: int fsamp
#output: string signalType

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

x = np.array(dat.iloc[:,0]).ravel()
threshold = (x.max() - x.min())/2
peaks, _ = find_peaks(x, height=threshold)

if(peaks >= npeaks*0.7 and peaks <= npeaks*1.3):
    signalType = 'ecg'
else:
    signalType = 'eda'
