#name: filtersPyphysio
#language: python
#input: dataframe data
#input: int fsamp
#input: string signalType
#input: dataframe paramsT
#output: dataframe newDf

# import packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# import the Signal classes
import pyphysio as ph


# convert to numpy
# ecg_data = np.array(ecg_data.ecg_data)
data = np.array(data.iloc[:,0])

# convert to signal class
sig = ph.EvenlySignal(values = data, sampling_freq = fsamp, signal_type = signalType)

for i in range(0,len(paramsT)):
    if paramsT['filter'][i] == 'IIR':
        sig = ph.IIRFilter(fp = paramsT['fp'][i], fs = paramsT['fs'][i], ftype = paramsT['ftype'][i])(sig)
    if paramsT['filter'][i] == 'normalize':
        sig = ph.Normalize(norm_method=paramsT['normMethod'][i])(sig)
    if paramsT['filter'][i] == 'resample':
        sig = sig.resample(fout=paramsT['fout'][i], kind=paramsT['kind'][i])

#fig = sig.plot()
newDf = pd.DataFrame({'time':range(0,len(sig)),'sig':sig})