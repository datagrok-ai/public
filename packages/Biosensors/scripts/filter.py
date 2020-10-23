#name: filtersPyphysio
#language: python
#input: dataframe ecg_data
#input: int fsamp = 2048
#input: string signalType = "ecg" {choices : ["ecg"]}
#input: dataframe paramsT
#output: graphics fig

# import packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# import the Signal classes
import pyphysio as ph


# convert to numpy
ecg_data = np.array(ecg_data.ecg_data)

# convert to signal class
ecg = ph.EvenlySignal(values = ecg_data, sampling_freq = fsamp, signal_type = signalType)

for i in range(0,len(paramsT)):
    if paramsT['filter'][i] == 'IIR':
        ecg = ph.IIRFilter(fp = paramsT['fp'][i], fs = paramsT['fs'][i], ftype = paramsT['ftype'][i])(ecg)
    if paramsT['filter'][i] == 'normalize':
        ecg = ph.Normalize(norm_method=paramsT['normMethod'][i])(ecg)
    if paramsT['filter'][i] == 'resample':
        ecg = ecg.resample(fout=paramsT['fout'][i], kind=paramsT['kind'][i])
        fsamp = 4096

fig = ecg.plot()