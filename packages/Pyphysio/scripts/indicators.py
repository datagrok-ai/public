#name: indicatorsPyphysio
#language: python
#input: dataframe ecg_data
#input: int fsamp = 2048
#input: string signalType = "ecg" {choices : ["ecg"]}
#input: dataframe paramsT
#input: string info = "Beat from ECG" {choices : ["Beat from ECG"]}
#input: string preset = "HRV" {choices : ["HRV"]}
#output: dataframe FD_HRV_df

# import packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# import the Signal classes
import pyphysio as ph
import pandas as pd

# convert to numpy
ecg_data = np.array(ecg_data.ecg_data)

# create label
label = np.zeros(1200)
label[300:600] = 1
label[900:1200] = 2
label = ph.EvenlySignal(label, sampling_freq = 10, signal_type = 'label')

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

if info == 'Beat from ECG':
    extracted = ph.BeatFromECG()(ecg)

if preset == 'HRV':
    HRV_FD = ph.preset_hrv_fd()
    fixed_length = ph.FixedSegments(step = 5, width = 10, labels = label)
    FD_HRV_ind, col_names = ph.fmap(fixed_length, ph.preset_hrv_fd(), extracted.resample(4))
    FD_HRV_df = pd.DataFrame(FD_HRV_ind, columns=col_names)