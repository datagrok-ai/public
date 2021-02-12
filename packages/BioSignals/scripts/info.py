#name: infoPyphysio
#language: python
#input: dataframe data
#input: int fsamp = 2048
#input: string signalType = "ecg" {choices : ["ecg","eda"]}
#input: dataframe paramsT
#input: string info = "Beat from ECG" {choices : ["Beat from ECG"]}
#output: dataframe newDf

# import packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# import the Signal classes
import pyphysio as ph
import pandas as pd

# convert to numpy
data = np.array(data.iloc[:,0])

# create label
label = np.zeros(1200)
label[300:600] = 1
label[900:1200] = 2
label = ph.EvenlySignal(label, sampling_freq = 10, signal_type = 'label')

# convert to signal class
sig = ph.EvenlySignal(values = data, sampling_freq = fsamp, signal_type = signalType)

for i in range(0,len(paramsT)):
    if paramsT['filter'][i] == 'IIR':
        sig = ph.IIRFilter(fp = paramsT['fp'][i], fs = paramsT['fs'][i], ftype = paramsT['ftype'][i])(sig)
    if paramsT['filter'][i] == 'normalize':
        sig = ph.Normalize(norm_method=paramsT['normMethod'][i])(sig)
    if paramsT['filter'][i] == 'resample':
        sig = sig.resample(fout=paramsT['fout'][i], kind=paramsT['kind'][i])
        fsamp = 4096

if(info == 'Beat from ECG'):
    extracted = ph.BeatFromECG()(sig)
    #fig, axs = plt.subplots(2)
    #fig.suptitle('Beat from ECG')
    #axs[2].plot(extracted,'.-')
    #axs[1].plot(label,'.-')
    newDf = pd.DataFrame({'time':range(0,len(extracted)),'extracted':extracted})
if(info == 'Phasic estimation'):
    extracted = ph.DriverEstim()(sig)
    phasic, tonic, _ = ph.PhasicEstim(delta=0.02)(extracted)
    newDf = pd.DataFrame({'time':range(0,len(phasic)),'phasic':phasic,'tonic':tonic})