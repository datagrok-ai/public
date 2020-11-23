#name: importPyphysio
#language: python
#input: dataframe data
#input: int fsamp = 2048
#input: string signalType = "ecg" {choices : ["ecg","eda"]}
#output: dataframe newDf

# import packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pyphysio as ph


# convert to numpy
# ecg_data = np.array(ecg_data.ecg_data)
data = np.array(data.iloc[:,0])

# convert to signal class
sig = ph.EvenlySignal(values = data, sampling_freq = fsamp, signal_type = signalType)

# render
#fig = sig.plot()
newDf = pd.DataFrame({'time':range(0,len(sig)),'sig':sig})