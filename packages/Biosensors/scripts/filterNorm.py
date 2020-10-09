#name: filterNormPyphysio
#language: python
#input: dataframe ecg_data
#input: int fsamp = 2048
#input: string signalType = "ecg" {choices : ["ecg"]}
#input: string norm_method = "standard" {choices : ["standard"]}
#output: graphics plt

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
# normalization : normalize data
ecg = ph.Normalize(norm_method=norm_method)(ecg)

plt = ecg.plot()