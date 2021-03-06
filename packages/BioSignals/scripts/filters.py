#name: filters
#language: python
#input: dataframe data
#input: int fsamp
#input: dataframe paramsT
#output: dataframe newDf

import numpy as np
import pyphysio as ph


sig = ph.EvenlySignal(values = np.array(data.iloc[:,0]), sampling_freq = fsamp)

i = len(paramsT) - 1
if paramsT['type'][i] == 'IIR':
    sig = ph.IIRFilter(fp = paramsT['fp'][i], fs = paramsT['fs'][i], ftype = paramsT['ftype'][i])(sig)
elif paramsT['type'][i] == 'FIR':
    # It doesn't work on the platform, and works in Jupyter notebook works, but returns:
    # Anaconda\lib\site-packages\scipy\signal\windows\windows.py:1214: RuntimeWarning: invalid value encountered in true_divide special.i0(beta))
    sig = ph.FIRFilter(fp = [paramsT['fp'][i]], fs = [paramsT['fs'][i]])(sig)
elif paramsT['type'][i] == 'normalize':
    sig = ph.Normalize(norm_method=paramsT['normMethod'][i])(sig)
elif paramsT['type'][i] == 'resample':
    sig = sig.resample(fout=paramsT['fout'][i], kind=paramsT['kind'][i])
elif paramsT['type'][i] == 'KalmanFilter':
    sig = ph.KalmanFilter(R=paramsT['R'][i], ratio=paramsT['ratio'][i])(sig)
elif paramsT['type'][i] == 'ImputeNAN':
    sig = ph.ImputeNAN(win_len=paramsT['win_len'][i], allnan=paramsT['allnan'][i])(sig)
elif paramsT['type'][i] == 'RemoveSpikes':
    sig = ph.RemoveSpikes(K=paramsT['K'][i], N=int(paramsT['N'][i]), dilate=paramsT['dilate'][i], D=paramsT['D'][i], method=paramsT['method'][i])(sig)
elif paramsT['type'][i] == 'DenoiseEDA':
    sig = ph.DenoiseEDA(win_len=paramsT['win_len'][i], threshold=paramsT['threshold'][i])(sig)
elif paramsT['type'][i] == 'ConvolutionalFilter':
    sig = ph.ConvolutionalFilter(win_len=paramsT['win_len'][i], irftype=str(paramsT['irftype'][i]))(sig)

newDf = pd.DataFrame({'time':range(0,len(sig)),'sig':sig})