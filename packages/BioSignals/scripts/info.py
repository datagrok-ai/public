#name: infoPyphysio
#language: python
#input: dataframe data
#input: int fsamp
#input: string signalType
#input: dataframe paramsT
#input: string info
#output: dataframe newDf

import numpy as np
import pyphysio as ph

data = np.array(data.iloc[:,0])

label = np.zeros(1200)
label[300:600] = 1
label[900:1200] = 2
label = ph.EvenlySignal(label, sampling_freq = 10, signal_type = 'label')

sig = ph.EvenlySignal(values = data, sampling_freq = fsamp, signal_type = signalType)

i = len(paramsT) - 1
if paramsT['filter'][i] == 'IIR':
    sig = ph.IIRFilter(fp = paramsT['fp'][i], fs = paramsT['fs'][i], ftype = paramsT['ftype'][i])(sig)
elif paramsT['filter'][i] == 'FIR':
    # It doesn't work on the platform, and works in Jupyter notebook works, but returns:
    # Anaconda\lib\site-packages\scipy\signal\windows\windows.py:1214: RuntimeWarning: invalid value encountered in true_divide special.i0(beta))
    sig = ph.FIRFilter(fp = [paramsT['fp'][i]], fs = [paramsT['fs'][i]])(sig)
elif paramsT['filter'][i] == 'normalize':
    sig = ph.Normalize(norm_method=paramsT['normMethod'][i])(sig)
elif paramsT['filter'][i] == 'resample':
    sig = sig.resample(fout=paramsT['fout'][i], kind=paramsT['kind'][i])
elif paramsT['filter'][i] == 'KalmanFilter':
    sig = ph.KalmanFilter(R=paramsT['R'][i], ratio=paramsT['ratio'][i])(sig)
elif paramsT['filter'][i] == 'ImputeNAN':
    sig = ph.ImputeNAN(win_len=paramsT['win_len'][i], allnan=paramsT['allnan'][i])(sig)
elif paramsT['filter'][i] == 'RemoveSpikes':
    sig = ph.RemoveSpikes(K=paramsT['K'][i], N=int(paramsT['N'][i]), dilate=paramsT['dilate'][i], D=paramsT['D'][i], method=paramsT['method'][i])(sig)
elif paramsT['filter'][i] == 'DenoiseEDA':
    sig = ph.DenoiseEDA(win_len=paramsT['win_len'][i], threshold=paramsT['threshold'][i])(sig)
elif paramsT['filter'][i] == 'ConvolutionalFilter':
    sig = ph.ConvolutionalFilter(win_len=paramsT['win_len'][i], irftype=str(paramsT['irftype'][i]))(sig)

if info == 'Beat from ECG':
    extracted = ph.BeatFromECG()(sig)
    newDf = pd.DataFrame({'time':range(0,len(extracted)),'extracted':extracted})
elif info == 'Phasic estimation' :
    extracted = ph.DriverEstim()(sig)
    phasic, tonic, _ = ph.PhasicEstim(delta=0.02)(extracted)
    newDf = pd.DataFrame({'time':range(0,len(phasic)),'phasic':phasic,'tonic':tonic})
elif info == 'Local energy':
    extracted = ph.Energy(win_len=2, win_step=2)(sig)
    newDf = pd.DataFrame({'time':range(0,len(extracted)),'extracted':extracted})
elif info == 'BeatFromBP':
    extracted = ph.BeatFromBP()(sig)
    newDf = pd.DataFrame({'time':range(0,len(extracted)),'extracted':extracted})