#name: indicatorsPyphysio
#language: python
#input: dataframe data
#input: int fsamp
#input: string signalType
#input: dataframe paramsT
#input: string info
#input: string preset
#output: dataframe out

import numpy as np
import pyphysio as ph

data = np.array(data.iloc[:,0])

label = np.zeros(1200)
label[300:600] = 1
label[900:1200] = 2
label = ph.EvenlySignal(label, sampling_freq = 10, signal_type = 'label')

sig = ph.EvenlySignal(values = data, sampling_freq = fsamp, signal_type = signalType)

for i in range(0,len(paramsT)):
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
        sig = ph.ImputeNAN(win_len=paramsT['win_len'][i], allnan=str(paramsT['allnan'][i]))(sig)
    elif paramsT['filter'][i] == 'RemoveSpikes':
        sig = ph.RemoveSpikes(K=paramsT['K'][i], N=int(paramsT['N'][i]), dilate=paramsT['dilate'][i], D=paramsT['D'][i], method=paramsT['method'][i])(sig)
    elif paramsT['filter'][i] == 'DenoiseEDA':
        sig = ph.DenoiseEDA(win_len=paramsT['win_len'][i], threshold=paramsT['threshold'][i])(sig)
    elif paramsT['filter'][i] == 'ConvolutionalFilter':
        sig = ph.ConvolutionalFilter(win_len=paramsT['win_len'][i], irftype=str(paramsT['irftype'][i]))(sig)

if info == 'Beat from ECG':
    extracted = ph.BeatFromECG()(sig)
elif info == 'Phasic estimation':
    extracted = ph.DriverEstim()(sig)
elif info == 'Local energy':
    extracted = ph.Energy(win_len=2, win_step=2)(sig)
elif info == 'BeatFromBP':
    extracted = ph.BeatFromBP()(sig)

if preset == 'HRV time domain':
    hrv_indicators = [ph.Mean(name='RRmean'), ph.StDev(name='RRstd'), ph.RMSSD(name='rmsSD'), ph.SDSD(name='sdsd'), ph.Triang(name='Triang'), ph.TINN(name='TINN')]
    fixed_length = ph.FixedSegments(step = 5, width = 10, labels = label)
    indicators, col_names = ph.fmap(fixed_length, hrv_indicators, extracted)
    RRmean = indicators[:, np.where(col_names == 'RRmean')[0]].ravel()
    RRstd = indicators[:, np.where(col_names == 'RRstd')[0]].ravel()
    rmsSD = indicators[:, np.where(col_names == 'rmsSD')[0]].ravel()
    sdsd = indicators[:, np.where(col_names == 'sdsd')[0]].ravel()
    Triang = indicators[:, np.where(col_names == 'Triang')[0]].ravel()
    tinn = indicators[:, np.where(col_names == 'TINN')[0]].ravel()
    out = pd.DataFrame({'time': range(len(RRmean)), 'RRmean': RRmean, 'RRstd': RRstd, 'rmsSD': rmsSD, 'sdsd': sdsd, 'Triang': Triang, 'TINN': tinn})
elif preset == 'HRV frequency domain':
    FD_HRV_ind, col_names = ph.fmap(fixed_length, ph.preset_hrv_fd(), extracted.resample(4))
    out = pd.DataFrame({'time': range(len(FD_HRV_ind)), 'VLF_Pow': FD_HRV_ind[:, 3], 'LF_Pow': FD_HRV_ind[:, 4],
                        'HF_Pow': FD_HRV_ind[:, 5], 'Total_Pow': FD_HRV_ind[:, 6]})
