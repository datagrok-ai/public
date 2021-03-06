#name: indicators
#language: python
#input: dataframe data
#input: int fsamp
#input: dataframe paramsT
#input: string info
#input: string preset
#output: dataframe out

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

if info == 'Beat from ECG':
    extracted = ph.BeatFromECG()(sig)
elif info == 'Phasic estimation':
    extracted = ph.DriverEstim()(sig)
elif info == 'Local energy':
    extracted = ph.Energy(win_len=2, win_step=2)(sig)
elif info == 'BeatFromBP':
    extracted = ph.BeatFromBP()(sig)

if preset == 'HRV time domain':

    hrv_indicators = [
        ph.Mean(name='RRmean'),
        ph.StDev(name='RRstd'),
        ph.RMSSD(name='rmsSD'),
        ph.SDSD(name='sdsd'),
        ph.Triang(name='Triang'),
        ph.TINN(name='TINN')
    ]

    fixed_length = ph.FixedSegments(
        step = paramsT['step'][i],
        width = paramsT['width'][i]
    )

    indicators, col_names = ph.fmap(
        fixed_length,
        hrv_indicators,
        extracted
    )

    out = pd.DataFrame({
        'time': range(len(indicators)),
        'RRmean': indicators[:, np.where(col_names == 'RRmean')[0]].ravel(),
        'RRstd': indicators[:, np.where(col_names == 'RRstd')[0]].ravel(),
        'rmsSD': indicators[:, np.where(col_names == 'rmsSD')[0]].ravel(),
        'sdsd': indicators[:, np.where(col_names == 'sdsd')[0]].ravel(),
        'Triang': indicators[:, np.where(col_names == 'Triang')[0]].ravel(),
        'TINN': indicators[:, np.where(col_names == 'TINN')[0]].ravel()
    })

elif preset == 'HRV frequency domain':

    FD_HRV_ind, col_names = ph.fmap(
        fixed_length,
        ph.preset_hrv_fd(),
        extracted.resample(paramsT['resampling_frequency'][i])
    )

    out = pd.DataFrame({
        'time': range(len(FD_HRV_ind)),
        'VLF_Pow': FD_HRV_ind[:, 3],
        'LF_Pow': FD_HRV_ind[:, 4],
        'HF_Pow': FD_HRV_ind[:, 5],
        'Total_Pow': FD_HRV_ind[:, 6]
    })

elif preset == 'HRV nonlinear domain':

    hrv_indicators = [
        ph.PoincareSD1(name='SD1'),
        ph.PoincareSD2(name='SD2'),
        ph.PoincareSD1SD2(name='SD1/SD2')
    ]

    fixed_length = ph.FixedSegments(
        step = paramsT['step'][i],
        width = paramsT['width'][i]
    )

    indicators, col_names = ph.fmap(
        fixed_length,
        hrv_indicators,
        extracted
    )

    out = pd.DataFrame({
        'time': range(len(indicators)),
        'SD1': indicators[:, np.where(col_names == 'SD1')[0]].ravel(),
        'SD2': indicators[:, np.where(col_names == 'SD2')[0]].ravel(),
        'SD1/SD2': indicators[:, np.where(col_names == 'SD1/SD2')[0]].ravel()
    })