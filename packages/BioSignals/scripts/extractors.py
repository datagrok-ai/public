#name: extractors
#language: python
#input: dataframe data
#input: int fsamp
#input: dataframe paramsT
#output: dataframe newDf

import numpy as np
import pyphysio as ph


sig = ph.EvenlySignal(values = np.array(data.iloc[:,0]), sampling_freq = fsamp)

i = len(paramsT) - 1

if paramsT['type'][i] == 'Beat from ECG':

    extracted = ph.BeatFromECG(
        bpm_max=paramsT['bpm_max'][i],
        delta=paramsT['delta'][i],
        k=paramsT['k'][i]
    )(sig)

    newDf = pd.DataFrame({
        'Index': range(len(extracted)),
        'RR intervals': extracted
    })

elif paramsT['type'][i] == 'Phasic estimation':

    extracted = ph.DriverEstim(
        t1=paramsT['t1'][i],
        t2=paramsT['t2'][i]
    )(sig)

    phasic, tonic, _ = ph.PhasicEstim(
        delta=paramsT['delta'][i],
        grid_size=paramsT['grid_size'][i],
        pre_max=paramsT['pre_max'][i],
        post_max=paramsT['post_max'][i]
    )(extracted)

    newDf = pd.DataFrame({
        'time': range(len(phasic)),
        'phasic': phasic,
        'tonic': tonic
    })

elif paramsT['type'][i] == 'Local energy':

    extracted = ph.Energy(
        win_len=paramsT['win_len'][i],
        win_step=paramsT['win_step'][i]
    )(sig)

    newDf = pd.DataFrame({
        'Index': range(len(extracted)),
        'extracted': extracted
    })

elif paramsT['type'][i] == 'BeatFromBP':

    extracted = ph.BeatFromBP(
        bpm_max=paramsT['bpm_max'][i],
        win_pre=paramsT['win_pre'][i],
        win_post=paramsT['win_post'][i]
    )(sig)

    newDf = pd.DataFrame({
        'Index': range(len(extracted)),
        'Interbeat intervals': extracted
    })