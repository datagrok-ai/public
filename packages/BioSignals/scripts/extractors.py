#name: extractors
#language: python
#input: dataframe data
#input: int fsamp
#input: string info
#output: dataframe newDf

import numpy as np
import pyphysio as ph


sig = ph.EvenlySignal(values = np.array(data.iloc[:,0]), sampling_freq = fsamp)

if info == 'Beat from ECG':
    extracted = ph.BeatFromECG()(sig)
    newDf = pd.DataFrame({
        'Index': range(len(extracted)),
        'RR intervals': extracted
    })
elif info == 'Phasic estimation':
    extracted = ph.DriverEstim()(sig)
    phasic, tonic, _ = ph.PhasicEstim(delta=0.02)(extracted)
    newDf = pd.DataFrame({
        'time': range(len(phasic)),
        'phasic': phasic,
        'tonic': tonic
    })
elif info == 'Local energy':
    extracted = ph.Energy(win_len=2, win_step=2)(sig)
    newDf = pd.DataFrame({
        'Index': range(len(extracted)),
        'extracted': extracted
    })
elif info == 'BeatFromBP':
    extracted = ph.BeatFromBP()(sig)
    newDf = pd.DataFrame({
        'Index': range(len(extracted)),
        'Interbeat intervals': extracted
    })