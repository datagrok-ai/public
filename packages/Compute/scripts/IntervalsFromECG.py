#name: IntervalsFromECG
#description: Identify the beats in an ECG signal and compute the IBIs
#language: python
#input: double samplingFrequency = 2048
#input: int bpmMax = 120 [Maximal expected heart rate (in beats per minute): (1 < bpm_max <= 400)]
#input: double delta = 0 [Threshold for the peak detection. By default it is computed from the signal (adaptive thresholding)]
#input: double k = 0.7 [Ratio at which the signal range is multiplied (when delta = 0), (0 < k < 1)]
#output: graphics Intervals

import pyphysio as ph
import matplotlib.pyplot as plt
from pyphysio.tests import TestData

ecg_data = TestData.ecg()

ecg = ph.EvenlySignal(values=ecg_data, sampling_freq=samplingFrequency, signal_type = 'ecg')

extracted = ph.BeatFromECG(bpm_max=bpmMax, delta=delta, k=k)(ecg)

Intervals = plt.plot(extracted)
plt.title('RR intervals')
plt.xlabel('Interval index')
plt.ylabel('Seconds')