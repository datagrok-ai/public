export function getRelevantMethods(signalType) {
  const dspPackageFilters = [
    'Moving Average Filter',
    'Exponential Filter',
    'Min Max Normalization',
    'Z-score Normalization',
    'Box Cox Transform',
    'Get Trend',
    'Detrend',
    'Fourier Filter',
    'Spectral Density',
    'Subsample',
    'Averaging Downsampling'
  ];
  const pyphysioFilters = [
    'IIRFilter',
    'FIRFilter',
    'normalize',
    'resample',
    'KalmanFilter',
    'ImputeNAN',
    'RemoveSpikes',
    'ConvolutionalFilter'
  ];
  let commonFilters = dspPackageFilters.concat(pyphysioFilters);
  let commonExtractors = ['Local energy'];
  let commonIndicators = [];
  switch (signalType) {
    case 'ECG':
      return {
        filters: commonFilters,
        extractors: commonExtractors.concat([
          'Beat from ECG'
        ]),
        indicators: commonIndicators.concat([
          'HRV time domain',
          'HRV frequency domain',
          'HRV nonlinear domain'
        ])
      };
    case 'EDA':
      return {
        filters: commonFilters.concat([
          'DenoiseEDA'
        ]),
        extractors: commonExtractors.concat([
          'Phasic estimation'
        ]),
        indicators: commonIndicators
      };
    case 'Accelerometer':
      return {
        filters: commonFilters,
        extractors: commonExtractors,
        indicators: commonIndicators
      };
    case 'EMG':
      return {
        filters: commonFilters,
        extractors: commonExtractors,
        indicators: commonIndicators
      };
    case 'EEG':
      return {
        filters: commonFilters,
        extractors: commonExtractors,
        indicators: commonIndicators
      };
    case 'ABP':
      return {
        filters: commonFilters,
        extractors: commonExtractors.concat([
          'BeatFromBP'
        ]),
        indicators: commonIndicators
      };
    case 'BVP(PPG)':
      return {
        filters: commonFilters,
        extractors: commonExtractors.concat([
          'BeatFromBP'
        ]),
        indicators: commonIndicators
      };
    case 'Respiration':
      return {
        filters: commonFilters,
        extractors: commonExtractors,
        indicators: commonIndicators
      };
  }
}