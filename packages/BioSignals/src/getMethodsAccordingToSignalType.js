export function getMethodsAccordingTo(signalType) {
    const dspPackageMethods = [
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
    const pyphysioMethods = [
        'IIR',
        'FIR',
        'normalize',
        'resample',
        'KalmanFilter',
        'ImputeNAN',
        'RemoveSpikes',
        'ConvolutionalFilter'
    ];
    let commonFilters = dspPackageMethods.concat(pyphysioMethods);
    let commonEstimators = ['Local energy'];
    let commonIndicators = [];
    switch (signalType) {
        case 'ECG':
            return {
                filters: commonFilters,
                estimators: commonEstimators.concat([
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
                estimators: commonEstimators.concat([
                    'Phasic estimation'
                ]),
                indicators: commonIndicators
            };
        case 'Accelerometer':
            return {
                filters: commonFilters,
                estimators: commonEstimators,
                indicators: commonIndicators
            };
        case 'EMG':
            return {
                filters: commonFilters,
                estimators: commonEstimators,
                indicators: commonIndicators
            };
        case 'EEG':
            return {
                filters: commonFilters,
                estimators: commonEstimators,
                indicators: commonIndicators
            };
        case 'ABP':
            return {
              filters: commonFilters,
              estimators: commonEstimators.concat([
                  'BeatFromBP'
              ]),
              indicators: commonIndicators
            };
        case 'BVP(PPG)':
            return {
              filters: commonFilters,
              estimators: commonEstimators.concat([
                  'BeatFromBP'
              ]),
              indicators: commonIndicators
            };
        case 'Respiration':
            return {
              filters: commonFilters,
              estimators: commonEstimators,
              indicators: commonIndicators
            };
    }
}