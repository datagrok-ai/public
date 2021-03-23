import * as ui from "datagrok-api/ui";

export function getFilterParameters(filterType) {
    switch (filterType) {
        case 'IIRFilter':
            let passFrequency = ui.floatInput('Pass frequency', '');
            let stopFrequency = ui.floatInput('Stop frequency', '');
            let ftype = ui.choiceInput('Filter type', '', ['butter', 'cheby1', 'cheby2', 'ellip']);
            return {'fp': passFrequency, 'fs': stopFrequency, 'ftype': ftype};
        case 'FIRFilter':
            let passFreq = ui.floatInput('Pass frequency', '');
            let stopFreq = ui.floatInput('Stop frequency', '');
            let windowType = ui.choiceInput('Window type', 'hamming', ['hamming']);
            return {'fp': passFreq, 'fs': stopFreq, 'windowType': windowType};
        case 'normalize':
            let normMethod = ui.choiceInput('norm_method', '', ['mean', 'standard', 'min', 'maxmin', 'custom']);
            return {'normMethod': normMethod};
        case 'resample':
            let fout = ui.intInput('Output sampling frequency', '');
            let kind = ui.choiceInput('Interpolation method', '', ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic']);
            return {'fout': fout, 'kind': kind};
        case 'KalmanFilter':
            let r = ui.floatInput('R', '');
            r.setTooltip("R should be positive");
            let ratio = ui.floatInput('Ratio', '');
            ratio.setTooltip("Ratio should be >1");
            return {'R': r, 'ratio': ratio};
        case 'ImputeNAN':
            let winLen = ui.floatInput('Window Length', '');
            winLen.setTooltip('Window Length should be positive');
            let allNan = ui.choiceInput('All NaN', '', ['zeros', 'nan']);
            return {'win_len': winLen, 'allnan': allNan};
        case 'RemoveSpikes':
            let K = ui.floatInput('K', 2);
            K.setTooltip("K should be positive");
            let N = ui.intInput('N', 1);
            N.setTooltip('N should be positive integer');
            let dilate = ui.floatInput('Dilate', 0);
            dilate.setTooltip('dilate should be >= 0.0');
            let D = ui.floatInput('D', 0.95);
            D.setTooltip('D should be >= 0.0');
            let method = ui.choiceInput('Method', '', ['linear', 'step']);
            return {'K': K, 'N': N, 'dilate': dilate, 'D': D, 'method': method};
        case 'DenoiseEDA':
            let winLen3 = ui.floatInput('Window Length', 2);
            winLen3.setTooltip('Window Length should be positive');
            let threshold = ui.floatInput('Threshold', '');
            threshold.setTooltip('Threshold should be positive');
            return {'win_len': winLen3, 'threshold': threshold};
        case 'ConvolutionalFilter':
            let irftype = ui.choiceInput('irftype', '', ['gauss', 'rect', 'triang', 'dgauss', 'custom']);
            irftype.setTooltip('Impulse response function (IRF)');
            let winLen1 = ui.floatInput('Window Length', 2);
            winLen1.setTooltip('Duration of the generated IRF in seconds (if irftype is not \'custom\')');
            return {'win_len': winLen1, 'irftype': irftype};
        case 'Moving Average Filter':
            let winLen2 = ui.floatInput('Window Length', '');
            winLen2.setTooltip('Order of filtering');
            return {'win_len': winLen2};
        case 'Exponential Filter':
            let filterRatio = ui.floatInput('Filter ratio', '');
            return {'filter_ratio': filterRatio};
        case 'Min Max Normalization':
            return {};
        case 'Z-score Normalization':
            return {};
        case 'Box Cox Transform':
            let lambda = ui.floatInput('lambda', '');
            let ofset = ui.floatInput('ofset', '');
            return {'lambda': lambda, 'ofset': ofset};
        case 'Get Trend':
            return {};
        case 'Detrend':
            return {};
        case 'Fourier Filter':
            let lowcut = ui.floatInput('lowcut', '');
            let hicut = ui.floatInput('hicut', '');
            let observationTime1 = ui.floatInput('Observation time', '');
            return {'lowcut': lowcut, 'hicut': hicut, 'observationTime': observationTime1};
        case 'Spectral Density':
            let observationTime = ui.floatInput('Observation time', '');
            return {'observationTime': observationTime};
        case 'Subsample':
            let subsampleSize = ui.floatInput('Subsample size', '');
            let offset1 = ui.floatInput('Offset', '');
            return {'subsampleSize': subsampleSize, 'offset': offset1};
        case 'Averaging Downsampling':
            let windowSize = ui.floatInput('Window size', '');
            let offset = ui.floatInput('Offset', '');
            return {'windowSize': windowSize, 'offset': offset};
    }
}