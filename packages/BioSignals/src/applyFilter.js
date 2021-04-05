import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

async function asample(data, col, parameters) {
  return await grok.functions.call('Dsp:asample',
    {
      'data': data,
      'col': col,
      'windowSize': parameters['windowSize'],
      'offset': parameters['offset']
    });
}

async function subsample(data, col, parameters) {
  return await grok.functions.call('Dsp:subsample',
    {
      'data': data,
      'col': col,
      'subsampleSize': parameters['subsampleSize'],
      'offset': parameters['offset']
    });
}

async function spectral_density(data, col, parameters) {
  return await grok.functions.call('Dsp:spectral_density',
    {
      'data': data,
      'col': col,
      'observationTime': parameters['observationTime']
    });
}

async function fourier_filter(data, column_to_filter, parameters) {
  return await grok.functions.call('Dsp:fourier_filter',
    {
      'data': data,
      'column_to_filter': column_to_filter,
      'lowcut': parameters['lowcut'],
      'hicut': parameters['hicut'],
      'observationTime': parameters['observationTime']
    });
}

async function remove_trend(data, col) {
  return await grok.functions.call('Dsp:remove_trend',
    {
      'data': data,
      'col': col
    });
}

async function get_trend(data, col) {
  return await grok.functions.call('Dsp:get_trend',
    {
      'data': data,
      'col': col
    });
}

async function box_cox_transform(data, columnToFilter, parameters) {
  return await grok.functions.call('Dsp:box_cox_transform',
    {
      'data': data,
      'column_to_filter': columnToFilter,
      'lambda': parameters['lambda'],
      'ofset': parameters['ofset']
    });
}

async function Zscore_transform(data, columnToFilter) {
  return await grok.functions.call('Dsp:Zscore_transform',
    {
      'data': data,
      'column_to_filter': columnToFilter
    });
}

async function MinMax_transform(data, columnToFilter) {
  return await grok.functions.call('Dsp:MinMax_transform',
    {
      'data': data,
      'column_to_filter': columnToFilter
    });
}

async function Exp_filter(data, columnToFilter, parameters) {
  return await grok.functions.call('Dsp:Exp_filter',
    {
      'data': data,
      'column_to_filter': columnToFilter,
      'filter_ratio': parameters['filter_ratio']
    });
}

async function SMA_filter(data, columnToFilter, parameters) {
  return await grok.functions.call('Dsp:SMA_filter',
    {
      'data': data,
      'column_to_filter': columnToFilter,
      'window_size': parameters['win_len']
    });
}

async function ConvolutionalFilter(data, columnToFilter, parameters) {
  return await grok.functions.call('BioSignals:ConvolutionalFilter',
    {
      'inputSignals': data,
      'column': columnToFilter,
      'sampling_frequency': parameters['inputSamplingFrequency'],
      'win_len': parameters['win_len'],
      'irftype': parameters['irftype']
    });
}

async function DenoiseEDA(data, columnToFilter, parameters) {
  return await grok.functions.call('BioSignals:DenoiseEDA',
    {
      'input_signals': data,
      'column': columnToFilter,
      'sampling_frequency': parameters['inputSamplingFrequency'],
      'win_len': parameters['win_len'],
      'threshold': parameters['threshold']
    });
}

async function IIRFilter(data, columnToFilter, parameters) {
  return await grok.functions.call('BioSignals:IIRFilter',
    {
      'input_signals': data,
      'column': columnToFilter,
      'sampling_frequency': parameters['inputSamplingFrequency'],
      'fs': parameters['fs'],
      'fp': parameters['fp'],
      'ftype': parameters['ftype']
    });
}

async function FIRFilter(data, columnToFilter, parameters) {
  return await grok.functions.call('BioSignals:FIRFilter',
    {
      'input_signals': data,
      'column': columnToFilter,
      'sampling_frequency': parameters['inputSamplingFrequency'],
      'fs': parameters['fs'],
      'fp': parameters['fp']
    });
}

async function normalize(data, columnToFilter, parameters) {
  return await grok.functions.call('BioSignals:Normalize',
    {
      'input_signals': data,
      'column': columnToFilter,
      'sampling_frequency': parameters['inputSamplingFrequency'],
      'norm_method': parameters['normMethod']
    });
}

async function resample(data, columnToFilter, parameters) {
  return await grok.functions.call('BioSignals:Resample',
    {
      'input_signals': data,
      'column': columnToFilter,
      'sampling_frequency': parameters['inputSamplingFrequency'],
      'fout': parameters['fout'],
      'kind': parameters['kind']
    });
}

async function KalmanFilter(data, columnToFilter, parameters) {
  return await grok.functions.call('BioSignals:KalmanFilter',
    {
      'input_signals': data,
      'column': columnToFilter,
      'sampling_frequency': parameters['inputSamplingFrequency'],
      'R': parameters['R'],
      'ratio': parameters['ratio']
    });
}

async function ImputeNAN(data, columnToFilter, parameters) {
  return await grok.functions.call('BioSignals:ImputeNAN',
    {
      'input_signals': data,
      'column': columnToFilter,
      'sampling_frequency': parameters['inputSamplingFrequency'],
      'allnan': parameters['allnan']
    });
}

async function RemoveSpikes(data, columnToFilter, parameters) {
  return await grok.functions.call('BioSignals:RemoveSpikes',
    {
      'input_signals': data,
      'column': columnToFilter,
      'sampling_frequency': parameters['inputSamplingFrequency'],
      'K': parameters['K'],
      'N': parameters['N'],
      'dilate': parameters['dilate'],
      'D': parameters['D'],
      'method': parameters['method']
    });
}

export async function applyFilter(t, parameters, i) {
  let pi = DG.TaskBarProgressIndicator.create('Calculating filter\'s output...');
  let nameOfLastFiltersOutput, plotFL, inputCase;
  for (let column of t.columns.byTags({'selectedFilterInput': undefined}))
    inputCase = t.columns.byName(column.name);
  const col = DG.Column.fromList('int', 'time', Array(inputCase.length).fill().map((_, idx) => idx));
  try {
    switch (parameters['type']) {
      case 'Moving Average Filter':
        await SMA_filter(t, inputCase, parameters);
        nameOfLastFiltersOutput = inputCase.name + ' SMA Filtered';
        plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
        return [plotFL, nameOfLastFiltersOutput];
      case 'Exponential Filter':
        await Exp_filter(t, inputCase, parameters);
        nameOfLastFiltersOutput = inputCase.name + ' Exponentially Filtered';
        plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
        return [plotFL, nameOfLastFiltersOutput];
      case 'Min Max Normalization':
        await MinMax_transform(t, inputCase);
        nameOfLastFiltersOutput = inputCase.name + ' Min Max Normalized';
        plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
        return [plotFL, nameOfLastFiltersOutput];
      case 'Z-score Normalization':
        await Zscore_transform(t, inputCase);
        nameOfLastFiltersOutput = inputCase.name + ' Z-score Normalized';
        plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
        return [plotFL, nameOfLastFiltersOutput];
      case 'Box Cox Transform':
        await box_cox_transform(t, inputCase, parameters);
        nameOfLastFiltersOutput = inputCase.name + ' Box Cox Transformed';
        plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
        return [plotFL, nameOfLastFiltersOutput];
      case 'Get Trend':
        await get_trend(t, inputCase);
        nameOfLastFiltersOutput = inputCase.name + ' Trend';
        plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
        return [plotFL, nameOfLastFiltersOutput];
      case 'Detrend':
        await remove_trend(t, inputCase);
        nameOfLastFiltersOutput = inputCase.name + ' Detrended';
        plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
        return [plotFL, nameOfLastFiltersOutput];
      case 'Fourier Filter':
        await fourier_filter(t, inputCase, parameters);
        nameOfLastFiltersOutput = inputCase.name + ' Fourier Filtered (L: ' + parameters['lowcut'] + '; H: ' + parameters['hicut'] + ')';
        plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
        return [plotFL, nameOfLastFiltersOutput];
      case 'Spectral Density':
        await spectral_density(t, inputCase, parameters);
        nameOfLastFiltersOutput = inputCase.name + ' Density';
        plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
        return [plotFL, nameOfLastFiltersOutput];
      case 'Subsample':
        await subsample(t, inputCase, parameters);
        nameOfLastFiltersOutput = inputCase.name + ' Subsample';
        plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
        return [plotFL, nameOfLastFiltersOutput];
      case 'Averaging Downsampling':
        await asample(t, inputCase, parameters);
        nameOfLastFiltersOutput = inputCase.name + ' Subsample';
        plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
        return [plotFL, nameOfLastFiltersOutput];
      case 'ConvolutionalFilter':
        nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + parameters['type'] + ')';
        plotFL = await ConvolutionalFilter(t, inputCase, parameters);
        return [plotFL, nameOfLastFiltersOutput];
      case 'DenoiseEDA':
        nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + parameters['type'] + ')';
        plotFL = await DenoiseEDA(t, inputCase, parameters);
        return [plotFL, nameOfLastFiltersOutput];
      case 'IIRFilter':
        nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + parameters['type'] + ')';
        plotFL = await IIRFilter(t, inputCase, parameters);
        return [plotFL, nameOfLastFiltersOutput];
      case 'FIRFilter':
        nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + parameters['type'] + ')';
        plotFL = await FIRFilter(t, inputCase, parameters);
        return [plotFL, nameOfLastFiltersOutput];
      case 'Normalize':
        nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + parameters['type'] + ')';
        plotFL = await normalize(t, inputCase, parameters);
        return [plotFL, nameOfLastFiltersOutput];
      case 'Resample':
        nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + parameters['type'] + ')';
        plotFL = await resample(t, inputCase, parameters);
        return [plotFL, nameOfLastFiltersOutput];
      case 'KalmanFilter':
        nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + parameters['type'] + ')';
        plotFL = await KalmanFilter(t, inputCase, parameters);
        return [plotFL, nameOfLastFiltersOutput];
      case 'ImputeNAN':
        nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + parameters['type'] + ')';
        plotFL = await ImputeNAN(t, inputCase, parameters);
        return [plotFL, nameOfLastFiltersOutput];
      case 'RemoveSpikes':
        nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + parameters['type'] + ')';
        plotFL = await RemoveSpikes(t, inputCase, parameters);
        return [plotFL, nameOfLastFiltersOutput];
    }
  } catch (e) {
    pi.close();
    alert(e);
    throw e;
  }
  pi.close();
}