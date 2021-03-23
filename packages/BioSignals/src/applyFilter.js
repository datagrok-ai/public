import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import {parametersToDataFrame} from "./package.js";

async function asample(data, col, windowSize, offset) {
  return await grok.functions.call('Dsp:asample',
    {
      'data': data,
      'col': col,
      'windowSize': windowSize,
      'offset': offset
    });
}

async function subsample(data, col, subsampleSize, offset) {
  return await grok.functions.call('Dsp:subsample',
    {
      'data': data,
      'col': col,
      'subsampleSize': subsampleSize,
      'offset': offset
    });
}

async function spectral_density(data, col, observationTime) {
  return await grok.functions.call('Dsp:spectral_density',
    {
      'data': data,
      'col': col,
      'observationTime': observationTime
    });
}

async function fourier_filter(data, column_to_filter, lowcut, hicut, observationTime) {
  return await grok.functions.call('Dsp:fourier_filter',
    {
      'data': data,
      'column_to_filter': column_to_filter,
      'lowcut': lowcut,
      'hicut': hicut,
      'observationTime': observationTime
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

async function box_cox_transform(data, columnToFilter, lambda, ofset) {
  return await grok.functions.call('Dsp:box_cox_transform',
    {
      'data': data,
      'column_to_filter': columnToFilter,
      'lambda': lambda,
      'ofset': ofset
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

async function Exp_filter(data, columnToFilter, filterRatio) {
  return await grok.functions.call('Dsp:Exp_filter',
    {
      'data': data,
      'column_to_filter': columnToFilter,
      'filter_ratio': filterRatio
    });
}

async function SMA_filter(data, columnToFilter, windowSize) {
  return await grok.functions.call('Dsp:SMA_filter',
    {
      'data': data,
      'column_to_filter': columnToFilter,
      'window_size': windowSize
    });
}

async function ConvolutionalFilter(data, paramsT, samplingFrequency) {
  return await grok.functions.call('BioSignals:ConvolutionalFilter',
    {
      'inputSignals': data,
      'sampling_frequency': samplingFrequency,
      'parameters': paramsT
    });
}

async function DenoiseEDA(data, paramsT, samplingFrequency) {
  return await grok.functions.call('BioSignals:DenoiseEDA',
    {
      'input_signals': data,
      'sampling_frequency': samplingFrequency,
      'parameters': paramsT
    });
}

async function IIRFilter(data, paramsT, samplingFrequency) {
  return await grok.functions.call('BioSignals:IIRFilter',
    {
      'input_signals': data,
      'sampling_frequency': samplingFrequency,
      'parameters': paramsT
    });
}

async function FIRFilter(data, paramsT, samplingFrequency) {
  return await grok.functions.call('BioSignals:FIRFilter',
    {
      'input_signals': data,
      'sampling_frequency': samplingFrequency,
      'parameters': paramsT
    });
}

async function normalize(data, paramsT, samplingFrequency) {
  return await grok.functions.call('BioSignals:normalize',
    {
      'input_signals': data,
      'sampling_frequency': samplingFrequency,
      'parameters': paramsT
    });
}

async function resample(data, paramsT, samplingFrequency) {
  return await grok.functions.call('BioSignals:resample',
    {
      'input_signals': data,
      'sampling_frequency': samplingFrequency,
      'parameters': paramsT
    });
}

async function KalmanFilter(data, paramsT, samplingFrequency) {
  return await grok.functions.call('BioSignals:KalmanFilter',
    {
      'input_signals': data,
      'sampling_frequency': samplingFrequency,
      'parameters': paramsT
    });
}

async function ImputeNAN(data, paramsT, samplingFrequency) {
  return await grok.functions.call('BioSignals:ImputeNAN',
    {
      'input_signals': data,
      'sampling_frequency': samplingFrequency,
      'parameters': paramsT
    });
}

async function RemoveSpikes(data, paramsT, samplingFrequency) {
  return await grok.functions.call('BioSignals:RemoveSpikes',
    {
      'input_signals': data,
      'sampling_frequency': samplingFrequency,
      'parameters': paramsT
    });
}

export async function applyFilter(i, col, inputCase, filterInputsList, filterOutputsObj, filterTypesList, filterParametersList, fsamp) {
  let t, nameOfLastFiltersOutput, plotFL;
  if (filterInputsList.length === 1) {
    t = DG.DataFrame.fromColumns([inputCase]);
  }
  else if (filterOutputsObj[filterInputsList[i-1].value].columns.byName(filterInputsList[i-1].value)) {
    t = DG.DataFrame.fromColumns([filterOutputsObj[filterInputsList[i-1].value].columns.byName(filterInputsList[i-1].value)]);
  }
  else {
    t = DG.DataFrame.fromColumns([filterOutputsObj[filterInputsList[i-1].value].columns.byName('sig')]);
  }

  const parametersTable = parametersToDataFrame(filterTypesList, filterParametersList);
  let currentlyChosenFilterType = filterTypesList[filterTypesList.length-1].value;
  switch (currentlyChosenFilterType) {
    case 'Moving Average Filter':
      await SMA_filter(t, inputCase, filterParametersList[i-1].win_len.value);
      nameOfLastFiltersOutput = inputCase.name + ' SMA Filtered';
      plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
      return [plotFL, nameOfLastFiltersOutput];
    case 'Exponential Filter':
      await Exp_filter(t, inputCase, filterParametersList[i-1].filter_ratio.value);
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
      await box_cox_transform(t, inputCase, filterParametersList[i-1].lambda.value, filterParametersList[i-1].ofset.value);
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
      await fourier_filter(t, inputCase, filterParametersList[i-1].lowcut.value, filterParametersList[i-1].hicut.value, filterParametersList[i-1].observationTime.value);
      nameOfLastFiltersOutput = inputCase.name + ' Fourier Filtered (L: ' + filterParametersList[i-1].lowcut.value + '; H: ' + filterParametersList[i-1].hicut.value + ')';
      plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
      return [plotFL, nameOfLastFiltersOutput];
    case 'Spectral Density':
      await spectral_density(t, inputCase, filterParametersList[i-1].observationTime.value);
      nameOfLastFiltersOutput = inputCase.name + ' Density';
      plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
      return [plotFL, nameOfLastFiltersOutput];
    case 'Subsample':
      await subsample(t, inputCase, filterParametersList[i-1].subsampleSize.value, filterParametersList[i-1].offset.value);
      nameOfLastFiltersOutput = inputCase.name + ' Subsample';
      plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
      return [plotFL, nameOfLastFiltersOutput];
    case 'Averaging Downsampling':
      await asample(t, inputCase, filterParametersList[i-1].windowSize.value, filterParametersList[i-1].offset.value);
      nameOfLastFiltersOutput = inputCase.name + ' Subsample';
      plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
      return [plotFL, nameOfLastFiltersOutput];
    case 'ConvolutionalFilter':
      nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + filterTypesList[i-1].value + ')';
      plotFL = await ConvolutionalFilter(t, parametersTable, fsamp);
      return [plotFL, nameOfLastFiltersOutput];
    case 'DenoiseEDA':
      nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + filterTypesList[i-1].value + ')';
      plotFL = await DenoiseEDA(t, parametersTable, fsamp);
      return [plotFL, nameOfLastFiltersOutput];
    case 'IIRFilter':
      nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + filterTypesList[i-1].value + ')';
      plotFL = await IIRFilter(t, parametersTable, fsamp);
      return [plotFL, nameOfLastFiltersOutput];
    case 'FIRFilter':
      nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + filterTypesList[i-1].value + ')';
      plotFL = await FIRFilter(t, parametersTable, fsamp);
      return [plotFL, nameOfLastFiltersOutput];
    case 'normalize':
      nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + filterTypesList[i-1].value + ')';
      plotFL = await normalize(t, parametersTable, fsamp);
      return [plotFL, nameOfLastFiltersOutput];
    case 'resample':
      nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + filterTypesList[i-1].value + ')';
      plotFL = await resample(t, parametersTable, fsamp);
      return [plotFL, nameOfLastFiltersOutput];
    case 'KalmanFilter':
      nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + filterTypesList[i-1].value + ')';
      plotFL = await KalmanFilter(t, parametersTable, fsamp);
      return [plotFL, nameOfLastFiltersOutput];
    case 'ImputeNAN':
      nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + filterTypesList[i-1].value + ')';
      plotFL = await ImputeNAN(t, parametersTable, fsamp);
      return [plotFL, nameOfLastFiltersOutput];
    case 'RemoveSpikes':
      nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + filterTypesList[i-1].value + ')';
      plotFL = await RemoveSpikes(t, parametersTable, fsamp);
      return [plotFL, nameOfLastFiltersOutput];
  }
}