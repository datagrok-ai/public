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

async function applyPyphysioFilter(data, paramsT, fsamp) {
  let f = await grok.functions.eval("BioSignals:filters");
  let call = f.prepare({
    'data': data,
    'fsamp': fsamp,
    'paramsT': paramsT
  });
  await call.call();
  return call.getParamValue('newDf');
}

export async function applyFilter(i, col, inputCase, filterInputsList, filterOutputsObj, filterTypesList, filterParametersList, fsamp) {
  let pi = DG.TaskBarProgressIndicator.create('Calculating and plotting filter\'s output...');
  let t;
  let nameOfLastFiltersOutput;
  if (filterInputsList.length === 1) {
    t = DG.DataFrame.fromColumns([inputCase]);
  }
  else if (filterOutputsObj[filterInputsList[i-1].value].columns.byName(filterInputsList[i-1].value)) {
    t = DG.DataFrame.fromColumns([filterOutputsObj[filterInputsList[i-1].value].columns.byName(filterInputsList[i-1].value)]);
  }
  else {
    t = DG.DataFrame.fromColumns([filterOutputsObj[filterInputsList[i-1].value].columns.byName('sig')]);
  }

  let plotFL;
  let currentlyChosenFilterType = filterTypesList[filterTypesList.length-1].value;
  switch (currentlyChosenFilterType) {
    case 'Moving Average Filter':
      await SMA_filter(t, inputCase, filterParametersList[i-1].win_len.value);
      nameOfLastFiltersOutput = inputCase.name + ' SMA Filtered';
      plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
      break;
    case 'Exponential Filter':
      await Exp_filter(t, inputCase, filterParametersList[i-1].filter_ratio.value);
      nameOfLastFiltersOutput = inputCase.name + ' Exponentially Filtered';
      plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
      break;
    case 'Min Max Normalization':
      await MinMax_transform(t, inputCase);
      nameOfLastFiltersOutput = inputCase.name + ' Min Max Normalized';
      plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
      break;
    case 'Z-score Normalization':
      await Zscore_transform(t, inputCase);
      nameOfLastFiltersOutput = inputCase.name + ' Z-score Normalized';
      plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
      break;
    case 'Box Cox Transform':
      await box_cox_transform(t, inputCase, filterParametersList[i-1].lambda.value, filterParametersList[i-1].ofset.value);
      nameOfLastFiltersOutput = inputCase.name + ' Box Cox Transformed';
      plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
      break;
    case 'Get Trend':
      await get_trend(t, inputCase);
      nameOfLastFiltersOutput = inputCase.name + ' Trend';
      plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
      break;
    case 'Detrend':
      await remove_trend(t, inputCase);
      nameOfLastFiltersOutput = inputCase.name + ' Detrended';
      plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
      break;
    case 'Fourier Filter':
      await fourier_filter(t, inputCase, filterParametersList[i-1].lowcut.value, filterParametersList[i-1].hicut.value, filterParametersList[i-1].observationTime.value);
      nameOfLastFiltersOutput = inputCase.name + ' Fourier Filtered (L: ' + filterParametersList[i-1].lowcut.value + '; H: ' + filterParametersList[i-1].hicut.value + ')';
      plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
      break;
    case 'Spectral Density':
      await spectral_density(t, inputCase, filterParametersList[i-1].observationTime.value);
      nameOfLastFiltersOutput = inputCase.name + ' Density';
      plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
      break;
    case 'Subsample':
      await subsample(t, inputCase, filterParametersList[i-1].subsampleSize.value, filterParametersList[i-1].offset.value);
      nameOfLastFiltersOutput = inputCase.name + ' Subsample';
      plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
      break;
    case 'Averaging Downsampling':
      await asample(t, inputCase, filterParametersList[i-1].windowSize.value, filterParametersList[i-1].offset.value);
      nameOfLastFiltersOutput = inputCase.name + ' Subsample';
      plotFL = DG.DataFrame.fromColumns([col, t.columns.byName(nameOfLastFiltersOutput)]);
      break;
    default:
      let paramsT = parametersToDataFrame(filterTypesList, filterParametersList);
      nameOfLastFiltersOutput = 'Output of Filter ' + i + ' (' + filterTypesList[i-1].value + ')';
      plotFL = await applyPyphysioFilter(t, paramsT, fsamp);
  }
  pi.close();
  return [plotFL, nameOfLastFiltersOutput]
}