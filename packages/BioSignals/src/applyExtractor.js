import * as grok from "datagrok-api/grok";

async function BeatfromECG(data, column, parameters) {
  return await grok.functions.call('BioSignals:BeatfromECG',
    {
      'data': data,
      'column': column,
      'sampling_frequency': parameters['inputSamplingFrequency'],
      'bpm_max': parameters['bpm_max'],
      'delta': parameters['delta'],
      'k': parameters['k']
    });
}

async function PhasicEstimation(data, column, parameters) {
  return await grok.functions.call('BioSignals:PhasicEstimation',
    {
      'data': data,
      'column': column,
      'sampling_frequency': parameters['inputSamplingFrequency'],
      't1': parameters['t1'],
      't2': parameters['t2'],
      'delta': parameters['delta'],
      'grid_size': parameters['grid_size'],
      'win_pre': parameters['win_pre'],
      'win_post': parameters['win_post']
    });
}

async function LocalEnergy(data, column, parameters) {
  return await grok.functions.call('BioSignals:LocalEnergy',
    {
      'data': data,
      'column': column,
      'sampling_frequency': parameters['inputSamplingFrequency'],
      'win_len': parameters['win_len'],
      'win_step': parameters['win_step']
    });
}

async function BeatFromBP(data, column, parameters) {
  return await grok.functions.call('BioSignals:BeatFromBP',
    {
      'data': data,
      'column': column,
      'sampling_frequency': parameters['inputSamplingFrequency'],
      'bpm_max': parameters['bpm_max'],
      'win_pre': parameters['win_pre'],
      'win_post': parameters['win_post']
    });
}

export async function applyExtractor(t, column, parameters) {
  switch (parameters['type']) {
    case 'LocalEnergy':
      return LocalEnergy(t, column, parameters);
    case 'BeatFromECG':
      return BeatfromECG(t, column, parameters);
    case 'PhasicEstimation':
      return PhasicEstimation(t, column, parameters);
    case 'BeatFromBP':
      return BeatFromBP(t, column, parameters);
  }
}