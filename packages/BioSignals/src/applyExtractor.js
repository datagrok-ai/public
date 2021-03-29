import * as grok from "datagrok-api/grok";

async function BeatfromECG(data, parameters) {
  return await grok.functions.call('BioSignals:BeatfromECG',
    {
      'data': data,
      'sampling_frequency': parameters['inputSamplingFrequency'],
      'bpm_max': parameters['bpm_max'],
      'delta': parameters['delta'],
      'k': parameters['k']
    });
}

async function PhasicEstimation(data, parameters) {
  return await grok.functions.call('BioSignals:PhasicEstimation',
    {
      'data': data,
      'sampling_frequency': parameters['inputSamplingFrequency'],
      't1': parameters['t1'],
      't2': parameters['t2'],
      'delta': parameters['delta'],
      'grid_size': parameters['grid_size'],
      'win_pre': parameters['win_pre'],
      'win_post': parameters['win_post']
    });
}

async function LocalEnergy(data, parameters) {
  return await grok.functions.call('BioSignals:LocalEnergy',
    {
      'data': data,
      'sampling_frequency': parameters['inputSamplingFrequency'],
      'win_len': parameters['win_len'],
      'win_step': parameters['win_step']
    });
}

async function BeatFromBP(data, parameters) {
  return await grok.functions.call('BioSignals:BeatFromBP',
    {
      'data': data,
      'sampling_frequency': parameters['inputSamplingFrequency'],
      'bpm_max': parameters['bpm_max'],
      'win_pre': parameters['win_pre'],
      'win_post': parameters['win_post']
    });
}

export async function applyExtractor(t, parameters) {
  switch (parameters['type']) {
    case 'Local energy':
      return LocalEnergy(t, parameters);
    case 'Beat from ECG':
      return BeatfromECG(t, parameters);
    case 'Phasic estimation':
      return PhasicEstimation(t, parameters);
    case 'BeatFromBP':
      return BeatFromBP(t, parameters);
  }
}