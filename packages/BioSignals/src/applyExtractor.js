import * as grok from "datagrok-api/grok";

async function BeatfromECG(data, paramsT, samplingFrequency) {
  return await grok.functions.call('BioSignals:BeatfromECG',
    {
      'data': data,
      'sampling_frequency': samplingFrequency,
      'bpm_max': paramsT.columns.byName('bpm_max').max,
      'delta': paramsT.columns.byName('delta').max,
      'k': paramsT.columns.byName('k').max
    });
}

async function PhasicEstimation(data, paramsT, samplingFrequency) {
  return await grok.functions.call('BioSignals:PhasicEstimation',
    {
      'data': data,
      'sampling_frequency': samplingFrequency,
      't1': paramsT.columns.byName('t1').max,
      't2': paramsT.columns.byName('t2').max,
      'delta': paramsT.columns.byName('delta').max,
      'grid_size': paramsT.columns.byName('grid_size').max,
      'win_pre': paramsT.columns.byName('win_pre').max,
      'win_post': paramsT.columns.byName('win_post').max
    });
}

async function LocalEnergy(data, paramsT, samplingFrequency) {
  return await grok.functions.call('BioSignals:LocalEnergy',
    {
      'data': data,
      'sampling_frequency': samplingFrequency,
      'win_len': paramsT.columns.byName('win_len').max,
      'win_step': paramsT.columns.byName('win_step').max
    });
}

async function BeatFromBP(data, paramsT, samplingFrequency) {
  return await grok.functions.call('BioSignals:BeatFromBP',
    {
      'data': data,
      'sampling_frequency': samplingFrequency,
      'bpm_max': paramsT.columns.byName('bpm_max').max,
      'win_pre': paramsT.columns.byName('win_pre').max,
      'win_post': paramsT.columns.byName('win_post').max
    });
}

export async function applyExtractor(t, samplingFrequency, parametersTable) {
  const extractorType = parametersTable.getCol('type').categories[0];
  switch (extractorType) {
    case 'Local energy': return LocalEnergy(t, parametersTable, samplingFrequency);
    case 'Beat from ECG': return BeatfromECG(t, parametersTable, samplingFrequency);
    case 'Phasic estimation': return PhasicEstimation(t, parametersTable, samplingFrequency);
    case 'BeatFromBP': return BeatFromBP(t, parametersTable, samplingFrequency);
  }
}