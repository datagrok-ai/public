import * as grok from "datagrok-api/grok";

async function hrvTimeDomain(t, parameters) {
  return await grok.functions.call('BioSignals:hrvTimeDomain',
    {
      'signal_values': t,
      'step': parameters['step'],
      'width': parameters['width']
    });
}

async function hrvFrequencyDomain(t, parameters) {
  return await grok.functions.call('BioSignals:hrvFrequencyDomain',
    {
      'signal_values': t,
      'resampling_frequency': parameters['inputSamplingFrequency'],
      'step': parameters['step'],
      'width': parameters['width']
    });
}

async function hrvNonlinearDomain(t, parameters) {
  return await grok.functions.call('BioSignals:hrvNonlinearDomain',
    {
      'signal_values': t,
      'step': parameters['step'],
      'width': parameters['width']
    });
}

export async function applyIndicator(t, parameters) {
  switch (parameters['type']) {
    case 'HRV time domain':
      return hrvTimeDomain(t, parameters);
    case 'HRV frequency domain':
      return hrvFrequencyDomain(t, parameters);
    case 'HRV nonlinear domain':
      return hrvNonlinearDomain(t, parameters);
  }
}