import * as grok from "datagrok-api/grok";

async function hrvTimeDomain(t, column, parameters) {
  return await grok.functions.call('BioSignals:hrvTimeDomain',
    {
      'signal_values': t,
      'column': column,
      'step': parameters['step'],
      'width': parameters['width']
    });
}

async function hrvFrequencyDomain(t, column, parameters) {
  return await grok.functions.call('BioSignals:hrvFrequencyDomain',
    {
      'signal_values': t,
      'column': column,
      'resampling_frequency': parameters['inputSamplingFrequency'],
      'step': parameters['step'],
      'width': parameters['width']
    });
}

async function hrvNonlinearDomain(t, column, parameters) {
  return await grok.functions.call('BioSignals:hrvNonlinearDomain',
    {
      'signal_values': t,
      'column': column,
      'step': parameters['step'],
      'width': parameters['width']
    });
}

export async function applyIndicator(t, column, parameters) {
  switch (parameters['type']) {
    case 'HRV time domain':
      return hrvTimeDomain(t, column, parameters);
    case 'HRV frequency domain':
      return hrvFrequencyDomain(t, column, parameters);
    case 'HRV nonlinear domain':
      return hrvNonlinearDomain(t, column, parameters);
  }
}