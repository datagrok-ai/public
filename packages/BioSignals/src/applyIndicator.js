import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

async function hrvTimeDomain(t, column, parameters) {
  return await grok.functions.call('BioSignals:HRVTimeDomain',
    {
      'signal_values': t,
      'column': column,
      'step': parameters['step'],
      'width': parameters['width']
    });
}

async function hrvFrequencyDomain(t, column, parameters) {
  return await grok.functions.call('BioSignals:HRVFrequencyDomain',
    {
      'signal_values': t,
      'column': column,
      'resampling_frequency': parameters['inputSamplingFrequency'],
      'step': parameters['step'],
      'width': parameters['width']
    });
}

async function hrvNonlinearDomain(t, column, parameters) {
  return await grok.functions.call('BioSignals:HRVNonlinearDomain',
    {
      'signal_values': t,
      'column': column,
      'step': parameters['step'],
      'width': parameters['width']
    });
}

export async function applyIndicator(t, column, parameters) {
  let pi = DG.TaskBarProgressIndicator.create('Calculating indicators...');
  try {
    switch (parameters['type']) {
      case 'HRVTimeDomain':
        return hrvTimeDomain(t, column, parameters);
      case 'HRVFrequencyDomain':
        return hrvFrequencyDomain(t, column, parameters);
      case 'HRVNonlinearDomain':
        return hrvNonlinearDomain(t, column, parameters);
    }
  } catch (e) {
    pi.close();
    alert(e);
    throw e;
  }
  pi.close();
}