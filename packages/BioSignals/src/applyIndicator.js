import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";

async function hrvTimeDomain(t, indicatorParametersDF) {
  return await grok.functions.call('BioSignals:hrvTimeDomain',
    {
      'signal_values': t,
      'step': indicatorParametersDF.columns.byName('step').max,
      'width': indicatorParametersDF.columns.byName('width').max
    });
}

async function hrvFrequencyDomain(t, indicatorParametersDF) {
  return await grok.functions.call('BioSignals:hrvFrequencyDomain',
    {
      'signal_values': t,
      'resampling_frequency': indicatorParametersDF.columns.byName('resampling_frequency').max,
      'step': indicatorParametersDF.columns.byName('step').max,
      'width': indicatorParametersDF.columns.byName('width').max
    });
}

async function hrvNonlinearDomain(t, indicatorParametersDF) {
  return await grok.functions.call('BioSignals:hrvNonlinearDomain',
    {
      'signal_values': t,
      'step': indicatorParametersDF.columns.byName('step').max,
      'width': indicatorParametersDF.columns.byName('width').max
    });
}

export async function applyIndicator(t, indicatorParametersDF, indicatorType) {
  switch (indicatorType) {
    case 'HRV time domain':
      return hrvTimeDomain(t, indicatorParametersDF);
    case 'HRV frequency domain':
      return hrvFrequencyDomain(t, indicatorParametersDF);
    case 'HRV nonlinear domain':
      return hrvNonlinearDomain(t, indicatorParametersDF);
  }
}