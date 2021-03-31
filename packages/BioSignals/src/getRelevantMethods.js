async function getScriptsNames(tag) {
  const scripts = await grok.dapi.scripts.filter(tag).list();
  return await scripts.map((script) => script.name);
}

export async function getRelevantMethods(signalType) {
  const dspPackageFilters = [
    'Moving Average Filter',
    'Exponential Filter',
    'Min Max Normalization',
    'Z-score Normalization',
    'Box Cox Transform',
    'Get Trend',
    'Detrend',
    'Fourier Filter',
    'Spectral Density',
    'Subsample',
    'Averaging Downsampling'
  ];
  const pyphysioFilters = await getScriptsNames('#filters');
  return {
    filters: dspPackageFilters.concat(pyphysioFilters),
    extractors: await getScriptsNames('#extractors'),
    indicators: await getScriptsNames('#indicators')
  };
}