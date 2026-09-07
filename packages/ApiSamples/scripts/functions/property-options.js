// Property.options is the union of the property's tags and its FuncParam
// options (tags win). For a FuncParam-backed property without tags it is the
// live options map, so writes decorate the registered function.

let brew = grok.functions.register({
  signature: 'string demoBrew(string flavor, double dose)',
  run: (flavor, dose) => `${flavor}: ${dose}`,
});

let dose = brew.inputs.find((p) => p.name === 'dose');
dose.options['units'] = 'mg';
dose.min = 0;
dose.max = 10;

grok.shell.info(`units: ${brew.inputs.find((p) => p.name === 'dose').options['units']}`);
