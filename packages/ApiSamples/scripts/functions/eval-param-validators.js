// Named validators on a function parameter, run against the call's current value.
// Passing validators are omitted; false maps to a "didn't pass" message;
// a string result becomes the message. Validator functions must be sync.

grok.functions.register({
  signature: 'string demoZipCode(string s)',
  run: (s) => /^\d{5}$/.test(s) ? null : 'Expected a 5-digit zip code',
});

let shipment = grok.functions.register({
  signature: 'string demoShipment(string zip)',
  run: (zip) => `shipping to ${zip}`,
});

shipment.inputs.find((p) => p.name === 'zip').options['validators'] = ['demoZipCode'];

let call = shipment.prepare({zip: '123'});
let results = await call.evalParamValidators('zip');
grok.shell.info(results.length === 0 ? 'valid' : results[0].message);

call.setParamValue('zip', '12345');
results = await call.evalParamValidators('zip');
grok.shell.info(results.length === 0 ? 'valid' : results[0].message);
