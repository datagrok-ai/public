// Evaluate a parameter's dynamic sources without building a form:
// choices (with the params they depend on), suggestions for the typed text,
// and a computed default.

grok.functions.register({
  signature: 'List<String> demoCities(String region)',
  run: (region) => [`${region}-north`, `${region}-south`],
});

let order = grok.functions.register({
  signature: 'string demoOrder(string region, string city, int qty)',
  run: (region, city, qty) => `${region}/${city} x${qty}`,
});

order.inputs.find((p) => p.name === 'city').options['choices'] = 'demoCities';
order.inputs.find((p) => p.name === 'qty').options['default'] = '2 + 2';

let call = order.prepare({region: 'EU'});
let choices = await call.evalParamChoices('city');
grok.shell.info(`cities: ${choices.items}, depends on: ${choices.dependsOn}`);
grok.shell.info(`default qty: ${await call.evalParamDefault('qty')}`);
