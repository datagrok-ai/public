grok.functions.register({
  signature: 'List<String> jsSuggestCountryName(String text)',
  isAsync: true,
  run: async function(text) {
    let response = await fetch('https://restcountries.eu/rest/v2/name/' + text);
    return response.status === 200 ? (await response.json()).map(country => country['name']) : [];
  }
});

grok.functions.register({
  signature: 'List<String> jsVeggies()',
  run: () => ['Artichoke', 'Cucumber', 'Cauliflower', 'Onion']
});

grok.functions.register({
  signature: 'List<String> jsSaltinessRange(double input)',
  run: (input) => input >= 0 && input <= 100 ? null : 'Saltiness is out of range'
});
