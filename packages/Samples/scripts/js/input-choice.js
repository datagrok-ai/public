//name: inputChoice
//help-url: https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation
//language: javascript
//input: string fromList = "France" {choices: ['France', 'Italy', 'Germany']}
//input: string fromFile = "France" {choices: OpenFile("System:AppData/Samples/countries.csv")}
//input: string fromQuery = "France" {choices: Samples:countries}
//output: string result

result = fromList + ' - ' + fromFile + ' - ' + fromQuery;