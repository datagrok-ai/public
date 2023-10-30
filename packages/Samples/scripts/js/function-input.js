//name: Function Input
//description: Embedding
//help-url: https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation
//language: javascript
//input: dataframe orders {category: Data; editor: Samples:OrdersByEmployee}
//input: int factor = 2 {category: Computation}
//output: int result

result = orders.rowCount * 2;