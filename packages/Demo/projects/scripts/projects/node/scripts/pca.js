//name: PCA Node
//description: Principal Component Analysis
//help-url: https://en.wikipedia.org/wiki/Principal_component_analysis
//language: nodejs
//tags: demo, hide-suggestions
//sample: cars.csv
//input: dataframe table {columns:numerical} [Input data table]
//input: column_list features {type:numerical;allowNulls:false} [Features columns]
//input: int components = 2 [Number of components]
//input: bool center = TRUE [Indicating whether the variables should be shifted to be zero centered]
//input: bool scale = TRUE [Indicating whether the variables should be scaled to have unit variance before the analysis takes place]
//output: dataframe result {action:join(table)} [PCA components]

const PCA = require('ml-pca');
const array = table.toArray();
const pca = new PCA(array, {'scale': scale, 'center': center});
const predicted = pca.predict(array, {'nComponents': components});
let result = new DataFrame(predicted, Array.from({length: components}, (v, i) => 'PCA' + i.toString()));
