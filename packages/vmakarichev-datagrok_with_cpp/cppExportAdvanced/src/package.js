/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

import { callWasm } from '../wasm/callWasm';

//tags: init
export async function init() {
  await initAnalyzer();
}

//top-menu: Tools | Data Science | PCA by Eigen
//tag: ML
//name: pca
//input: dataframe table
//input: column_list columns
//input: int componentsCount
//output: dataframe result 
export function pca(table, columns, componentsCount) {
  return callWasm(Analyzer, 'principalComponentAnalysis', [columns, componentsCount]);
}

//name: error
//input: dataframe df
//input: column col1
//input: column col2
//output: double mad 
export function error(df, col1, col2) {
  return callWasm(Analyzer, 'error', [col1, col2]);
}

//top-menu: Tools | Data Science | MVA (PLS) by Eigen
//name: pls
//tag: ML
//input: dataframe table
//input: column_list features
//input: column predict
//input: int componentsCount
export function pls(table, features, predict, componentsCount) {
  let callOutput = callWasm(Analyzer, 'partialLeastSquareRegression', [features, predict, componentsCount]);

  // CREATING VISUALIZATION

  let start = new Date().getTime();
  
  let dfView = grok.shell.getTableView(table.name);

  // 1. Predicted vs Reference scatter plot

  let prediction = callOutput[0];
  prediction.name = predict.name + '(predicted)';

  let dfReferencePrediction = DG.DataFrame.fromColumns([predict, prediction]);
  dfReferencePrediction.name = 'Reference vs. Predicted';
  
  dfView.addViewer(DG.Viewer.scatterPlot(dfReferencePrediction, 
    { title: dfReferencePrediction.name,
      x: predict.name,
      y: prediction.name,
      showRegressionLine: true,
      markerType: 'circle'
     }));

  // 2. Regression Coefficients Bar Chart
  let regressionCoefficients = callOutput[1];
  regressionCoefficients.name = 'regression coefficient';

  let namesOfPredictors = [];
  for(let col of features)
    namesOfPredictors.push(col.name); 
  
  let  predictorNamesColumn = DG.Column.fromStrings('feature', namesOfPredictors);  

  let dfRegrCoefs = DG.DataFrame.fromColumns([predictorNamesColumn, regressionCoefficients]);
  dfRegrCoefs.name = 'Regression Coefficients';
    
  dfView.addViewer(DG.Viewer.barChart(dfRegrCoefs, 
    {title: dfRegrCoefs.name, split: 'feature', 
     value: 'regression coefficient', valueAggrType: 'avg'}));
  

  // 3. Scores Scatter Plot

  let scoresColumns = [];
  
  let xScores = callOutput[2];
  for(let i = 0; i < xScores.length; i++) {
    xScores[i].name = `x.score.t${i+1}`;
    scoresColumns.push(xScores[i]);
  }

  let yScores = callOutput[3];
  for(let i = 0; i < yScores.length; i++) {
    yScores[i].name = `y.score.u${i+1}`;
    scoresColumns.push(yScores[i]);
  }  

  let scores = DG.DataFrame.fromColumns(scoresColumns);
  scores.name = 'Scores';
  //grok.shell.addTableView(scores);

  dfView.addViewer(DG.Viewer.scatterPlot(scores, 
    { title: scores.name,
      x: xScores[0].name,
      y: yScores[0].name,      
      markerType: 'circle'
     }));

  // 4. Loading Scatter Plot
  let loadingCols = [];

  let loadingLabels = [];
  for(let col of features)
    loadingLabels.push(col.name);

  loadingCols.push(DG.Column.fromStrings('labels', loadingLabels));

  let xLoadings = callOutput[4];
  for(let i = 0; i < xLoadings.length; i++) {
    xLoadings[i].name = `x.loading.p${i+1}`;
    loadingCols.push(xLoadings[i]);
  }

  let dfLoadings = DG.DataFrame.fromColumns(loadingCols);
  dfLoadings.name = 'Loadings';
  
  dfView.addViewer(DG.Viewer.scatterPlot(dfLoadings, 
    { title: dfLoadings.name,
      x: xLoadings[0].name,
      y: xLoadings[0].name,      
      markerType: 'circle',
      labels: 'labels'
     }));

  let finish = new Date().getTime();

  console.log(`Time for creating viewers is ${finish - start} ms.`);
}

