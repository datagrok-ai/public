// This is a demo-package that illustartes wasm-functions call:
// in the main stream and in webworkers. 

/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

// Imports for call wasm runtime-system: in the main stream and in webworkers
import { callWasm } from '../wasm/callWasm';
import { getCppInput, getResult } from '../wasm/callWasmForWebWorker';

//tags: init
export async function init() {
  await initCppLib();
}

//name: sum
//input: int a
//input: int b
//output: int res 
export function sum(a, b) {
  return callWasm(CppLib, 'sum', [a, b]);
}

//name: sum
//input: int a
//input: int b
export function sumInWebWorker(a, b) {
  var worker = new Worker(new URL('../wasm/sumWorker.js', import.meta.url));
  worker.postMessage(getCppInput(CppLib['sum'].arguments,[a, b]));
  worker.onmessage = function(e) {
    let output = getResult(CppLib['sum'], e.data);

    // The following is added manually
    alert(`${a} + ${b} = ${output}`);
  }
}

//name: maxFloatCol
//input: dataframe df
//input: column col
//output: double max 
export function maxFloatCol(df, col) {
  return callWasm(CppLib, 'maxFloatCol', [col]);
}

//name: maxFloatCol
//input: dataframe df
//input: column col
export function maxFloatColInWebWorker(df, col) {
  var worker = new Worker(new URL('../wasm/maxFloatColWorker.js', import.meta.url));
  worker.postMessage(getCppInput(CppLib['maxFloatCol'].arguments,[col]));
  worker.onmessage = function(e) {
    let output = getResult(CppLib['maxFloatCol'], e.data);

    // The following is added manually
    alert(`Max of the column '${col.name}' of the dataframe '${df.name}' is ${output}`);
  }
}

//name: maxIntCol
//input: dataframe df
//input: column col
//output: int max 
export function maxIntCol(df, col) {
  return callWasm(CppLib, 'maxIntCol', [col]);
}

//name: maxIntCol
//input: dataframe df
//input: column col
export function maxIntColInWebWorker(df, col) {
  var worker = new Worker(new URL('../wasm/maxIntColWorker.js', import.meta.url));
  worker.postMessage(getCppInput(CppLib['maxIntCol'].arguments,[col]));
  worker.onmessage = function(e) {
    let output = getResult(CppLib['maxIntCol'], e.data);

    // The following is added manually
    alert(`Max of the column '${col.name}' of the dataframe '${df.name}' is ${output}`);
  }
}

//name: addFloatCols
//input: dataframe df
//input: column col1
//input: column col2
//output: column sum
export function addFloatCols(df, col1, col2) {
  // this line was generated automatically
  //return callWasm(CppLib, 'addFloatCols', [col1, col2]);

  // The following is added manually
  let res = callWasm(CppLib, 'addFloatCols', [col1, col2]);  
  res.name = `sum of ${col1.name} and ${col2.name}`;
  df.columns.add(res);
}

//name: addFloatCols
//input: dataframe df
//input: column col1
//input: column col2
export function addFloatColsInWebWorker(df, col1, col2) {
  var worker = new Worker(new URL('../wasm/addFloatColsWorker.js', import.meta.url));
  worker.postMessage(getCppInput(CppLib['addFloatCols'].arguments,[col1, col2]));
  worker.onmessage = function(e) {
    let output = getResult(CppLib['addFloatCols'], e.data);

    // The following is added manually
    output.name = `sum of ${col1.name} and ${col2.name}`;
    df.columns.add(output);
  }
}

//name: addIntCols
//input: dataframe df
//input: column col1
//input: column col2
//output: column sum
export function addIntCols(df, col1, col2) {
  // this line was generated automatically
  //return callWasm(CppLib, 'addIntCols', [col1, col2]);

  // The following is added manually
  let res = callWasm(CppLib, 'addIntCols', [col1, col2]);  
  res.name = `sum of ${col1.name} and ${col2.name}`;
  df.columns.add(res);
}

//name: addIntCols
//input: dataframe df
//input: column col1
//input: column col2
export function addIntColsInWebWorker(df, col1, col2) {
  var worker = new Worker(new URL('../wasm/addIntColsWorker.js', import.meta.url));
  worker.postMessage(getCppInput(CppLib['addIntCols'].arguments,[col1, col2]));
  worker.onmessage = function(e) {
    let output = getResult(CppLib['addIntCols'], e.data);

    // Provide output usage!
    output.name = `sum of ${col1.name} and ${col2.name}`;
    df.columns.add(output);
  }
}

//name: doubledInts
//input: dataframe table
//input: column_list cols
//output: dataframe result 
export function doubledInts(table, cols) {
  // this line was generated automatically
  //return callWasm(CppLib, 'doubledInts', [cols]);

  // The following is added manually
  let res = callWasm(CppLib, 'doubledInts', [cols]);
  res.name = 'Doubled_ints';
  return res;
}

//name: doubledInts
//input: dataframe table
//input: column_list cols
export function doubledIntsInWebWorker(table, cols) {
  var worker = new Worker(new URL('../wasm/doubledIntsWorker.js', import.meta.url));
  worker.postMessage(getCppInput(CppLib['doubledInts'].arguments,[cols]));
  worker.onmessage = function(e) {
    let output = getResult(CppLib['doubledInts'], e.data);

    // The following is added manually
    output.name = 'Doubled_ints';
    grok.shell.addTableView(output);
  }
}

//name: doubledFloats
//input: dataframe table
//input: column_list cols
//output: dataframe result 
export function doubledFloats(table, cols) {
  // this line was generated automatically
  //return callWasm(CppLib, 'doubledFloats', [cols]);

  // The following is added manually
  let res = callWasm(CppLib, 'doubledFloats', [cols]);
  res.name = 'Doubled_floats';
  return res;
}

//name: doubledFloats
//input: dataframe table
//input: column_list cols
export function doubledFloatsInWebWorker(table, cols) {
  var worker = new Worker(new URL('../wasm/doubledFloatsWorker.js', import.meta.url));
  worker.postMessage(getCppInput(CppLib['doubledFloats'].arguments,[cols]));
  worker.onmessage = function(e) {
    let output = getResult(CppLib['doubledFloats'], e.data);

    // The following is added manually
    output.name = 'Doubled_floats';
    grok.shell.addTableView(output);
  }
}

//name: pca
//input: dataframe table
//input: column_list columns
//input: int componentsCount
//output: dataframe result 
export function pca(table, columns, componentsCount) {
  return callWasm(CppLib, 'principalComponentAnalysis', [columns, componentsCount]);
}

//name: pca
//input: dataframe table
//input: column_list columns
//input: int componentsCount
export function pcaInWebWorker(table, columns, componentsCount) {
  var worker = new Worker(new URL('../wasm/principalComponentAnalysisWorker.js', import.meta.url));
  worker.postMessage(getCppInput(CppLib['principalComponentAnalysis'].arguments,[columns, componentsCount]));
  worker.onmessage = function(e) {
    let output = getResult(CppLib['principalComponentAnalysis'], e.data);

    // The following is added manually
    output.name = 'Principal components';
    grok.shell.addTableView(output);
  }
}

//name: error
//input: dataframe df
//input: column col1
//input: column col2
//output: double mad 
export function error(df, col1, col2) {
  return callWasm(CppLib, 'error', [col1, col2]);
}

//name: error
//input: dataframe df
//input: column col1
//input: column col2
export function errorInWebWorker(df, col1, col2) {
  var worker = new Worker(new URL('../wasm/errorWorker.js', import.meta.url));
  worker.postMessage(getCppInput(CppLib['error'].arguments,[col1, col2]));
  worker.onmessage = function(e) {
    let output = getResult(CppLib['error'], e.data);

    // The following is added manually
    alert(`Deviation of ${col1.name} and ${col2.name} is ${output}`);
  }
}

//name: pls
//input: dataframe table
//input: column_list features
//input: column predict
//input: int componentsCount
export function pls(table, features, predict, componentsCount) {
  // this line was generated automatically
  //return callWasm(CppLib, 'partialLeastSquareRegression', [features, predict, componentsCount]);

  // ALL the following is added manually
  let callOutput = callWasm(CppLib, 'partialLeastSquareRegression', [features, predict, componentsCount]);

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
}

//name: pls
//input: dataframe table
//input: column_list features
//input: column predict
//input: int componentsCount
export function plsInWebWorker(table, features, predict, componentsCount) {
  var worker = new Worker(new URL('../wasm/partialLeastSquareRegressionWorker.js', import.meta.url));
  worker.postMessage(getCppInput(CppLib['partialLeastSquareRegression'].arguments,[features, predict, componentsCount]));
  worker.onmessage = function(e) {
    let output = getResult(CppLib['partialLeastSquareRegression'], e.data);

    // ALL the following is added manually
    let dfView = grok.shell.getTableView(table.name);

    // 1. Predicted vs Reference scatter plot
  
    let prediction = output[0];
    prediction.name = predict.name + '(predicted)';
  
    let dfReferencePrediction = DG.DataFrame.fromColumns([predict, prediction]);
    dfReferencePrediction.name = 'Reference vs. Predicted';

    //grok.shell.addTableView(dfReferencePrediction);
    
    dfView.addViewer(DG.Viewer.scatterPlot(dfReferencePrediction, 
      { title: dfReferencePrediction.name,
        x: predict.name,
        y: prediction.name,
        showRegressionLine: true,
        markerType: 'circle'
       }));

    // 2. Regression Coefficients Bar Chart
    let regressionCoefficients = output[1];
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
  
    let xScores = output[2];
    for(let i = 0; i < xScores.length; i++) {
      xScores[i].name = `x.score.t${i+1}`;
      scoresColumns.push(xScores[i]);
    }

    let yScores = output[3];
    for(let i = 0; i < yScores.length; i++) {
      yScores[i].name = `y.score.u${i+1}`;
      scoresColumns.push(yScores[i]);
    }  

    let scores = DG.DataFrame.fromColumns(scoresColumns);
    scores.name = 'Scores';    

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

    let xLoadings = output[4];
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
  }
}

