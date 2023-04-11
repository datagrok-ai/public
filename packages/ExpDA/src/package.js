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
  await initEDALib();
}

//top-menu: Tools | Data Science | PCA
//name: PCA
//input: dataframe Table
//input: column_list features
//input: int components = 2
//output: dataframe PCA {action:join(Table)}
export function PCA(Table, features, components) {
  // check components count
  if(components <= 0) {
    let bal = new DG.Balloon();
    bal.error('components must be positive');
    return;
  }

  // call wasm computation of PCA
  let pca = callWasm(EDALib, 'principalComponentAnalysis', [features, components]);  

  // rename PCA-columns
  for(const col of pca.columns.toList())
    col.name = 'PCA' + col.name;

  return pca;
}

//name: PCA
//input: dataframe Table
//input: column_list features
//input: int components = 2
export function PCAInWebWorker(Table, features, components) {
  // check components count
  if(components <= 0) {
    let bal = new DG.Balloon();
    bal.error('components must be positive');
    return;
  }

  var worker = new Worker(new URL('../wasm/principalComponentAnalysisWorker.js', import.meta.url));
  worker.postMessage(getCppInput(EDALib['principalComponentAnalysis'].arguments,
    [features, components]));

  worker.onmessage = function(e) {
    let pca = getResult(EDALib['principalComponentAnalysis'], e.data);

    // rename PCA-columns and add them to Table
    for(const col of pca.columns.toList()) {
      col.name = 'PCA' + col.name;
      Table.columns.add(col);
    }

    /*output.name = 'Principal components';
    grok.shell.addTableView(output);*/
  }
}

//name: error
//input: dataframe df
//input: column col1
//input: column col2
//output: double mad 
export function error(df, col1, col2) {
  return callWasm(EDALib, 'error', [col1, col2]);
}

//name: error
//input: dataframe df
//input: column col1
//input: column col2
export function errorInWebWorker(df, col1, col2) {
  var worker = new Worker(new URL('../wasm/errorWorker.js', import.meta.url));
  worker.postMessage(getCppInput(EDALib['error'].arguments,[col1, col2]));
  worker.onmessage = function(e) {
    let output = getResult(EDALib['error'], e.data);

    // Provide output usage!
    alert(`Deviation of ${col1.name} and ${col2.name} is ${output}`);
  }
}

//top-menu: Tools | Data Science | PLS
//name: pls
//input: dataframe table
//input: column_list features
//input: column predict
//input: int components = 3
export function pls(table, features, predict, components) {
  // check components count
  if(components <= 0) {
    let bal = new DG.Balloon();
    bal.error('components must be positive');
    return;
  }

   // ALL the following is added manually
   let callOutput = callWasm(EDALib, 'partialLeastSquareRegression', 
     [features, predict, components]);

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
//input: int components = 3
export function plsInWebWorker(table, features, predict, components) {
  // check components count
  if(components <= 0) {
    let bal = new DG.Balloon();
    bal.error('components must be positive');
    return;
  }

  var worker = new Worker(new URL('../wasm/partialLeastSquareRegressionWorker.js', import.meta.url));
  worker.postMessage(getCppInput(EDALib['partialLeastSquareRegression'].arguments,
    [features, predict, components]));

  worker.onmessage = function(e) {
    let output = getResult(EDALib['partialLeastSquareRegression'], e.data);

    // Provide output usage!

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

