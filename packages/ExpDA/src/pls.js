// pls.js

// Partial least-squares regression (PLS)

// Imports for call wasm runtime-system: in the main stream and in webworkers
import { callWasm } from '../wasm/callWasm';
import { getCppInput, getResult } from '../wasm/callWasmForWebWorker';

// Inputs checkers
import {checkComponenets, isProcExpensive} from './utils';

const PLS_MAX = 10000000;
const PLS_METHOD = 'partialLeastSquareRegression';

// Add results of computations: UI is provided here
function addPLSresults(table, features, predict, callOutput) {

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
} // addPLSresults

// Get results of PLS-computations
function getPLS(table, features, predict, components) { 
  addPLSresults(table, features, predict, callWasm(EDALib, PLS_METHOD, [features, predict, components]));
}

// Get results of PLS-computations: in webworker computations
function getPLScomputedInWebWorker(table, features, predict, components) {
  var worker = new Worker(new URL('../wasm/partialLeastSquareRegressionWorker.js', import.meta.url));
  
  worker.postMessage(getCppInput(EDALib[PLS_METHOD].arguments,
    [features, predict, components]));

  worker.onmessage = function(e) {
    addPLSresults(table, features, predict, getResult(EDALib[PLS_METHOD], e.data));    
  }
}

// Perform PLS-analysis
export function performPLS(table, features, predict, components) {
  checkComponenets(features, components);
  
  if (isProcExpensive(table, features, PLS_MAX))
    ui.dialog({title:'Warning'})
      .add(ui.span(['Computation may require much time. Continue it?']))
      .onCancel(() => { return;})
      .onOK(() => { getPLScomputedInWebWorker(table, features, predict, components); })
      .show();
  else
    getPLS(table, features, predict, components);
}