// Exploratory data analysis (EDA) tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_principalComponentAnalysisInWebWorker,
  _partialLeastSquareRegressionInWebWorker} from '../wasm/EDAAPI';

import {checkComponenets} from './utils';

const PCA_COL_NAME = 'PCA';

// Principal components analysis (PCA)
export async function computePCA(table: DG.DataFrame, features: DG.ColumnList, components: number): Promise<DG.DataFrame> 
{
  checkComponenets(features, components);

  let _output: any;
  let _promise = _principalComponentAnalysisInWebWorker(table, features, components);

  await _promise.then(
    _result => {  
       _output = _result;

       // rename columns
       for (const col of _output.columns.toList()) 
         col.name = PCA_COL_NAME + col.name;              
    },
    _error => {  throw new Error (`Error: ${_error}`); }
  );

  return _output;  
} 

// Add results of computations: custom UI is provided here
function addPLSresults(table: DG.DataFrame, features: DG.ColumnList, predict: DG.Column, callOutput: any) {

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
  for (const col of features)
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
  for (let i = 0; i < xScores.length; i++) {
    xScores[i].name = `x.score.t${i+1}`;
    scoresColumns.push(xScores[i]);
  }
 
  let yScores = callOutput[3];
  for (let i = 0; i < yScores.length; i++) {
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
  for (let col of features)
    loadingLabels.push(col.name);
 
  loadingCols.push(DG.Column.fromStrings('labels', loadingLabels));
 
  let xLoadings = callOutput[4];
  for (let i = 0; i < xLoadings.length; i++) {
    xLoadings[i].name = `x.loading.p${i+1}`;
    loadingCols.push(xLoadings[i]);
  }
 
  let dfLoadings = DG.DataFrame.fromColumns(loadingCols);
  dfLoadings.name = 'Loadings';
   
  dfView.addViewer(DG.Viewer.scatterPlot(dfLoadings, 
    { title: dfLoadings.name,
      x: xLoadings[0].name,
      y: xLoadings[xLoadings.length - 1].name,      
      markerType: 'circle',
      labels: 'labels'
     }));
} // addPLSresults

// Partial least square regression (PLS)
export async function computePLS(table: DG.DataFrame, features: DG.ColumnList, predict: DG.Column, components: number): Promise<void> 
{
  checkComponenets(features, components);

  let _promise = _partialLeastSquareRegressionInWebWorker(table, features, predict, components);

  await _promise.then(
    _result => addPLSresults(table, features, predict, _result),    
    _error => {  throw new Error (`Error: ${_error}`); }
  );
} 