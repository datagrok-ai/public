// pca.js

// Principal components analysis (PCA)

// Imports for call wasm runtime-system: in the main stream and in webworkers
import { callWasm } from '../wasm/callWasm';
import { getCppInput, getResult } from '../wasm/callWasmForWebWorker';

// Inputs checkers
import {checkComponenets, isProcExpensive} from './utils';

const PCA_MAX = 200000;
const PCA_COL_NAME = 'PCA';
const PCA_METHOD = 'principalComponentAnalysis';

// Add principal components to the table
function addPCAtoTable(table, pca) {  
  for (const col of pca.columns.toList()) {
    col.name = PCA_COL_NAME + col.name;
    table.columns.add(col);
  }
}

// Get principal components
export function getPCA(table, features, components) {
  addPCAtoTable(table, callWasm(EDALib, PCA_METHOD, [features, components]));
}

// Get principal components: in webworker computations
export function getPCAcomputedInWebWorker(table, features, components) {  
  var worker = new Worker(new URL('../wasm/principalComponentAnalysisWorker.js', import.meta.url));
 
  worker.postMessage(getCppInput(EDALib[PCA_METHOD].arguments,
    [features, components]));
  
  worker.onmessage = function(e) {    
    addPCAtoTable(table, getResult(EDALib[PCA_METHOD], e.data));
  }
}

// // Perform PCA-computation
export function performPCA(table, features, components) {
  checkComponenets(features, components);

  if (isProcExpensive(table, features, PCA_MAX))
    ui.dialog({title:'Warning'})
      .add(ui.span(['Computation may require much time. Continue it?']))
      .onCancel(() => { return;})
      .onOK(() => { getPCAcomputedInWebWorker(table, features, components); })
      .show();
  else
    getPCA(table, features, components); 
}