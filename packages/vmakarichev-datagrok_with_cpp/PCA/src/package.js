/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

// Tools for calling exported C/C++-functions
import {cppFuncWrapper, Arg, ArgColumn, ArgNewColumn, ArgColumns, ArgNewColumns} 
       from '../wasm/cppRuntimeSystem';

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export async function init() {
  await initEigenPCA();  
}

//top-menu: Tools | Data Science | Principal Component Analysis by Eigen
//name: pca
//tags: ml
//input: dataframe table
//input: column_list columns
//input: int numOfPrincipalComponents
//output: dataframe principalComponents
export function pca(table, columns, numOfPrincipalComponents) {

  // create arguments for exported C/C++-function call
  let data = columns.toList();
  let argData = new ArgColumns(data, 'f32');  
  let argNumComp = new Arg(numOfPrincipalComponents);
  let argPrincipalComponents = new ArgNewColumns('f32', data[0].length, numOfPrincipalComponents);    
  let args = [argData, argNumComp, argPrincipalComponents];
  
  // call exported C/C++-function: result code is obtained
  // let resultCode = cppFuncWrapper(EigenPCA, 'principalComponentAnalysis', 'num', args);  
  // console.log('Exported C/C++-function result code: ' + resultCode);

  // call exported function
  cppFuncWrapper(EigenPCA, 'principalComponentAnalysis', 'num', args);

  return DG.DataFrame.fromColumns(argPrincipalComponents.data); 
}
