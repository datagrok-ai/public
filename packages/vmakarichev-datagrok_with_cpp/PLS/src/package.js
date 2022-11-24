/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

import { callWasm } from '../wasm/callWasm';

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export async function init() {
  await initEigenPLS();  
}

//name: pls
//input: dataframe df
//input: column_list predictorColumns
//input: column responceColumn
//input: int componentsCount = 2
export function pls(df, predictorColumns, responceColumn, componentsCount) {
  
  let prediction = callWasm(EigenPLS, 'partialLeastSquareRegression', 
    [predictorColumns, responceColumn, componentsCount]);

  prediction.name = responceColumn.name + '(pred)';

  df.columns.add(prediction);  
  
}

