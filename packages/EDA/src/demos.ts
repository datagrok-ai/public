/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {computePLS} from './EDAtools';
import {addPLSvisualization} from './EDAui';

// Demo multivariate analysis (PLS)
export async function demoPLS(rowCount: number, colCount: number, componentsCount: number): Promise<void> {
  // check inputs
  if ((rowCount <= 0) || (colCount <= 0) || (componentsCount <= 0) || (componentsCount > colCount)) {
    const bal = new DG.Balloon;
    bal.error('Incorrect inputs.');
    return;
  }

  // further, custom interface is provided
  
  const PREDICT = 'Reference';
  
  const bigDemoTable = grok.data.testData('random walk', rowCount, colCount);   
  bigDemoTable.name = `${rowCount} x ${colCount}`;
  
  for (const col of bigDemoTable.columns)
    col.name = 'Feature ' + col.name;
  bigDemoTable.columns.byIndex(0).name = PREDICT;
  
  grok.shell.addTableView(bigDemoTable);
  let predict = bigDemoTable.columns.byName(PREDICT);
  let features = bigDemoTable.columns.remove(PREDICT);      
  
  const plsResults = await computePLS(bigDemoTable, features, predict, componentsCount);
  
  addPLSvisualization(bigDemoTable, features, predict, plsResults);
  
  bigDemoTable.columns.add(predict);
}
