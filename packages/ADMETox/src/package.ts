/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { addPredictions } from './admet-analysis/admet-calculation';
import { getModelsSingle } from './admet-analysis/admet-calculation';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: ADME/Tox
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function admetWidget(smiles: string): DG.Widget<any> {
  return getModelsSingle(smiles);
}

//name: ADME/Tox...
//input: column col {semType: Molecule}
export async function admetoxCalculators(col: DG.Column) {
  let table = col.dataFrame;
  return await addPredictions(col, table);
}
