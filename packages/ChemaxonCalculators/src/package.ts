/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { getCalculationsSingle } from './chemaxon-calculators/chemaxon-calculators';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: getLogP
export async function getLogP() {
  let data = await grok.data.query(`${_package.name}:logPCalculation`, {
    'atomIncrements': false,
    'inputFormat': 'smiles',
    'method': 'CHEMAXON',
    'structure': 'CCC'
  });
  grok.shell.addTableView(data);
}

//name: getHLB
export async function getHLB() {
  let data = await grok.data.query(`${_package.name}:HLBNumberCalculation`, {
    'inputFormat': 'smiles',
    'structure': 'CCC'
  });
  grok.shell.addTableView(data); 
}

//name: getCharge 
export async function getCharge() {
  let data = await grok.data.query(`${_package.name}:HLBNumberCalculation`, {
    'inputFormat': 'smiles',
    'ph': 7.4,
    'structure': 'CCC'
  });
  grok.shell.addTableView(data); 
}

//name: Chemaxon Calculators
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function admetWidget(smiles: string): DG.Widget {
  return getCalculationsSingle(smiles);
}


