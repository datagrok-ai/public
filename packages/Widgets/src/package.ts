/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import {RadioButtonFilter} from './filters/radio-button-filter';
import {MultiValueFilter} from './filters/multi-value-filter';
import {SmilesLengthWidget, TimeWidget} from './widgets';

export const _package = new DG.Package();

//name: Single Choice
//description: A filter that lets you select exactly one category
//tags: filter
//output: filter result
export function radioButtonFilter() {
  return new RadioButtonFilter();
}

//name: Multi Choice
//description: A filter that works with columns of multi-value cells (such as lists of identifiers)
//tags: filter
//output: filter result
export function multiValueFilter() {
  return new MultiValueFilter();
}

//name: TimeWidget
//description: Shows current time
//output: widget result
export function timeWidget() {
  return new TimeWidget();
}

//name: SmilesLengthWidgetPanel
//input: string smiles = CN1C=NC2=C1C(=O)N(C(=O)N2C)C {semType: Molecule}
//tags: panel
//output: widget result
export function smilesLengthWidgetPanel(smiles: string) {
  return new SmilesLengthWidget().apply({smiles: smiles});
}

//name: SmilesLengthWidget
//output: widget result
export function smilesLengthWidget() {
  return new SmilesLengthWidget();
}
