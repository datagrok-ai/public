import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as C from '../utils/constants';
import {getMonomericMols} from '../calculations/monomerLevelMols';
import {updateDivInnerHTML} from '../utils/ui-utils';
import {delay} from '@datagrok-libraries/utils/src/test';
import {NOTATION} from '@datagrok-libraries/bio';

export const MONOMER_MOLS_COL = 'monomeric-mols';

export const enum MONOMERIC_COL_TAGS {
  MONOMERIC_MOLS = 'monomeric-mols',
  LAST_INVALIDATED_VERSION = 'last-invalidated-version',
  MONOMERS_DICT = 'monomers-dict'
}

/**
 * Searches substructure in each row of Macromolecule column
 *
 * @param {DG.column} col Column with 'Macromolecule' semantic type
 */
export function substructureSearchDialog(col: DG.Column): void {
  const units = col.getTag(DG.TAGS.UNITS);
  const separator = col.getTag(C.TAGS.SEPARATOR);
  // const notations = [NOTATION.FASTA, NOTATION.SEPARATOR, NOTATION.HELM];

  const substructureInput = ui.textInput('Substructure', '');

  const editHelmLink = ui.link('Edit helm', async () => {
    updateDivInnerHTML(inputsDiv, grid.root);
    await ui.tools.waitForElementInDom(grid.root);
    setTimeout(() => {
      grid.cell('substr_helm', 0).element.children[0].dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter'}));
    }, 100);
  });

  const df = DG.DataFrame.create(1);
  df.columns.addNewString('substr_helm').init((i) => '');
  df.col('substr_helm')!.semType = col.semType;
  df.col('substr_helm')!.setTag(DG.TAGS.UNITS, NOTATION.HELM);
  const grid = df.plot.grid();
  const separatorInput = ui.textInput('Separator', separator);

  const inputsDiv = ui.div();

  const inputs = units === NOTATION.HELM ? ui.divV([editHelmLink]) :
    units === NOTATION.SEPARATOR ? ui.inputs([substructureInput, separatorInput]) :
      ui.inputs([substructureInput]);

  updateDivInnerHTML(inputsDiv, inputs);

  ui.dialog('Substructure search')
    .add(ui.divV([
      ui.divText(`Notation: ${units}`),
      inputsDiv
    ]))
    .onOK(async () => {
      let substructure = units === NOTATION.HELM ? df.get('substr_helm', 0) : substructureInput.value;
      if (units === NOTATION.SEPARATOR && separatorInput.value !== separator && separatorInput.value !== '')
        substructure = substructure.replaceAll(separatorInput.value, separator);
      const matchesColName = `Matches: ${substructure}`;
      const colExists = col.dataFrame.columns.names()
        .filter((it) => it.toLocaleLowerCase() === matchesColName.toLocaleLowerCase()).length > 0;
      if (!colExists) {
        let matches: DG.BitSet;
        if (units === NOTATION.HELM)
          matches = await helmSubstructureSearch(substructure, col);
        else
          matches = linearSubstructureSearch(substructure, col);
        col.dataFrame.columns.add(DG.Column.fromBitSet(matchesColName, matches));
      } else { grok.shell.warning(`Search ${substructure} is already performed`); }
    })
    .show();
}

export function linearSubstructureSearch(substructure: string, col: DG.Column, separator?: string): DG.BitSet {
  const re = separator ? prepareSubstructureRegex(substructure, separator) : substructure;
  const resultArray = DG.BitSet.create(col.length);
  for (let i = 0; i < col.length; i++) {
    const macromolecule = col.get(i);
    if (macromolecule.match(re) || macromolecule === substructure)
      resultArray.set(i, true, false);
  }
  return resultArray;
}

function prepareSubstructureRegex(substructure: string, separator: string) {
  const char = `${separator}`.replace(/[\-\[\]\/\{\}\(\)\*\+\?\.\\\^\$\|]/g, '\\$&');
  const startsWithSep = substructure.charAt(0) === separator;
  const endsWithSep = substructure.charAt(substructure.length - 1) === separator;
  const substrWithoutSep = substructure.replace(new RegExp(`^${char}|${char}$`, 'g'), '');
  const re = startsWithSep ? endsWithSep ? `${char}${substrWithoutSep}${char}` :
      `${char}${substrWithoutSep}${char}|${char}${substrWithoutSep}$` :
    endsWithSep ? `^${substrWithoutSep}${char}|${char}${substrWithoutSep}${char}` :
      `^${substrWithoutSep}${char}|${char}${substrWithoutSep}${char}|${char}${substrWithoutSep}$`;
  return re;
}

export async function helmSubstructureSearch(substructure: string, col: DG.Column): Promise<DG.BitSet> {
  if (col.version !== col.temp[MONOMERIC_COL_TAGS.LAST_INVALIDATED_VERSION])
    await invalidateMols(col, true);
  const substructureCol = DG.Column.string('helm', 1).init((i) => substructure);
  substructureCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
  const substructureMolsCol =
    await getMonomericMols(substructureCol, true, col.temp[MONOMERIC_COL_TAGS.MONOMERS_DICT]);
  const matchesCol = await grok.functions.call('Chem:searchSubstructure', {
    molStringsColumn: col.temp[MONOMERIC_COL_TAGS.MONOMERIC_MOLS],
    molString: substructureMolsCol.get(0),
    molBlockFailover: '',
  });
  return matchesCol.get(0);
}

export async function invalidateMols(col: DG.Column, pattern: boolean) {
  const progressBar = DG.TaskBarProgressIndicator.create(`Invalidating molfiles for ${col.name}`);
  await delay(10);
  const monomersDict = new Map();
  const monomericMolsCol = await getMonomericMols(col, pattern, monomersDict);
  col.temp[MONOMERIC_COL_TAGS.MONOMERIC_MOLS] = monomericMolsCol;
  col.temp[MONOMERIC_COL_TAGS.MONOMERS_DICT] = monomersDict;
  col.temp[MONOMERIC_COL_TAGS.LAST_INVALIDATED_VERSION] = col.version;
  progressBar.close();
}
