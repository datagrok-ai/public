import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import * as C from '../utils/constants';
import {getMonomericMols} from '../calculations/monomerLevelMols';
import {BitSet} from 'datagrok-api/dg';
import {updateDivInnerHTML} from '../utils/ui-utils';

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
  df.col('substr_helm')!.setTag(DG.TAGS.UNITS, bio.NOTATION.HELM);
  const grid = df.plot.grid();
  const separatorInput = ui.textInput('Separator', separator);

  const inputsDiv = ui.div();

  const inputs = units === bio.NOTATION.HELM ? ui.divV([editHelmLink]) :
    units === bio.NOTATION.SEPARATOR ? ui.inputs([substructureInput, separatorInput]) :
      ui.inputs([substructureInput]);

  updateDivInnerHTML(inputsDiv, inputs);

  ui.dialog('Substructure search')
    .add(ui.divV([
      ui.divText(`Notation: ${units}`),
      inputsDiv
    ]))
    .onOK(async () => {
      let substructure = units === bio.NOTATION.HELM ? df.get('substr_helm', 0) : substructureInput.value;
      if (units === bio.NOTATION.SEPARATOR && separatorInput.value !== separator && separatorInput.value !== '')
        substructure = substructure.replaceAll(separatorInput.value, separator);
      const matchesColName = `Matches: ${substructure}`;
      const colExists = col.dataFrame.columns.names()
        .filter((it) => it.toLocaleLowerCase() === matchesColName.toLocaleLowerCase()).length > 0;
      if (!colExists) {
        let matches: BitSet;
        if (units === bio.NOTATION.HELM)
          matches = await helmSubstructureSearch(substructure, col);
        else
          matches = linearSubstructureSearch(substructure, col);
        col.dataFrame.columns.add(DG.Column.fromBitSet(matchesColName, matches));
      } else { grok.shell.warning(`Search ${substructure} is already performed`); }
    })
    .show();
}

export function linearSubstructureSearch(substructure: string, col: DG.Column): DG.BitSet {
  const lowerCaseSubstr = substructure.toLowerCase();
  const resultArray = DG.BitSet.create(col.length);
  for (let i = 0; i < col.length; i++) {
    const macromolecule = col.get(i).toLowerCase();
    if (macromolecule.indexOf(lowerCaseSubstr) !== -1)
      resultArray.set(i, true, false);
  }
  return resultArray;
}

async function helmSubstructureSearch(substructure: string, col: DG.Column): Promise<BitSet> {
  const helmColWithSubstructure = DG.Column.string('helm', col.length + 1)
    .init((i) => i === col.length ? substructure : col.get(i));
  helmColWithSubstructure.setTag(DG.TAGS.UNITS, bio.NOTATION.HELM);
  const monomericMolsCol = await getMonomericMols(helmColWithSubstructure, true);
  const molSubstructure = monomericMolsCol.get(col.length);
  const monomericMolsDf = DG.DataFrame.fromColumns([monomericMolsCol]);
  monomericMolsDf.rows.removeAt(col.length);
  const matchesCol = await grok.functions.call('Chem:searchSubstructure', {
    molStringsColumn: monomericMolsDf.columns.byIndex(0),
    molString: molSubstructure,
    molBlockFailover: '',
  });
  return matchesCol.get(0);
}
