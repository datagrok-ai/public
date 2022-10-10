import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/units-handler';
import * as C from '../utils/constants';

/**
 * Searches substructure in each row of Macromolecule column
 *
 * @param {DG.column} col Column with 'Macromolecule' semantic type
 */
export function substructureSearchDialog(col: DG.Column): void {
  const units = col.getTag(DG.TAGS.UNITS);
  const separator = col.getTag(C.TAGS.SEPARATOR);
  const notations = [NOTATION.FASTA, NOTATION.SEPARATOR];

  const substructureInput = ui.textInput('Substructure', '');
  const notationInput = ui.choiceInput('Notation', units, notations);
  const separatorInput = ui.textInput('Separator', separator);

  // hide the separator input for non-SEPARATOR target notations
  const toggleSeparator = () => {
    if (notationInput.value !== NOTATION.SEPARATOR)
      separatorInput.root.hidden = true;
    else
      separatorInput.root.hidden = false;
  };

  toggleSeparator();

  notationInput.onChanged(() => {
    toggleSeparator();
  });

  ui.dialog('Substructure search')
    .add(ui.inputs([
      substructureInput,
      notationInput,
      separatorInput
    ]))
    .onOK(() => {
      let substructure = substructureInput.value;
      if (notationInput.value !== NOTATION.FASTA && separatorInput.value !== separator)
        substructure = substructure.replaceAll(separatorInput.value, separator);
      const matchesColName = `Matches: ${substructure}`;
      const colExists = col.dataFrame.columns.names()
        .filter((it) => it.toLocaleLowerCase() === matchesColName.toLocaleLowerCase()).length > 0;
      if (!colExists) {
        const matches = substructureSearch(substructure, col);
        col.dataFrame.columns.add(DG.Column.fromBitSet(matchesColName, matches));
      } else { grok.shell.warning(`Search ${substructure} is already performed`); }
    })
    .show();
}

export function substructureSearch(substructure: string, col: DG.Column): DG.BitSet {
  const lowerCaseSubstr = substructure.toLowerCase();
  const resultArray = DG.BitSet.create(col.length);
  for (let i = 0; i < col.length; i++) {
    const macromolecule = col.get(i).toLowerCase();
    if (macromolecule.indexOf(lowerCaseSubstr) !== -1)
      resultArray.set(i, true, false);
  }
  return resultArray;
}
