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
      const resColumn = substructureSearch(substructure, col);
      if (!col.dataFrame.columns.names().includes(resColumn.name))
        col.dataFrame.columns.add(resColumn);
      else
        grok.shell.warning(`Search ${substructure} is already performed`);
    })
    .show();
}

export function substructureSearch(substructure: string, col: DG.Column): DG.Column<boolean> {
  const lowerCaseSubstr = substructure.toLowerCase();
  const resultArray = new Array<boolean>(col.length);
  for (let i = 0; i < col.length; i++) {
    const macromolecule = col.get(i).toLowerCase();
    resultArray[i] = macromolecule.indexOf(lowerCaseSubstr) !== -1 ? true : false;
  }
  return DG.Column.fromList('bool', `Matches: ${substructure}`, resultArray);
}
