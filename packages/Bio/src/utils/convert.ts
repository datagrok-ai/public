import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {NotationConverter, NOTATION} from './notation-converter';

/**
 * Converts notations of a Macromolecule column
 *
 * @param {DG.column} col Column with 'Macromolecule' semantic type
 */
export function convert(col: DG.Column): void {
  const converter = new NotationConverter(col);
  const current: NOTATION = converter.sourceNotation;
  //TODO: read all notations
  const notations = [
    NOTATION.FASTA,
    NOTATION.SEPARATOR,
    NOTATION.HELM
  ];
  const separatorArray = ['-', '.', '/'];
  const filteredNotations = notations.filter((e) => e !== current);
  const targetNotationInput = ui.choiceInput('Convert to', filteredNotations[0], filteredNotations);

  const separatorInput = ui.choiceInput('Choose separator', separatorArray[0], separatorArray);

  ui.dialog('Convert sequence notation')
    .add(ui.div([
      ui.h1('Current notation: ' + current),
      targetNotationInput.root,
      // TODO: conditional separator input
      separatorInput.root
    ]))
    .onOK(() => {
      //TODO: create new converted column
      const targetNotation = targetNotationInput.value as NOTATION;
      const separator = separatorInput.value!;
      const newColumn = converter.convert(targetNotation, separator);
      col.dataFrame.columns.add(newColumn);
    })
    .show();
}
