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
  const filtered = notations.filter((e) => e !== current);
  const choices = ui.choiceInput('Convert to', filtered[0], filtered);

  ui.dialog('Convert sequence')
    .add(
      ui.div([
        ui.h1('current notation'),
        ui.div(current),
        choices.root
      ])
    )
    .onOK(() => {
      //TODO: create new converted column
      converter.targetNotation = choices.value!;
      const newColumn = converter.convert();
      col.dataFrame.columns.add(newColumn);
    })
    .show();
}
