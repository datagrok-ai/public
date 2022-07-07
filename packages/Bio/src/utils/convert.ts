import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {NotationConverter} from './notation-converter';

/**
 * Converts notations of a Macromolecule column
 *
 * @param {DG.column} col Column with 'Macromolecule' semantic type
 */
export function convert(col: DG.Column): void {
  const current = col.tags[DG.TAGS.UNITS];
  //TODO: read all notations
  const units = [
    'fasta',
    'separator',
    'HELM'
  ];
  const choices = ui.choiceInput('convert to', '', units.filter((e) => e !== current));

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
      const converter = new NotationConverter(col, choices.value!);
      const newColumn = converter.convert();
      col.dataFrame.columns.add(newColumn);
    })
    .show();
}
