import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {Subscription} from 'rxjs';
import {NotationConverter, NOTATION} from '@datagrok-libraries/bio/src/utils/notation-converter';


let convertDialog: DG.Dialog | null = null;
let convertDialogSubs: Subscription[] = [];

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

  if (convertDialog == null) {
    convertDialog = ui.dialog('Convert sequence notation')
      .add(ui.div([
        ui.h1('Current notation: ' + current),
        targetNotationInput.root,
        // TODO: conditional separator input
        separatorInput.root
      ]))
      .onOK(async () => {
        const targetNotation = targetNotationInput.value as NOTATION;
        const separator: string | null = separatorInput.value;

        await convertDo(col, targetNotation, separator);
      })
      .show();

    convertDialogSubs.push(convertDialog.onClose.subscribe((value) => {
      convertDialogSubs.forEach((s) => {s.unsubscribe(); });
      convertDialogSubs = [];
      convertDialog = null;
    }));
  }
}

export async function convertDo(srcCol: DG.Column, targetNotation: NOTATION, separator: string | null): Promise<DG.Column> {
  const converter = new NotationConverter(srcCol);
  const newColumn = converter.convert(targetNotation, separator);
  srcCol.dataFrame.columns.add(newColumn);
  await grok.data.detectSemanticTypes(srcCol.dataFrame);
  return newColumn;
}

