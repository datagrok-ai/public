import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import $ from 'cash-dom';
import {Subscription} from 'rxjs';

import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';


let convertDialog: DG.Dialog | null = null;
let convertDialogSubs: Subscription[] = [];

/**
 * Converts notations of a Macromolecule column
 *
 * @param {DG.column} col Column with 'Macromolecule' semantic type
 */
export function convert(col?: DG.Column): void {
  let srcCol = col ?? grok.shell.t.columns.bySemType('Macromolecule')!;
  if (!srcCol)
    throw new Error('No column with Macromolecule semantic type found');
  let converterSh = SeqHandler.forColumn(srcCol);
  let currentNotation: NOTATION = converterSh.notation;
  const dialogHeader = ui.divText(
    'Current notation: ' + currentNotation,
    {
      style: {
        'text-align': 'center',
        'font-weight': 'bold',
        'font-size': '14px',
        'padding': '5px',
      },
    },
  );
  const notations = [
    NOTATION.FASTA,
    NOTATION.SEPARATOR,
    NOTATION.HELM,
  ];
  const toggleColumn = (newCol: DG.Column) => {
    if (newCol.semType !== DG.SEMTYPE.MACROMOLECULE) {
      targetColumnInput.value = srcCol;
      return;
    }

    srcCol = newCol;
    converterSh = SeqHandler.forColumn(srcCol);
    currentNotation = converterSh.notation;
    if (currentNotation === NOTATION.HELM)
      separatorInput.value = '/'; // helm monomers can have - in the name like D-aThr;
    dialogHeader.textContent = 'Current notation: ' + currentNotation;
    filteredNotations = notations.filter((e) => e !== currentNotation);
    targetNotationInput = ui.input.choice('Convert to', {value: filteredNotations[0], items: filteredNotations,
      onValueChanged: toggleSeparator});
    toggleSeparator();
    convertDialog?.clear();
    convertDialog?.add(ui.div([
      dialogHeader,
      targetColumnInput.root,
      targetNotationInput.root,
      separatorInput.root
    ]));
  };

  const targetColumnInput = ui.input.column('Column', {table: grok.shell.t, value: srcCol,
    onValueChanged: (value) => toggleColumn(value)});

  const separatorArray = ['-', '.', '/'];
  let filteredNotations = notations.filter((e) => e !== currentNotation);

  const separatorInput = ui.input.choice('Separator', {value: separatorArray[0], items: separatorArray});

  // hide the separator input for non-SEPARATOR target notations
  const toggleSeparator = () => {
    if (targetNotationInput.value !== NOTATION.SEPARATOR)
      $(separatorInput.root).hide();
    else
      $(separatorInput.root).show();
  };
  let targetNotationInput = ui.input.choice('Convert to', {value: filteredNotations[0], items: filteredNotations,
    onValueChanged: toggleSeparator});

  // set correct visibility on init
  toggleSeparator();

  targetNotationInput.onChanged.subscribe(() => {
    toggleSeparator();
  });

  if (convertDialog == null) {
    convertDialog = ui.dialog('Convert Sequence Notation')
      .add(ui.div([
        dialogHeader,
        targetColumnInput.root,
        targetNotationInput.root,
        separatorInput.root,
      ]))
      .onOK(async () => {
        const targetNotation = targetNotationInput.value as NOTATION;
        const separator: string | undefined = targetNotation === NOTATION.SEPARATOR ? separatorInput.value! : undefined;

        await convertDo(srcCol, targetNotation, separator);
      })
      .show({x: 350, y: 100});

    convertDialogSubs.push(convertDialog.onClose.subscribe((_: any) => {
      convertDialogSubs.forEach((s) => { s.unsubscribe(); });
      convertDialogSubs = [];
      convertDialog = null;
    }));
  }
}

/** Creates a new column with converted sequences and detects its semantic type
 * @param {DG.Column} srcCol Column with 'Macromolecule' semantic type
 * @param {NOTATION} targetNotation Target notation
 * @param {string | null} separator Separator for SEPARATOR notation
 */
export async function convertDo(srcCol: DG.Column, targetNotation: NOTATION, separator?: string): Promise<DG.Column> {
  const converterSh = SeqHandler.forColumn(srcCol);
  const newColumn = converterSh.convert(targetNotation, separator);
  srcCol.dataFrame.columns.add(newColumn);

  // Call detector directly to escape some error on detectSemanticTypes
  const semType = await grok.functions.call('Bio:detectMacromolecule', {col: newColumn});
  if (semType)
    newColumn.semType = semType;

  // call to calculate 'cell.renderer' tag
  await grok.data.detectSemanticTypes(srcCol.dataFrame);

  return newColumn;
}
