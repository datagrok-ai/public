/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import $ from 'cash-dom';
import {Subscription} from 'rxjs';

import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {NotationProviderBase} from '@datagrok-libraries/bio/src/utils/macromolecule/types';


let convertDialog: DG.Dialog | null = null;
let convertDialogSubs: Subscription[] = [];

/**
 * Converts notations of a Macromolecule column
 *
 * @param {DG.Column<string>} col Column with 'Macromolecule' semantic type
 * @param {ISeqHelper} seqHelper
 */
export async function convert(col: DG.Column<string> | undefined, seqHelper: ISeqHelper): Promise<void> {
  let srcCol = col ?? grok.shell.t.columns.bySemType('Macromolecule')!;
  if (!srcCol)
    throw new Error('No column with Macromolecule semantic type found');
  let converterSh = seqHelper.getSeqHandler(srcCol);
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
    NOTATION.BILN
  ];

  const notationProviderConstructors = await NotationProviderBase.getProviderConstructors();

  const toggleColumn = (newCol: DG.Column) => {
    srcCol = newCol;
    converterSh = seqHelper.getSeqHandler(srcCol);
    currentNotation = converterSh.notation;
    if (currentNotation === NOTATION.HELM)
      separatorInput.value = '/'; // helm monomers can have - in the name like D-aThr;
    dialogHeader.textContent = 'Current notation: ' + currentNotation;
    filteredNotations = notations.filter((e) => e !== currentNotation);
    if (currentNotation === NOTATION.CUSTOM)
      filteredNotations = [NOTATION.HELM];
    if (currentNotation === NOTATION.HELM) {
      // add custom notations that
      notationProviderConstructors.forEach((c) => {
        if (c.implementsFromHelm)
          filteredNotations.unshift(c.notationName as NOTATION); // hack :)
      });
    }
    targetNotationInput = ui.input.choice('Convert to', {
      value: filteredNotations[0], items: filteredNotations,
      onValueChanged: toggleSeparator
    });
    toggleSeparator();
    convertDialog?.clear();
    convertDialog?.add(ui.div([
      dialogHeader,
      targetColumnInput.root,
      targetNotationInput.root,
      separatorInput.root
    ]));
  };

  const targetColumnInput = ui.input.column('Column', {
    table: grok.shell.t, value: srcCol,
    filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE,
    nullable: false,
    onValueChanged: (value) => toggleColumn(value)
  });

  const separatorArray = ['-', '.', '/'];
  let filteredNotations = notations.filter((e) => e !== currentNotation);
  if (currentNotation === NOTATION.CUSTOM)
    filteredNotations = [NOTATION.HELM];
  if (currentNotation === NOTATION.HELM) {
    // add custom notations that
    notationProviderConstructors.forEach((c) => {
      if (c.implementsFromHelm)
        filteredNotations.unshift(c.notationName as NOTATION); // hack :)
    });
  }

  const separatorInput = ui.input.choice('Separator', {value: separatorArray[0], items: separatorArray});

  // hide the separator input for non-SEPARATOR target notations
  const toggleSeparator = () => {
    if (targetNotationInput.value !== NOTATION.SEPARATOR)
      $(separatorInput.root).hide();
    else
      $(separatorInput.root).show();
  };
  let targetNotationInput = ui.input.choice('Convert to', {
    value: filteredNotations[0], items: filteredNotations,
    onValueChanged: toggleSeparator
  });

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
        let notationProviderConstructor: typeof NotationProviderBase | undefined = undefined;
        if (!notations.includes(targetNotation) && notationProviderConstructors.find((c) => c.notationName === targetNotation))
          notationProviderConstructor = notationProviderConstructors.find((c) => c.notationName === targetNotation)!;

        await convertDo(srcCol, seqHelper, targetNotation, separator, notationProviderConstructor);
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
export async function convertDo(srcCol: DG.Column, seqHelper: ISeqHelper, targetNotation: NOTATION, separator?: string, notationProviderConstructor?: typeof NotationProviderBase): Promise<DG.Column> {
  if (notationProviderConstructor) {
    const newColName = srcCol.dataFrame.columns.getUnusedName(`${notationProviderConstructor.notationName}(${srcCol.name})`);
    const newCol = DG.Column.string(newColName, srcCol.length);
    newCol.init((i) => {
      const seq = srcCol.get(i);// we know for sure (in Macron accent) that it is helm
      try {
        return notationProviderConstructor.convertFromHelm(seq, {});
      } catch (e) {
        console.error(`Error converting sequence at row ${i}: ${e instanceof Error ? e.message : e}`);
        return '';
      }
    });
    srcCol.dataFrame.columns.add(newCol);
    const semType = await grok.functions.call('Bio:detectMacromolecule', {col: newCol});
    if (semType)
      newCol.semType = semType;
    await grok.data.detectSemanticTypes(srcCol.dataFrame);
    return newCol;
  } else {
    const converterSh = seqHelper.getSeqHandler(srcCol);
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
}
