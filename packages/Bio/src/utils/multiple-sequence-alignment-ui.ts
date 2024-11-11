import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {ColumnInputOptions} from '@datagrok-libraries/utils/src/type-declarations';
import {delay} from '@datagrok-libraries/utils/src/test';
import {ALPHABET, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {MsaWarning, runKalign} from './multiple-sequence-alignment';
import {pepseaMethods, runPepsea} from './pepsea';
import {checkInputColumn} from './check-input-column';
import {MultipleSequenceAlignmentUIOptions} from './types';
import {kalignVersion, msaDefaultOptions} from './constants';

import '../../css/msa.css';
import {_package} from '../package';

export async function multipleSequenceAlignmentUI(
  options: MultipleSequenceAlignmentUIOptions, seqHelper: ISeqHelper,
): Promise<DG.Column> {
  return new Promise(async (resolve, reject) => {
    options.clustersCol ??= null;
    options.pepsea ??= {};
    options.pepsea.method ??= msaDefaultOptions.pepsea.method;
    options.pepsea.gapOpen ??= msaDefaultOptions.pepsea.gapOpen;
    options.pepsea.gapExtend ??= msaDefaultOptions.pepsea.gapExtend;

    const table = options.col?.dataFrame ?? grok.shell.t;
    const seqCol = options.col ?? table.columns.bySemType(DG.SEMTYPE.MACROMOLECULE);
    if (seqCol == null) {
      const errMsg: string = `Multiple Sequence Alignment analysis requires a dataset with a macromolecule column.`;
      grok.shell.warning(errMsg);
      reject(new MsaWarning(ui.divText(errMsg)));
      return; // Prevents creating the MSA dialog
    }

    // UI for PepSea alignment
    const methodInput = ui.input.choice('Method', {value: options.pepsea.method, items: pepseaMethods});
    methodInput.setTooltip('Alignment method');

    // UI for Kalign alignment
    const terminalGapInput = ui.input.float('Terminal gap', {value: options?.kalign?.terminalGap});
    terminalGapInput.setTooltip('Penalty for opening a gap at the beginning or end of the sequence');
    const kalignVersionDiv = ui.p(`Kalign version: ${kalignVersion}`, 'kalign-version');

    // shared UI
    const gapOpenInput = ui.input.float('Gap open', {value: options.pepsea.gapOpen});
    gapOpenInput.setTooltip('Gap opening penalty at group-to-group alignment');
    const gapExtendInput = ui.input.float('Gap extend', {value: options.pepsea.gapExtend});
    gapExtendInput.setTooltip('Gap extension penalty to skip the alignment');

    const msaParamsDiv = ui.inputs([gapOpenInput, gapExtendInput, terminalGapInput]);
    const msaParamsButton = ui.button('Alignment parameters', () => {
      msaParamsDiv.hidden = !msaParamsDiv.hidden;
      [gapOpenInput, gapExtendInput, terminalGapInput].forEach((input) => {
        input.root.style.removeProperty('max-width');
        input.captionLabel.style.removeProperty('max-width');
      });
    }, 'Adjust alignment parameters such as penalties for opening and extending gaps');
    msaParamsButton.classList.add('msa-params-button');
    msaParamsDiv.hidden = true;
    msaParamsButton.prepend(ui.icons.settings(() => null));
    const pepseaInputRootStyles: CSSStyleDeclaration[] = [methodInput.root.style];
    const kalignInputRootStyles: CSSStyleDeclaration[] = [terminalGapInput.root.style, kalignVersionDiv.style];

    let performAlignment: (() => Promise<DG.Column<string> | null>) | undefined;

    let prevSeqCol = seqCol;
    const colInput = ui.input.column(
      'Sequence', {
        table: table, value: seqCol, onValueChanged: async (value: DG.Column<any>) => {
          if (!value || value.semType !== DG.SEMTYPE.MACROMOLECULE) {
            okBtn.disabled = true;
            await delay(0); // to
            colInput.value = prevSeqCol as DG.Column<string>;
            return;
          }
          prevSeqCol = value;
          okBtn.disabled = false;
          performAlignment = await onColInputChange(
            colInput.value, table, seqHelper, pepseaInputRootStyles, kalignInputRootStyles,
            methodInput, clustersColInput, gapOpenInput, gapExtendInput, terminalGapInput,
          );
        }, filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE
      } as ColumnInputOptions
    ) as DG.InputBase<DG.Column<string>>;
    colInput.setTooltip('Sequences column to use for alignment');
    const clustersColInput = ui.input.column('Clusters', {table: table, value: options.clustersCol!});
    clustersColInput.nullable = true;

    const dlg = ui.dialog('MSA')
      .add(colInput)
      .add(clustersColInput)
      .add(methodInput)
      .add(msaParamsDiv)
      .add(msaParamsButton)
      .add(kalignVersionDiv)
      .onOK(async () => { await onDialogOk(colInput, table, performAlignment, resolve, reject); });
    const okBtn = dlg.getButton('OK');

    colInput.fireChanged(); // changes okBtn
    //if column is specified (from tests), run alignment and resolve with the result
    if (options.col) {
      performAlignment = await onColInputChange(
        options.col, table, seqHelper, pepseaInputRootStyles, kalignInputRootStyles,
        methodInput, clustersColInput, gapOpenInput, gapExtendInput, terminalGapInput,
      );
      await onDialogOk(colInput, table, performAlignment, resolve, reject);
      return; // Prevents show the dialog
    }

    dlg.show();
  });
}

async function onDialogOk(
  colInput: DG.InputBase<DG.Column<any>>,
  table: DG.DataFrame,
  performAlignment: (() => Promise<DG.Column<string> | null>) | undefined,
  resolve: (value: DG.Column<any>) => void,
  reject: (reason: any) => void,
): Promise<void> {
  let msaCol: DG.Column<string> | null = null;
  const pi = DG.TaskBarProgressIndicator.create('Analyze for MSA ...');
  try {
    colInput.fireChanged();
    if (colInput.value.semType !== DG.SEMTYPE.MACROMOLECULE)
      throw new Error('Chosen column has to be of Macromolecule semantic type');
    if (performAlignment === undefined) // value can only be undefined when column can't be processed with either method
      throw new Error('Invalid column format');
    msaCol = await performAlignment(); // progress
    if (msaCol == null)
      return reject('PepSeA container has not started');

    table.columns.add(msaCol);
    await grok.data.detectSemanticTypes(table);

    resolve(msaCol);
  } catch (err: any) {
    reject(err);
  } finally {
    pi.close();
  }
}


async function onColInputChange(
  col: DG.Column<string>, table: DG.DataFrame, seqHelper: ISeqHelper,
  pepseaInputRootStyles: CSSStyleDeclaration[], kalignInputRootStyles: CSSStyleDeclaration[],
  methodInput: DG.InputBase<string | null>, clustersColInput: DG.InputBase<DG.Column<any> | null>,
  gapOpenInput: DG.InputBase<number | null>, gapExtendInput: DG.InputBase<number | null>,
  terminalGapInput: DG.InputBase<number | null>,
): Promise<(() => Promise<DG.Column<string> | null>) | undefined> {
  try {
    if (col.semType !== DG.SEMTYPE.MACROMOLECULE)
      return;
    const unusedName = table.columns.getUnusedName(`msa(${col.name})`);

    if (checkInputColumn(col, col.name, seqHelper,
      [NOTATION.FASTA, NOTATION.SEPARATOR], [ALPHABET.DNA, ALPHABET.RNA, ALPHABET.PT])[0]
    ) { // Kalign - natural alphabets. if the notation is separator, convert to fasta and then run kalign
      switchDialog(pepseaInputRootStyles, kalignInputRootStyles, 'kalign');
      gapOpenInput.value = null;
      gapExtendInput.value = null;
      terminalGapInput.value = null;
      const potentialColSh = seqHelper.getSeqHandler(col);
      const performCol: DG.Column<string> = potentialColSh.isFasta() ? col :
        potentialColSh.convert(NOTATION.FASTA);
      return async () => await runKalign(performCol, false, unusedName, clustersColInput.value);
    } else if (checkInputColumn(col, col.name, seqHelper, [NOTATION.HELM], [])[0]) {
      // PepSeA branch - Helm notation or separator notation with unknown alphabets
      switchDialog(pepseaInputRootStyles, kalignInputRootStyles, 'pepsea');
      gapOpenInput.value ??= msaDefaultOptions.pepsea.gapOpen;
      gapExtendInput.value ??= msaDefaultOptions.pepsea.gapExtend;

      return async () => {
        return runPepsea(col, unusedName, methodInput.value!,
          gapOpenInput.value!, gapExtendInput.value!, clustersColInput.value);
      };
    } else if (checkInputColumn(col, col.name, seqHelper, [NOTATION.SEPARATOR], [ALPHABET.UN])[0]) {
      //if the column is separator with unknown alphabet, it might be helm. check if it can be converted to helm
      const potentialColSh = seqHelper.getSeqHandler(col);
      const helmCol = potentialColSh.convert(NOTATION.HELM);
      switchDialog(pepseaInputRootStyles, kalignInputRootStyles, 'pepsea');
      gapOpenInput.value ??= msaDefaultOptions.pepsea.gapOpen;
      gapExtendInput.value ??= msaDefaultOptions.pepsea.gapExtend;
      // convert to helm and assign alignment function to PepSea

      return async () => {
        return runPepsea(helmCol, unusedName, methodInput.value!,
          gapOpenInput.value!, gapExtendInput.value!, clustersColInput.value);
      };
    } else {
      gapOpenInput.value = null;
      gapExtendInput.value = null;
      terminalGapInput.value = null;
      switchDialog(pepseaInputRootStyles, kalignInputRootStyles, 'kalign');
      return;
    }
  } catch (err: any) {
    const errMsg: string = err instanceof Error ? err.message : err.toString();
    grok.shell.error(errMsg);
    _package.logger.error(errMsg);
  }
}

type MSADialogType = 'kalign' | 'pepsea';

function switchDialog(
  pepseaInputRootStyles: CSSStyleDeclaration[], kalignInputRootStyles: CSSStyleDeclaration[], dialogType: MSADialogType,
) {
  if (dialogType === 'kalign') {
    for (const inputRootStyle of pepseaInputRootStyles)
      inputRootStyle.display = 'none';
    for (const inputRootStyle of kalignInputRootStyles)
      inputRootStyle.removeProperty('display');
  } else {
    for (const inputRootStyle of kalignInputRootStyles)
      inputRootStyle.display = 'none';
    for (const inputRootStyle of pepseaInputRootStyles)
      inputRootStyle.removeProperty('display');
  }
}
