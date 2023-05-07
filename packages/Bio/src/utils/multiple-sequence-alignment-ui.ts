import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {ALPHABET, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {runKalign} from './multiple-sequence-alignment';
import {pepseaMethods, runPepsea} from './pepsea';
import {checkInputColumnUI} from './check-input-column';
import {NotationConverter} from '@datagrok-libraries/bio/src/utils/notation-converter';
import {_package} from '../package';

export class MsaWarning extends Error {
  constructor(message: string, options?: ErrorOptions) {
    super(message, options);
  }
}

export async function multipleSequenceAlignmentUI(
  col: DG.Column<string> | null = null,
  pepseaMethod: typeof pepseaMethods[number] = pepseaMethods[0]
): Promise<DG.Column> {
  return new Promise(async (resolve, reject) => {
    const table = col?.dataFrame ?? grok.shell.t;
    const seqCol = col ?? table.columns.bySemType(DG.SEMTYPE.MACROMOLECULE);
    if (seqCol == null) {
      const errMsg = `MSAError: dataset doesn't conain any Macromolecule column`;
      grok.shell.warning(errMsg);
      reject(new MsaWarning(errMsg));
    }

    // UI
    const methodInput = ui.choiceInput('Method', pepseaMethod, pepseaMethods);
    methodInput.setTooltip('Alignment method');
    const gapOpenInput = ui.floatInput('Gap open', 1.53);
    gapOpenInput.setTooltip('Gap opening penalty at group-to-group alignment');
    const gapExtendInput = ui.floatInput('Gap extend', 0);
    gapExtendInput.setTooltip('Gap extension penalty to skip the alignment');
    const inputRootStyles = [methodInput.root.style, gapOpenInput.root.style, gapExtendInput.root.style];
    let performAlignment: (() => Promise<DG.Column<string>>) | undefined;

    // TODO: allow only macromolecule colums to be chosen
    const colInput = ui.columnInput('Sequence', table, seqCol, () => {
      performAlignment = onColInputChange(
        colInput.value,
        table,
        inputRootStyles,
        methodInput,
        clustersColInput,
        gapOpenInput,
        gapExtendInput
      );
    }
    ) as DG.InputBase<DG.Column<string>>;
    colInput.setTooltip('Sequences column to use for alignment');
    const clustersColInput = ui.columnInput('Clusters', table, null);
    clustersColInput.nullable = true;
    colInput.fireChanged();
    //if column is specified (from tests), run alignment and resolve with the result
    if (col) {
      performAlignment = onColInputChange(
        col,
        table,
        inputRootStyles,
        methodInput,
        clustersColInput,
        gapOpenInput,
        gapExtendInput
      );

      await onDialogOk(colInput, table, performAlignment, resolve, reject);
      return;
    }
    const dlg = ui.dialog('MSA')
      .add(colInput)
      .add(clustersColInput)
      .add(methodInput)
      .add(gapOpenInput)
      .add(gapExtendInput)
      .onOK(async () => {
        await onDialogOk(colInput, table, performAlignment, resolve, reject);
      })
      .show();
  });
}

async function onDialogOk(
  colInput: DG.InputBase< DG.Column<any>>,
  table: DG.DataFrame,
  performAlignment: (() => Promise<DG.Column<string>>) | undefined,
  resolve: (value: DG.Column<any>) => void,
  reject: (reason: any) => void
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
      return grok.shell.warning('Wrong column format');

    table.columns.add(msaCol);
    await grok.data.detectSemanticTypes(table);

    resolve(msaCol);
  } catch (err: any) {
    const errMsg: string = err instanceof Error ? err.message : err.toString();
    grok.shell.error(errMsg);
    reject(err);
  } finally {
    pi.close();
  }
}


function onColInputChange(
  col: DG.Column<string>,
  table: DG.DataFrame,
  inputRootStyles: CSSStyleDeclaration[],
  methodInput: DG.InputBase<string | null>,
  clustersColInput: DG.InputBase<DG.Column<any> | null>,
  gapOpenInput: DG.InputBase<number | null>,
  gapExtendInput: DG.InputBase<number | null>
): (() => Promise<DG.Column<string>>) | undefined {
  try {
    if (col.semType !== DG.SEMTYPE.MACROMOLECULE)
      return;
    const unusedName = table.columns.getUnusedName(`msa(${col.name})`);

    if (checkInputColumnUI(col, col.name,
      [NOTATION.FASTA, NOTATION.SEPARATOR], [ALPHABET.DNA, ALPHABET.RNA, ALPHABET.PT], false)
    ) { // Kalign - natural alphabets. if the notation is separator, convert to fasta and then run kalign
      for (const inputRootStyle of inputRootStyles)
        inputRootStyle.display = 'none';
      const potentialColNC = new NotationConverter(col);
      const performCol: DG.Column<string> = potentialColNC.isFasta() ? col :
        potentialColNC.convert(NOTATION.FASTA);
      return async () => await runKalign(performCol, false, unusedName, clustersColInput.value);
    } else if (checkInputColumnUI(col, col.name,
      [NOTATION.HELM], [], false)
    ) { // PepSeA branch - Helm notation or separator notation with unknown alphabets
      for (const inputRootStyle of inputRootStyles)
        inputRootStyle.removeProperty('display');

      return async () => await runPepsea(col, unusedName, methodInput.value!,
          gapOpenInput.value!, gapExtendInput.value!, clustersColInput.value);
    } else {
      for (const inputRootStyle of inputRootStyles)
        inputRootStyle.display = 'none';

      return;
    }
  } catch (err: any) {
    const errMsg: string = err instanceof Error ? err.message : err.toString();
    grok.shell.error(errMsg);
    _package.logger.error(errMsg);
  }
}
