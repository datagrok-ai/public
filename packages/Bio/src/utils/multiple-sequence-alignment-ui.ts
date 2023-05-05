import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {ALPHABET, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {runKalign} from './multiple-sequence-alignment';
import {pepseaMethods, runPepsea} from './pepsea';
import {checkInputColumnUI} from './check-input-column';

export function multipleSequenceAlignmentUI(col: DG.Column<string> | null = null) {
  const table = col?.dataFrame ?? grok.shell.t;
  const seqCol = col ?? table.columns.bySemType(DG.SEMTYPE.MACROMOLECULE);
  if (seqCol == null) {
    grok.shell.warning(`MSAError: dataset doesn't conain any Macromolecule column`);
    return;
  }

  let performAlignment: () => Promise<DG.Column<string> | null> = async () => null;
  const methodInput = ui.choiceInput('Method', pepseaMethods[0], pepseaMethods);
  methodInput.setTooltip('Alignment method');
  const gapOpenInput = ui.floatInput('Gap open', 1.53);
  gapOpenInput.setTooltip('Gap opening penalty at group-to-group alignment');
  const gapExtendInput = ui.floatInput('Gap extend', 0);
  gapExtendInput.setTooltip('Gap extension penalty to skip the alignment');
  const inputRootStyles = [methodInput.root.style, gapOpenInput.root.style, gapExtendInput.root.style];

  const colInput = ui.columnInput('Sequence', table, seqCol, () => {
    const potentialCol = colInput.value;
    const unusedName = table.columns.getUnusedName(`msa(${potentialCol.name})`);

    if (checkInputColumnUI(
      potentialCol, potentialCol.name, [NOTATION.FASTA], [ALPHABET.DNA, ALPHABET.RNA, ALPHABET.PT], false)) {
      for (const inputRootStyle of inputRootStyles)
        inputRootStyle.display = 'none';

      performAlignment = () => runKalign(potentialCol, false, unusedName, clustersColInput.value);
    } else if (checkInputColumnUI(potentialCol, potentialCol.name, [NOTATION.HELM], [], false)) {
      for (const inputRootStyle of inputRootStyles)
        inputRootStyle.display = 'initial';

      performAlignment = () => runPepsea(potentialCol, unusedName, methodInput.value!, gapOpenInput.value!,
        gapExtendInput.value!, clustersColInput.value);
    } else {
      for (const inputRootStyle of inputRootStyles)
        inputRootStyle.display = 'none';

      performAlignment = async () => null;
    }
  }) as DG.InputBase<DG.Column<string>>;
  colInput.setTooltip('Sequences column to use for alignment');
  colInput.fireChanged();

  const clustersColInput = ui.columnInput('Clusters', table, null);
  clustersColInput.nullable = true;

  let msaCol: DG.Column<string> | null = null;
  ui.dialog('MSA')
    .add(colInput)
    .add(clustersColInput)
    .add(methodInput)
    .add(gapOpenInput)
    .add(gapExtendInput)
    .onOK(async () => {
      const pi = DG.TaskBarProgressIndicator.create('Analyze for MSA ...');
      try {
        colInput.fireChanged();
        msaCol = await performAlignment(); // progress
        if (msaCol == null)
          return grok.shell.warning('Wrong column format');

        table.columns.add(msaCol);
        await grok.data.detectSemanticTypes(table);
      } catch (err: any) {
        const errMsg: string = err instanceof Error ? err.message : err.toString();
        grok.shell.error(errMsg);
      } finally {
        pi.close();
      }
    })
    .show();
}
