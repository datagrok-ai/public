import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {startAnalysis} from '../widgets/peptides';
import * as C from '../utils/constants';
import {scaleActivity} from '../utils/misc';
import {ALIGNMENT, ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {PeptidesModel} from '../model';
import {delay} from '@datagrok-libraries/utils/src/test';

export async function macromoleculeSarFastaDemoUI(): Promise<void> {
  grok.shell.windows.showContextPanel = true;
  const alignedSequenceCol = 'AlignedSequence';
  const simpleActivityColName = 'IC50';
  const simpleTable = DG.DataFrame.fromCsv(await _package.files.readAsText('aligned.csv'));
  grok.shell.addTableView(simpleTable);
  await delay(0);
  const simpleActivityCol = simpleTable.getCol(simpleActivityColName);
  const simpleAlignedSeqCol = simpleTable.getCol(alignedSequenceCol);
  simpleAlignedSeqCol.semType = DG.SEMTYPE.MACROMOLECULE;
  simpleAlignedSeqCol.setTag(C.TAGS.ALPHABET, ALPHABET.PT);
  simpleAlignedSeqCol.meta.units = NOTATION.FASTA;
  simpleAlignedSeqCol.setTag(bioTAGS.aligned, ALIGNMENT.SEQ_MSA);
  const simpleScaledCol = scaleActivity(simpleActivityCol, C.SCALING_METHODS.MINUS_LG);
  const clustersCol = DG.Column.string('Cluster', simpleTable.rowCount).init('0');
  const alignedCol: DG.Column<string> = await grok.functions.call('Bio:alignSequences', {
    sequenceCol: simpleAlignedSeqCol,
    clustersCol,
  });
  const model: PeptidesModel = await startAnalysis(simpleActivityCol, alignedCol, null, simpleTable, simpleScaledCol,
    C.SCALING_METHODS.MINUS_LG) as PeptidesModel;
  model.modifyWebLogoSelection({monomerOrCluster: 'D', positionOrClusterType: '13'});
}
