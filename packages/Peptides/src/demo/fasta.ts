import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {startAnalysis} from '../widgets/peptides';
import * as C from '../utils/constants';
import {scaleActivity} from '../utils/misc';
import {ALIGNMENT, ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';

export async function macromoleculeSarFastaDemoUI(): Promise<void> {
  const alignedSequenceCol = 'AlignedSequence';
  const simpleActivityColName = 'IC50';
  const simpleTable = DG.DataFrame.fromCsv(await _package.files.readAsText('aligned.csv'));
  const simpleActivityCol = simpleTable.getCol(simpleActivityColName);
  const simpleAlignedSeqCol = simpleTable.getCol(alignedSequenceCol);
  simpleAlignedSeqCol.semType = DG.SEMTYPE.MACROMOLECULE;
  simpleAlignedSeqCol.setTag(C.TAGS.ALPHABET, ALPHABET.PT);
  simpleAlignedSeqCol.setTag(DG.TAGS.UNITS, NOTATION.FASTA);
  simpleAlignedSeqCol.setTag(bioTAGS.aligned, ALIGNMENT.SEQ_MSA);
  const simpleScaledCol = scaleActivity(simpleActivityCol, '-lg');

  const model = await startAnalysis(simpleActivityCol, simpleAlignedSeqCol, null, simpleTable, simpleScaledCol, '-lg');
  model?.analysisView.addViewer('WebLogo');
}
