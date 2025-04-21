import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {startAnalysis} from '../widgets/peptides';
import * as C from '../utils/constants';
import {scaleActivity} from '../utils/misc';
import {ALIGNMENT, ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {PeptidesModel} from '../model';
import {awaitCheck, delay} from '@datagrok-libraries/utils/src/test';
import {MCLSettings} from '../utils/types';

export async function macromoleculeSarFastaDemoUI(): Promise<void> {
  grok.shell.windows.showContextPanel = true;
  const alignedSequenceCol = 'AlignedSequence';
  const simpleActivityColName = 'IC50';
  const simpleTable = DG.DataFrame.fromCsv(await _package.files.readAsText('aligned.csv'));
  simpleTable.name = 'Simple peptides';
  if (!simpleTable.id)
    simpleTable.id = `Simple-peptides-analysis-table-id-${Math.random()}-${Date.now()}`;
  const view = grok.shell.addTableView(simpleTable);
  await delay(50);
  const grid = view.grid;
  await new Promise((res) => {
    let timer: any = null;
    const sub = grid?.onAfterDrawContent?.subscribe(() => {
      if (timer)
        clearTimeout(timer);
      sub.unsubscribe();
      res(undefined);
    });
    timer = setTimeout(() => {
      sub.unsubscribe();
      res(undefined);
    }, 3000);
  });
  const simpleActivityCol = simpleTable.getCol(simpleActivityColName);
  const simpleAlignedSeqCol = simpleTable.getCol(alignedSequenceCol);
  simpleAlignedSeqCol.semType = DG.SEMTYPE.MACROMOLECULE;
  simpleAlignedSeqCol.setTag(C.TAGS.ALPHABET, ALPHABET.PT);
  simpleAlignedSeqCol.meta.units = NOTATION.FASTA;
  simpleAlignedSeqCol.setTag(bioTAGS.aligned, ALIGNMENT.SEQ_MSA);
  const simpleScaledCol = scaleActivity(simpleActivityCol, C.SCALING_METHODS.MINUS_LG);
  const mclSettings = new MCLSettings();
  mclSettings.threshold = 94;
  const model: PeptidesModel = await startAnalysis(simpleActivityCol, simpleAlignedSeqCol,
    null, simpleTable, simpleScaledCol,
    C.SCALING_METHODS.MINUS_LG, {addMCL: true, useEmbeddingsClusters: true, mclSettings: mclSettings}) as PeptidesModel;
  // model.modifyWebLogoSelection({monomerOrCluster: 'D', positionOrClusterType: '13'});
}
