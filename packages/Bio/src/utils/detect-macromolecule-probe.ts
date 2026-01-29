import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


import {_package} from '../package';
import {delay} from '@datagrok-libraries/test/src/test';

type IDetectorReport = { categoriesSample: any[], rejectReason: string };

type IDetectorDebugStore = { last: IDetectorReport };

export async function detectMacromoleculeProbeDo(
  csv: string, colName: string | undefined, probeCount: number
): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create(`detectMacromolecule probe ...`);
  try {
    let progressLast = 0;
    const store: IDetectorDebugStore = await grok.functions.call('Bio:detectMacromoleculeEnableStore');
    let failCount: number = 0;
    for (let i = 0; i < probeCount; ++i) {
      const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
      const seqCol = colName ? df.getCol(colName) : df.columns.byIndex(0);

      const detectRes: string = await grok.functions.call('Bio:detectMacromolecule', {col: seqCol});
      if (detectRes !== DG.SEMTYPE.MACROMOLECULE) {
        ++failCount;
        console.warn(`Reject reason: ${store.last.rejectReason}`);
      }
      const progress = i / probeCount;
      if ((progress - progressLast) >= 0.1) {
        progressLast = progress;
        pi.update(100 * progress, `detectMacromolecule probe ${failCount}/${i}/${probeCount} ...`);
        await delay(0);
      }
    }
    if (failCount > 0)
      grok.shell.warning(`detectMacromolecule failed ${failCount} / ${probeCount}`);
    else
      grok.shell.info(`detectMacromolecule success ${probeCount}`);
  } finally {
    pi.close();
  }
}
