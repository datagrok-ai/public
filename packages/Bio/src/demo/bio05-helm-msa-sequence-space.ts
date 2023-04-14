import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package, sequenceSpaceTopMenu} from '../package';
import {step} from './utils';

import {IWebLogoViewer} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {pepseaMethods, runPepsea} from '../utils/pepsea';
import {StringMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';

const helmFn: string = 'samples/sample_HELM.csv';

export async function demoBio05UI(funcPath: string): Promise<void> {
  let view: DG.TableView;
  let df: DG.DataFrame;
  let helmCol: DG.Column<string>;
  let msaHelmCol: DG.Column<string>;
  let wlViewer: DG.Viewer & IWebLogoViewer;
  let ssViewer: DG.ScatterPlotViewer;

  const helmColName: string = 'HELM';
  const msaHelmColName: string = 'msa(HELM)';

  try {
    await step(`Loading peptides notation 'HELM'.`, async () => {
      view = grok.shell.addTableView(df = await _package.files.readCsv(helmFn));
      view.path = view.basePath = funcPath;
    })();

    await step('MSA on non-natural aminoacids with PepSeA.', async () => {
      helmCol = df.getCol(helmColName);
      const method: string = pepseaMethods[0];
      const gapOpen: number = 1.53;
      const gapExtend: number = 0;
      msaHelmCol = await runPepsea(helmCol, msaHelmColName, method, gapOpen, gapExtend, undefined);
      df.columns.add(msaHelmCol);
    })();

    await step('Composition analysis on MSA results', async () => {
      wlViewer = await df.plot.fromType('WebLogo', {
        sequenceColumnName: msaHelmColName
      }) as DG.Viewer & IWebLogoViewer;
      view.dockManager.dock(wlViewer, DG.DOCK_TYPE.DOWN, null, 'Composition analysis', 0.2);
    })();

    await step('Building sequence space.', async () => {
      const method: string = 'UMAP';
      ssViewer = (await sequenceSpaceTopMenu(df, msaHelmCol,
        'UMAP', StringMetricsNames.Levenshtein, true)) as DG.ScatterPlotViewer;
      view.dockManager.dock(ssViewer, DG.DOCK_TYPE.RIGHT, null, 'Sequence Space', 0.35);
    })();
  } catch (err: any) {
    if (err instanceof Error)
      _package.logger.error(err.message, undefined, err.stack);
    else
      _package.logger.error(err.toString());
  }
}
