import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package, sequenceSpaceTopMenu} from '../package';
import {handleError} from './utils';

import {IWebLogoViewer} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {pepseaMethods, runPepsea} from '../utils/pepsea';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/reduce-dimensionality';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';

const helmFn: string = 'samples/HELM.csv';

export async function demoBio05UI(): Promise<void> {
  let view: DG.TableView;
  let df: DG.DataFrame;
  let helmCol: DG.Column<string>;
  let msaHelmCol: DG.Column<string>;
  let wlViewer: DG.Viewer & IWebLogoViewer;
  let ssViewer: DG.ScatterPlotViewer;

  const helmColName: string = 'HELM';
  const msaHelmColName: string = 'msa(HELM)';
  const dimRedMethod: DimReductionMethods = DimReductionMethods.UMAP;

  try {
    const demoScript = new DemoScript(
      'Helm, MSA, Sequence Space',
      'MSA and composition analysis on Helm data');
    await demoScript
      .step(`Load peptides with non-natural aminoacids in 'HELM' notation`, async () => {
        view = grok.shell.addTableView(df = await _package.files.readCsv(helmFn));

        grok.shell.windows.showContextPanel = false;
        grok.shell.windows.showProperties = false;
      }, {
        description: 'Load dataset with macromolecules of \'Helm\' notation.',
        delay: 2000,
      })
      .step('Align peptides with non-natural aminoacids with PepSeA', async () => {
        helmCol = df.getCol(helmColName);
        const method: string = pepseaMethods[0];
        const gapOpen: number = 1.53;
        const gapExtend: number = 0;
        msaHelmCol = (await runPepsea(helmCol, msaHelmColName, method, gapOpen, gapExtend, undefined))!;
        df.columns.add(msaHelmCol);
        await grok.data.detectSemanticTypes(df);
      }, {
        // eslint-disable-next-line max-len
        description: 'Multiple sequence alignment (MSA) performed with PepSeA tool operating on non-natural aminoacids as well.',
        delay: 2000,
      })
      .step('Build sequence space', async () => {
        ssViewer = (await sequenceSpaceTopMenu(df, msaHelmCol,
          dimRedMethod, MmDistanceFunctionsNames.LEVENSHTEIN, true)) as DG.ScatterPlotViewer;
        view.dockManager.dock(ssViewer, DG.DOCK_TYPE.RIGHT, null, 'Sequence Space', 0.35);
      }, {
        description: 'Reduce sequence space dimensionality to display on 2D representation.',
        delay: 2000,
      })
      .step('Analyse sequence composition', async () => {
        wlViewer = await df.plot.fromType('WebLogo', {
          sequenceColumnName: msaHelmColName,
          positionWidth: 40,
          maxHeight: 50,
        }) as DG.Viewer & IWebLogoViewer;
        view.dockManager.dock(wlViewer, DG.DOCK_TYPE.DOWN, null, 'Composition analysis', 0.2);
      }, {
        // eslint-disable-next-line max-len
        description: 'Composition analysis allows to reveal functional features of sequences like motifs, or variable loops.',
        delay: 2000,
      })
      .start();
  } catch (err: any) {
    handleError(err);
  }
}
