import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IWebLogoViewer} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {awaitStatus, DockerContainerStatus} from '@datagrok-libraries/bio/src/utils/docker';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';

import {Pepsea, pepseaMethods, runPepsea} from '../utils/pepsea';
import {sequenceSpaceTopMenu} from '../package';
import {handleError} from './utils';

import {_package} from '../package';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';

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

  const pepseaDcId = (await Pepsea.getDockerContainer()).id;
  // // region For test: Stop container to test auto-start
  // await grok.dapi.docker.dockerContainers.stop(pepseaDcId);
  // await Pepsea.awaitStatus(pepseaDcId, 'stopped', 15000);
  // // endregion
  const pepseaDcPromise: Promise<DG.DockerContainer> = Pepsea.getDockerContainer();
  let pepseaDcStatus: DockerContainerStatus;
  let pepseaDcStartPromise: Promise<void>;

  try {
    const demoScript = new DemoScript('Helm, MSA, Sequence Space', 'MSA and composition analysis on Helm data',
      false, {autoStartFirstStep: true});
    await demoScript
      .step(`Load peptides with non-natural aminoacids in 'HELM' notation`, async () => {
        [pepseaDcStatus, df] = await Promise.all([
          (async () => { return (await pepseaDcPromise).status; })(),
          _package.files.readCsv(helmFn)
        ]);
        view = grok.shell.addTableView(df);

        grok.shell.windows.showContextPanel = false;
        grok.shell.windows.showProperties = false;

        if (pepseaDcStatus === 'started' || pepseaDcStatus === 'checking') {
          _package.logger.debug(
            `demoBio05UI(), PepSeA ('${Pepsea.dcName}') docker container status = '${pepseaDcStatus}'.`);
          pepseaDcStartPromise = Promise.resolve();
        } else {
          _package.logger.warning(
            `demoBio05UI(), PepSeA ('${Pepsea.dcName}') docker container is trying to start...`);

          await grok.dapi.docker.dockerContainers.run(pepseaDcId);
          pepseaDcStartPromise = awaitStatus(pepseaDcId, 'started', 30000, _package.logger);
        }
      }, {
        description: 'Load dataset with macromolecules of \'Helm\' notation.',
        delay: 2000,
      })
      .step('Align peptides with non-natural aminoacids with PepSeA', async () => {
        const pi = DG.TaskBarProgressIndicator.create('MSA by PepSeA ...');
        try {
          // TODO: Show splash if pepseaDcStartPromise is not resolved still
          await pepseaDcStartPromise; // throws timeout
          // Hide splash

          helmCol = df.getCol(helmColName);
          const method: string = pepseaMethods[0];
          const gapOpen: number = 1.53;
          const gapExtend: number = 0;
          msaHelmCol = (await runPepsea(helmCol, msaHelmColName, method, gapOpen, gapExtend, undefined))!;
          if (!msaHelmCol)
            throw new Error(`Empty MSA result.`);
          df.columns.add(msaHelmCol);
          await grok.data.detectSemanticTypes(df);
        } finally {
          pi.close();
        }
      }, {
        // eslint-disable-next-line max-len
        description: 'Multiple sequence alignment (MSA) performed with PepSeA tool operating on non-natural aminoacids as well.',
        delay: 2000,
      })
      .step('Build sequence space', async () => {
        const preprocessingFunc = DG.Func.find({package: 'Bio', name: 'macromoleculePreprocessingFunction'})[0];
        ssViewer = (await sequenceSpaceTopMenu(df, msaHelmCol,
          dimRedMethod, MmDistanceFunctionsNames.LEVENSHTEIN, true, preprocessingFunc, undefined, undefined, true)) as DG.ScatterPlotViewer;
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
