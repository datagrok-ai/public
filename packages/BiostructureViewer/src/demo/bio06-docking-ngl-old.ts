import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';
import {handleError} from './utils';
import {INglViewer} from '@datagrok-libraries/bio/src/viewers/ngl-gl-viewer';

import {_package} from '../package';
import {awaitGrid} from '../tests/utils';
import {NglViewer} from '../viewers/ngl-viewer';
import {delay} from '@datagrok-libraries/test/src/test';

const ligandsDataFn: string = 'samples/1bdq-obs-pred.sdf';
const structureDataFn: string = 'samples/1bdq-wo-ligands.pdb';

export async function demoBio06NoScriptOld(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('Demo Docking Conformations ...');
  try {
    const sdfBytes: Uint8Array = await _package.files.readAsBytes(ligandsDataFn);
    const df = (await grok.functions.call(
      'Chem:importSdf', {bytes: sdfBytes}))[0];
    const view = grok.shell.addTableView(df);

    const pdbStr: string = await _package.files.readAsText(structureDataFn);

    const viewer = await df.plot.fromType('NGL', {
      pdb: pdbStr,
      ligandColumnName: 'molecule',
      showSelectedRowsLigands: true,
      showCurrentRowLigand: true,
      showMouseOverRowLigand: true,
    }) as unknown as NglViewer;
    view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Biostructure', 0.62);

    grok.shell.windows.showHelp = true;
    // TODO: Dependency on datagrok-api ^1.15.0
    // @ts-ignore
    if (grok.shell.windows.help) {
      // @ts-ignore
      grok.shell.windows.help.showHelp(viewer.helpUrl);
    }

    await Promise.all([awaitGrid(view.grid), viewer.awaitRendered()]);
  } finally {
    pi.close();
  }
}

export async function demoBio06UI(): Promise<void> {
  let view: DG.TableView;
  let df: DG.DataFrame;
  let viewer: DG.Viewer & INglViewer;

  try {
    await new DemoScript(
      'Docking NGL',
      'Docking ligands along the structure'
    )
      .step('Loading ligands', async () => {
        grok.shell.windows.showContextPanel = false;
        grok.shell.windows.showProperties = false;

        const sdfBytes: Uint8Array = await _package.files.readAsBytes(ligandsDataFn);
        df = (await grok.functions.call(
          'Chem:importSdf', {bytes: sdfBytes}))[0];

        view = grok.shell.addTableView(df);
      }, {
        description: `Load dataset wih ligand structures (sdf).`,
        delay: 2000,
      })
      .step('NGL viewer with structure', async () => {
        const pdbStr: string = await _package.files.readAsText(structureDataFn);

        viewer = await df.plot.fromType('NGL', {
          pdb: pdbStr,
          ligandColumnName: 'molecule',
          showSelectedRowsLigands: true,
          showCurrentRowLigand: true,
          showMouseOverRowLigand: true,
        });
        view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Biostructure', 0.62);
      }, {
        description: `Add NGL viewer to display PDB. ` +
          `Try to change the current, selected ligands to see them along the structure.`,
        delay: 2000,
      })
      /*
      .step('Tracking selected, current', async () => {
        df.selection.init((rowI: number) => [3, 5, 6].includes(rowI));
        df.currentRowIdx = 1;
      }, {
        description: `Display selected and current ligands along the structure with Datagrok colors.`,
        delay: 2000,
      })
      .step('Tracking mouse over', async () => {
        df.selection.init((rowI) => false);
        df.currentRowIdx = -1;

        for (let rowI: number = 0; rowI < 5; rowI++) {
          df.mouseOverRowIdx = rowI;
          await delay(300);
        }
      }, {
        description: `Display mouse overed ligands along the structure.`,
        delay: 1000,
      })
      /**/
      .start();
  } catch (err: any) {
    handleError(err);
  }
}
