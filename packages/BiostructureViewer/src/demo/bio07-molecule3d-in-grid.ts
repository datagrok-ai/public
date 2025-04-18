import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IBiostructureViewer} from '@datagrok-libraries/bio/src/viewers/molstar-viewer';
import {handleError} from './utils';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';

import {_package} from '../package';
import {MolstarViewer} from '../viewers/molstar-viewer';
import {awaitGrid} from '../tests/utils';

const pdbCsvFn: string = 'pdb_data.csv';
const pdbColName: string = 'pdb';

export async function demoBio07NoScript(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('Demo Proteins ...');
  try {
    grok.shell.windows.showContextPanel = false;
    grok.shell.windows.showProperties = false;

    let csv: string;
    let df: DG.DataFrame;
    try {
      csv = await _package.files.readAsText(pdbCsvFn);
      if (!csv)
        throw new Error('Data file is empty.');
      df = DG.DataFrame.fromCsv(csv);
    } catch (err: any) {
      grok.shell.warning(`Error reading file '${pdbCsvFn}': ${err.toString()}`);
      // fallback on empty data
      const idCol = DG.Column.fromStrings('id', []);
      const pdbIdCol = DG.Column.fromStrings('pdb_id', []);
      const pdbCol = DG.Column.fromStrings('pdb', []);
      pdbCol.semType = DG.SEMTYPE.MOLECULE3D;
      df = DG.DataFrame.fromColumns([idCol, pdbIdCol, pdbCol]);
    }

    const view = grok.shell.addTableView(df);
    view.grid.columns.byName('id')!.width = 0;

    df.currentCell = df.cell(0, pdbColName);
    const pdbStr: string = df.currentCell.value;
    const viewer = (await df.plot.fromType('Biostructure', {
      pdb: pdbStr,
      biostructureIdColumnName: pdbColName,
    })) as unknown as MolstarViewer;
    view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Biostructure', 0.5);

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

export async function demoBio07UI(): Promise<void> {
  let view: DG.TableView;
  let df: DG.DataFrame;
  let viewer: DG.Viewer & IBiostructureViewer;

  try {
    await new DemoScript(
      'Molecule3D in Grid',
      'View structures PDB in grid',
    )
      .step('Loading structures', async () => {
        grok.shell.windows.showContextPanel = false;
        grok.shell.windows.showProperties = false;

        df = await _package.files.readCsv(pdbCsvFn);
        view = grok.shell.addTableView(df);
        view.grid.columns.byName('id')!.width = 0;
      }, {
        description: 'Load dataset with structures (PDB).',
        delay: 2000,
      })
      .step('Biostructure viewer', async () => {
        df.currentCell = df.cell(0, pdbColName);
        const pdbStr: string = df.currentCell.value;
        viewer = (await df.plot.fromType('Biostructure', {
          pdb: pdbStr,
        }));
        view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Biostructure', 0.5);
      }, {
        description: `Add Biostructure viewer`,
        delay: 2000,
      })
      .step('Tracking PDB cell', async () => {
        df.currentCell = df.cell(2, pdbColName);
        const pdbStr: string = df.currentCell.value;
        viewer.setOptions({pdb: pdbStr});
      }, {
        description: `'Molecule3D' cell renderer handle mouse click displaying data with the BiostructureViewer.`,
        delay: 2000,
      })
      .step('Tracking PDB cell', async () => {
        df.currentCell = df.cell(1, pdbColName);
        const pdbStr: string = df.currentCell.value;
        viewer.setOptions({pdb: pdbStr});
      }, {
        description: `'Molecule3D' cell renderer handle mouse click displaying data with the BiostructureViewer.`,
        delay: 2000,
      })
      .start();
  } catch (err: any) {
    handleError(err);
  }
}
