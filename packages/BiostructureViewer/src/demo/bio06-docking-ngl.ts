import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';
import {handleError} from './utils';
import {INglViewer} from '@datagrok-libraries/bio/src/viewers/ngl-gl-viewer';

import {awaitGrid} from '../tests/utils';
import {NglViewer} from '../viewers/ngl-viewer';
import {delay} from '@datagrok-libraries/test/src/test';

import {_package} from '../package';

const ligandsDataFn: string = 'CHEMBL2366517/ic50.csv';
const structureDataFn: string = 'samples/1bdq-wo-ligands.pdb';

const colNameDesc: { [colName: string]: { name: string, description: string, semType?: string } } = {
  'mol': {name: 'Mol', description: 'Pose, docking result', semType: DG.SEMTYPE.MOLECULE},
  'id': {name: 'id', description: 'Molecule ChEMBL ID'},
  'MW': {name: 'MW', description: 'Molecular Weight'},
  'affinity': {name: 'Affinity', description: 'Estimated Free Energy of Binding'},
  'intermolecular': {name: 'Intermolecular', description: 'Final Intermolecular Energy'},
  'electrostatic': {name: 'Electrostatic', description: 'Electrostatic Energy'},
  'ligand-fixed': {name: 'Ligand-Fixed', description: 'Moving Ligand-Fixed Receptor'},
  'ligand-moving': {name: 'Ligand-Moving', description: 'Moving Ligand-Moving Receptor'},
  'total internal': {name: 'Total Internal', description: 'Final Total Internal Energy'},
  'torsional free': {name: 'Torsional Free', description: 'Torsional Free Energy'},
  'unbound systems': {name: 'Unbound System\'s', description: 'Unbound System\'s Energy'},
};

export async function demoBio06NoScript(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('Demo Docking Conformations ...');
  try {
    const [dfCsv, layoutStr] = await Promise.all([
      // Read through function to cache
      grok.functions.call(`${_package.name}:readAsTextDapi`,
        {file: 'System:AppData/BiostructureViewer/CHEMBL2366517/ic50.pose.flt2.mol.log-log.2.csv'}),
      // the structure for the viewer is stored within the layout
      _package.files.readAsText('CHEMBL2366517/demo-docking-conformations.layout'),
    ]);
    if (!dfCsv)
      throw new Error('Empty data');

    const df = DG.DataFrame.fromCsv(dfCsv);
    for (const [colName, colConfig] of Object.entries(colNameDesc)) {
      const col = df.col(colName);
      if (col) {
        col.setTag('friendlyName', colConfig.name);
        col.setTag('description', colConfig.description);
        if (colConfig.semType)
          col.semType = colConfig.semType;
      }
    }

    // df.getCol('mol').temp['regenerate-coords'] = 'true'; // TODO: flat and beautiful molecules
    const view = grok.shell.addTableView(df);

    const layout = DG.ViewLayout.fromJson(layoutStr);
    view.loadLayout(layout);
    // _package.files.readAsText('CHEMBL2366517/demo-docking-conformations.layout').then((layoutStr: string) => {
    //   const layout = DG.ViewLayout.fromJson(layoutStr);
    //   view.loadLayout(layout);
    // });
    await awaitGrid(view.grid, 5000);
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
