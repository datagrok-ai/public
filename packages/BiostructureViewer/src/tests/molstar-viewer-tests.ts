import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Base64} from 'js-base64';

import {category, expect/*, expect*/, test} from '@datagrok-libraries/utils/src/test';

import {MolstarViewerApp} from '../apps/molstar-viewer-app';
import {awaitGrid} from './utils';
import {buildDataJson, MolstarViewer, PROPS as msvPROPS} from '../viewers/molstar-viewer/molstar-viewer';

import {_package} from '../package-test';


category('MolstarViewer', () => {
  test('ligands-Molecule', async () => {
    const [sdfCnt, pdbData] = await Promise.all([
      _package.files.readAsBytes('samples/1bdq-obs-pred.sdf'),
      _package.files.readAsBytes('samples/1bdq-wo-ligands.pdb'),
    ]);
    const df: DG.DataFrame = (await grok.functions.call('Chem:importSdf', {bytes: sdfCnt}))[0];
    const ligandCol = df.getCol('molecule');

    const view = grok.shell.addTableView(df);
    await awaitGrid(view.grid);

    const viewer = (await df.plot.fromType('Biostructure', {
      [msvPROPS.dataJson]: buildDataJson(pdbData, 'pdb'),
      [msvPROPS.ligandColumnName]: ligandCol.name,
      [msvPROPS.showCurrentRowLigand]: true,
      [msvPROPS.showSelectedRowsLigands]: true,
    })) as unknown as MolstarViewer;
    view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Biostructure', 0.4);

    // The currentRowIdx is set to 0 by default at start
    // df.currentRowIdx = -1;
    // await awaitGrid(view.grid);
    // await Promise.all([viewer.viewPromise, viewer.ligandsPromise]);
    // expect(viewer.ligands.current === null, true, 'No current ligand expected');
    // expect(viewer.ligands.selected.length, 0);

    df.currentRowIdx = 0;
    await Promise.all([awaitGrid(view.grid), viewer.viewPromise]);
    expect(viewer.ligands.current !== null, true, 'The current ligand expected');
    expect(viewer.ligands.selected.length, 0);

    df.selection.init((rowI) => rowI === 1 || rowI === 2);
    await Promise.all([awaitGrid(view.grid), viewer.viewPromise]);
    expect(viewer.ligands.current != null, true, 'The current ligand expected');
    expect(viewer.ligands.selected.length, 2);
  });

  test('ligands-Molecule3D', async () => {
    const [pdbqtCnt, targetData] = await Promise.all([
      _package.files.readAsText('docking/ligand_out.pdbqt'),
      _package.files.readAsBytes('docking/3SWZ.pdbqt'),
    ]);
    const df: DG.DataFrame = (await grok.functions.call('BiostructureViewer:importPdbqt',
      {fileContent: pdbqtCnt, test: true}))[0];
    const ligandCol = df.getCol('molecule');

    const view = grok.shell.addTableView(df);
    await awaitGrid(view.grid);

    const viewer = (await df.plot.fromType('Biostructure', {
      [msvPROPS.dataJson]: buildDataJson(targetData, 'pdbqt'),
      [msvPROPS.ligandColumnName]: ligandCol.name,
      [msvPROPS.showCurrentRowLigand]: true,
      [msvPROPS.showSelectedRowsLigands]: true,
    })) as unknown as MolstarViewer;
    view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Biostructure', 0.4);

    df.currentRowIdx = 0;
    await Promise.all([awaitGrid(view.grid), viewer.viewPromise]);
    expect(viewer.ligands.current !== null, true, 'The current ligand expected');
    expect(viewer.ligands.selected.length, 0);

    df.selection.init((rowI) => rowI === 1 || rowI === 2);
    await Promise.all([awaitGrid(view.grid), viewer.viewPromise]);
    expect(viewer.ligands.current != null, true, 'The current ligand expected');
    expect(viewer.ligands.selected.length, 2);
  });
});
