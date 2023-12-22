import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, delay, expect/*, expect*/, test} from '@datagrok-libraries/utils/src/test';
import {BiostructureData, BiostructureDataJson} from '@datagrok-libraries/bio/src/pdb/types';

import {awaitGrid} from './utils';
import {DebounceIntervals, MolstarViewer, PROPS as msvPROPS} from '../viewers/molstar-viewer/molstar-viewer';

import {_package} from '../package-test';


category('MolstarViewer', () => {
  test('ligands-Molecule', async () => {
    const [sdfCnt, pdbData]: [Uint8Array, BiostructureData] = await Promise.all([
      _package.files.readAsBytes('samples/1bdq-obs-pred.sdf'),
      (async () => {
        return {binary: true, ext: 'pdb', data: await _package.files.readAsBytes('samples/1bdq-wo-ligands.pdb')};
      })(),
    ]);
    const df: DG.DataFrame = (await grok.functions.call('Chem:importSdf', {bytes: sdfCnt}))[0];
    const ligandCol = df.getCol('molecule');

    const view = grok.shell.addTableView(df);
    await awaitGrid(view.grid);

    const viewer = (await df.plot.fromType('Biostructure', {
      [msvPROPS.dataJson]: BiostructureDataJson.fromData(pdbData),
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
    await Promise.all([awaitGrid(view.grid), viewer.awaitRendered()]);
    await delay(DebounceIntervals.ligands * 2.5); // await for debounce onRebuildViewLigands
    expect(viewer.ligands.current !== null, true, 'The current ligand expected.');
    expect(viewer.ligands.current!.rowIdx, 0, 'The current ligand of rowIdx = 0.');
    expect(viewer.ligands.selected.length, 0);

    df.selection.init((rowI) => rowI === 1 || rowI === 2);
    await delay(DebounceIntervals.ligands * 2.5); // await for debounce onRebuildViewLigands
    await Promise.all([awaitGrid(view.grid), viewer.awaitRendered()]);
    expect(viewer.ligands.current != null, true, 'The current ligand expected.');
    expect(viewer.ligands.current!.rowIdx, 0, 'The current ligand of rowIdx = 0.');
    expect(viewer.ligands.selected.length, 2);

    df.currentRowIdx = 3;
    await delay(DebounceIntervals.ligands * 2.5); // await for debounce onRebuildViewLigands
    await Promise.all([awaitGrid(view.grid), viewer.awaitRendered()]);
    expect(viewer.ligands.current != null, true, 'The current ligand expected.');
    expect(viewer.ligands.current!.rowIdx, 3, 'The current ligand of rowIdx = 3.');
    expect(viewer.ligands.selected.length, 2);
  });

  test('ligands-Molecule3D', async () => {
    const [pdbqtCnt, targetData]: [string, BiostructureData] = await Promise.all([
      _package.files.readAsText('docking/ligand_out.pdbqt'),
      (async () => {
        return {binary: true, ext: 'pdbqt', data: await _package.files.readAsBytes('docking/3SWZ.pdbqt')};
      })(),
    ]);
    const df: DG.DataFrame = (await grok.functions.call('BiostructureViewer:importPdbqt',
      {fileContent: pdbqtCnt, test: true}))[0];
    const ligandCol = df.getCol('molecule');

    const view = grok.shell.addTableView(df);
    await awaitGrid(view.grid);

    const viewer = (await df.plot.fromType('Biostructure', {
      [msvPROPS.dataJson]: BiostructureDataJson.fromData(targetData),
      [msvPROPS.ligandColumnName]: ligandCol.name,
      [msvPROPS.showCurrentRowLigand]: true,
      [msvPROPS.showSelectedRowsLigands]: true,
    })) as unknown as MolstarViewer;
    view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Biostructure', 0.4);

    df.currentRowIdx = 0;
    await delay(50); // await for debounce onRebuildViewLigands
    await Promise.all([awaitGrid(view.grid), viewer.awaitRendered()]);
    expect(viewer.ligands.current !== null, true, 'The current ligand expected');
    expect(viewer.ligands.current!.rowIdx, 0, 'The current ligand of rowIdx = 0.');
    expect(viewer.ligands.selected.length, 0);

    df.selection.init((rowI) => rowI === 2 || rowI === 3);
    await delay(50); // await for debounce onRebuildViewLigands
    await Promise.all([awaitGrid(view.grid), viewer.awaitRendered()]);
    expect(viewer.ligands.current != null, true, 'The current ligand expected');
    expect(viewer.ligands.current!.rowIdx, 0, 'The current ligand of rowIdx = 0.');
    expect(viewer.ligands.selected.length, 2);

    df.currentRowIdx = 1;
    await delay(50); // await for debounce onRebuildViewLigands
    await Promise.all([awaitGrid(view.grid), viewer.awaitRendered()]);
    expect(viewer.ligands.current != null, true, 'The current ligand expected');
    expect(viewer.ligands.current!.rowIdx, 1, 'The current ligand of rowIdx = 1.');
    expect(viewer.ligands.selected.length, 2);
  });
});
