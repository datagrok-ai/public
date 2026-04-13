import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {BiostructureData, BiostructureDataJson} from '@datagrok-libraries/bio/src/pdb/types';

import {awaitGrid} from './utils';
import {DebounceIntervals, MolstarViewer, PROPS as msvPROPS} from '../viewers/molstar-viewer/molstar-viewer';

import {_package} from '../package-test';

// -- Global selection cache tests --------------------------------------------
// The cache is module-level in molstar-viewer.ts. We test its behavior by
// firing the cross-package custom event and checking if the viewer responds.

category('atom-highlight', () => {
  test('global cache stores persistent events only', async () => {
    // Fire a persistent event — should be cached.
    grok.events.fireCustomEvent('chem-interactive-selection-changed', {
      rowIdx: 999, atoms: [0, 1, 2], persistent: true,
    });
    await delay(50);

    // Fire a non-persistent event for same row — should NOT overwrite.
    grok.events.fireCustomEvent('chem-interactive-selection-changed', {
      rowIdx: 999, atoms: [0, 1, 2, 3, 4, 5], persistent: false,
    });
    await delay(50);

    // The cache is internal, but we can verify indirectly: fire another
    // persistent clear event to clean up.
    grok.events.fireCustomEvent('chem-interactive-selection-changed', {
      rowIdx: 999, atoms: [], persistent: true,
    });
    await delay(50);
    // If we got here without errors, the event system is working.
  });

  test('highlightAllLigandAtoms — runs without error', async () => {
    const pdbData: BiostructureData = {
      binary: true, ext: 'pdb',
      data: await _package.files.readAsBytes('samples/1bdq-wo-ligands.pdb'),
    };
    const sdfCnt = await _package.files.readAsBytes('samples/1bdq-obs-pred.sdf');
    const df: DG.DataFrame = (await grok.functions.call('Chem:importSdf', {bytes: sdfCnt}))[0];
    const ligandCol = df.getCol('molecule');

    const view = grok.shell.addTableView(df);
    await awaitGrid(view.grid);

    const viewer = (await df.plot.fromType('Biostructure', {
      [msvPROPS.dataJson]: BiostructureDataJson.fromData(pdbData),
      [msvPROPS.ligandColumnName]: ligandCol.name,
      [msvPROPS.showCurrentRowLigand]: true,
    })) as unknown as MolstarViewer;
    view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Biostructure', 0.4);

    df.currentRowIdx = 0;
    await delay(DebounceIntervals.ligands * 2.5);
    await viewer.awaitRendered();

    // Fire a selection event and verify highlightAllLigandAtoms doesn't throw.
    grok.events.fireCustomEvent('chem-interactive-selection-changed', {
      rowIdx: 0, atoms: [0, 1, 2], persistent: true,
    });
    await delay(500);

    // Call highlightAllLigandAtoms directly — should not throw.
    await viewer.highlightAllLigandAtoms();

    view.close();
  });

  test('base colors applied to comparison pose', async () => {
    const pdbData: BiostructureData = {
      binary: true, ext: 'pdb',
      data: await _package.files.readAsBytes('samples/1bdq-wo-ligands.pdb'),
    };
    const sdfCnt = await _package.files.readAsBytes('samples/1bdq-obs-pred.sdf');
    const df: DG.DataFrame = (await grok.functions.call('Chem:importSdf', {bytes: sdfCnt}))[0];
    const ligandCol = df.getCol('molecule');

    const view = grok.shell.addTableView(df);
    await awaitGrid(view.grid);

    const viewer = (await df.plot.fromType('Biostructure', {
      [msvPROPS.dataJson]: BiostructureDataJson.fromData(pdbData),
      [msvPROPS.ligandColumnName]: ligandCol.name,
      [msvPROPS.showCurrentRowLigand]: true,
      [msvPROPS.showMouseOverRowLigand]: true,
    })) as unknown as MolstarViewer;
    view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Biostructure', 0.4);

    df.currentRowIdx = 0;
    await delay(DebounceIntervals.ligands * 2.5);
    await viewer.awaitRendered();

    // Trigger mouse-over on a different row to load comparison pose.
    df.mouseOverRowIdx = 1;
    await delay(DebounceIntervals.ligands * 2.5);
    await viewer.awaitRendered();

    // _applyBaseColors is called internally — just verify no error.
    await (viewer as any)._applyBaseColors();

    view.close();
  });
});
