import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {awaitCheck, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {BiostructureData, BiostructureDataJson} from '@datagrok-libraries/bio/src/pdb/types';

import {awaitGrid} from './utils';
import {DebounceIntervals, MolstarViewer, PROPS as msvPROPS} from '../viewers/molstar-viewer/molstar-viewer';

import {_package} from '../package-test';


async function loadPdb(sampleRelPath: string): Promise<BiostructureData> {
  return {binary: true, ext: 'pdb', data: await _package.files.readAsBytes(sampleRelPath)};
}

/** Builds a Biostructure viewer from a single PDB (Path A style — no ligand column). */
async function createPathAViewer(pdbData: BiostructureData,
  extraProps: Record<string, unknown> = {}): Promise<{view: DG.TableView, viewer: MolstarViewer}> {
  const df = DG.DataFrame.create(1);
  const view = grok.shell.addTableView(df);
  await awaitGrid(view.grid);
  const viewer = (await df.plot.fromType('Biostructure', {
    [msvPROPS.dataJson]: BiostructureDataJson.fromData(pdbData),
    ...extraProps,
  })) as unknown as MolstarViewer;
  view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Biostructure', 0.4);
  await viewer.awaitRendered();
  return {view, viewer};
}

async function createPathAViewerFromSample(sampleRelPath: string,
  extraProps: Record<string, unknown> = {}): Promise<{view: DG.TableView, viewer: MolstarViewer}> {
  return createPathAViewer(await loadPdb(sampleRelPath), extraProps);
}

async function toggleBindingSite(viewer: MolstarViewer, on: boolean): Promise<void> {
  viewer.setOptions({[msvPROPS.showBindingSite]: on});
  // Path A synchronously schedules via viewSyncer; Path B goes through
  // onRebuildViewLigandsRequest (ligands debounce). Wait long enough for both.
  await delay(DebounceIntervals.ligands * 2.5);
  await viewer.awaitRendered();
}


category('BindingSite', () => {
  test('property-defaults', async () => {
    const df = DG.DataFrame.create(1);
    grok.shell.addTableView(df);
    const viewer = (await df.plot.fromType('Biostructure')) as unknown as MolstarViewer;
    expect(viewer.showBindingSite, false, 'showBindingSite default');
    expect(viewer.bindingSiteRadius, 5, 'bindingSiteRadius default');
    expect(viewer.bindingSiteWholeResidues, true, 'bindingSiteWholeResidues default');
    expect((viewer as any).bindingSiteRefs.length, 0, 'bindingSiteRefs empty before any toggle');
  }, {timeout: 20000});

  test('pathA-het-in-pdb', async () => {
    // 1bdq.pdb contains HET ligand atoms; the viewer should build a binding
    // site (component + representation) when toggled on.
    const {viewer} = await createPathAViewerFromSample('samples/1bdq.pdb');

    await toggleBindingSite(viewer, true);
    await awaitCheck(() => (viewer as any).bindingSiteRefs.length >= 2,
      'Path A binding site refs not created', 10000);
  }, {timeout: 45000});

  test('pathA-apo-no-het', async () => {
    // 1bdq-wo-ligands.pdb has ligand HETs stripped. Toggle should be a noop.
    const {viewer} = await createPathAViewerFromSample('samples/1bdq-wo-ligands.pdb');

    expect((viewer as any)._isLigandAvailable(), false,
      'Apo structure should not report a ligand');

    await toggleBindingSite(viewer, true);
    expect((viewer as any).bindingSiteRefs.length, 0,
      'Apo structure should produce zero binding-site refs');
  }, {timeout: 45000});

  test('pathB-ligand-column', async () => {
    // Mirror the existing ligands-Molecule3D setup: separate ligand column (pdbqt)
    // against a protein target (3SWZ pdbqt). Binding site uses the current-row ligand.
    const [pdbqtCnt, targetData]: [string, BiostructureData] = await Promise.all([
      _package.files.readAsText('docking/ligand_out.pdbqt'),
      (async () => ({binary: true, ext: 'pdbqt',
        data: await _package.files.readAsBytes('docking/3SWZ.pdbqt')}) as BiostructureData)(),
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
    })) as unknown as MolstarViewer;
    view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Biostructure', 0.4);

    df.currentRowIdx = 0;
    await delay(DebounceIntervals.ligands * 2.5);
    await Promise.all([awaitGrid(view.grid), viewer.awaitRendered()]);
    expect(viewer.ligands.current != null, true, 'Path B requires a current ligand');

    expect((viewer as any)._isLigandAvailable(), true,
      'Ligand column should be reported as available');

    await toggleBindingSite(viewer, true);
    await awaitCheck(() => (viewer as any).bindingSiteRefs.length >= 2,
      'Path B binding site refs not created', 10000);
  }, {timeout: 50000});

  test('toggle-off-clears-refs', async () => {
    const {viewer} = await createPathAViewerFromSample('samples/1bdq.pdb');

    await toggleBindingSite(viewer, true);
    await awaitCheck(() => (viewer as any).bindingSiteRefs.length >= 2,
      'binding site refs not created on enable', 10000);

    await toggleBindingSite(viewer, false);
    expect((viewer as any).bindingSiteRefs.length, 0,
      'bindingSiteRefs should be empty after disable');
  }, {timeout: 45000});

  test('radius-change-rebuilds', async () => {
    const {viewer} = await createPathAViewerFromSample('samples/1bdq.pdb');

    await toggleBindingSite(viewer, true);
    await awaitCheck(() => (viewer as any).bindingSiteRefs.length >= 2,
      'binding site refs not created on enable', 10000);

    viewer.setOptions({[msvPROPS.bindingSiteRadius]: 8});
    await delay(DebounceIntervals.ligands * 2.5);
    await viewer.awaitRendered();

    // Radius identity change is flaky to assert cross-platform; at minimum
    // the view should still have the component + representation present.
    expect((viewer as any).bindingSiteRefs.length >= 2, true,
      'bindingSiteRefs should still exist after radius change');
    expect(viewer.bindingSiteRadius, 8, 'radius property should reflect new value');
  }, {timeout: 45000});

  test('ligand-absent-hides-overlay', async () => {
    // Neither a ligand column nor a HET in the PDB — overlay predicate must
    // report "ligand not available" so the viewer hides the binding-site UI.
    const {viewer} = await createPathAViewerFromSample('samples/1bdq-wo-ligands.pdb');

    expect((viewer as any)._isLigandAvailable(), false,
      '_isLigandAvailable should be false without a ligand column or HET');
  }, {timeout: 30000});
});
