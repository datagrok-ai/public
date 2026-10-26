/* eslint-disable max-len, max-lines-per-function */
/**
 * The "Chemical engine" setting (src/engine/README.md): RDKit mode never fetches the Crux file; switching both ways and
 * while the file loads; a failed load (simulated, and a missing file) and a failed worker job fall back; the renderer,
 * `canvasMol`, `drawMolecule` and the SVG export in both engines, the renderer's caches per engine, the Crux heap over
 * thousands of renders, the evicted molecules freed; descriptors and the scaffold tree from Crux against the chem-chem
 * container and the Python script; the substructure search's Crux engine against RDKit; fingerprints, their caches,
 * structural alerts and pharmacophore features from Crux against the RDKit pool. The tests of the switch itself call
 * `setChemicalEngine`; the others switch through ./engine-utils.
 */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {after, awaitCheck, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {ensureContainerRunning} from '@datagrok-libraries/test/src/test-container-utils';
import {_package} from '../package-test';
import {PackageFunctions, _package as chemPackage} from '../package';
import {RDKitCellRenderer} from '../rendering/rdkit-cell-renderer';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {ChemicalEngine, chemSetting, cruxReady, engineInEffect, isCruxRequested, onEngineSwitch, requestedEngine,
  setChemicalEngine} from '../engine/chem-engine';
import {CRUX_WASM_FILE, cruxHeapBytes, cruxLoadFailed, cruxLoadRequests, cruxWasmUrl, holdCruxForTests, noteCruxError,
  resetCruxLoaderForTests, setCruxFetchForTests, setCruxLoadFailureForTests, setCruxParseTrapForTests}
  from '../engine/crux-module';
import {calculateDescriptors, forgetCruxTreeForTests, getDescriptorsTree, lastDescriptorTiming,
  setCruxDescriptorGapForTests, setDescriptorTimingForTests} from '../descriptors/descriptors-engine';
import {cruxFallbackCount, cruxOrRdKitAsync} from '../engine/fallback';
import * as container from '../docker/api';
import {SubstructureSearchEngine, getCruxService, getSubstructureSearchEngine, isCruxSearch,
  resetCruxServiceForTests, setSubstructureSearchEngine} from '../crux/crux-searches';
import {Fingerprint} from '../utils/chem-common';
import {FILTER_TYPES, chemGetDiversities, chemGetFingerprints, chemGetSimilarities, chemSubstructureSearchLibrary}
  from '../chem-searches';
import {alertsByEngine, fingerprintsByEngine, inchisByEngine, structuralAlertsByEngine} from '../engine/worker-jobs';
import {RuleSet, runStructuralAlertsDetection} from '../panels/structural-alerts';
import {PHARM_FAMILY_NAMES, PharmFamilySet, runPharmacophoreDetection} from '../panels/pharmacophore-features';
import {CONTAINER_TIMEOUT, readDataframe} from './utils';
import {expectNoFallback, inBothEngines, inEngine, withCruxUnreadable} from './engine-utils';
import {isNoneValue} from '../crux/crux-api';
import {ALIGN_BY_SCAFFOLD_LAYOUT_PERSISTED_TAG, SubstructureSearchType} from '../constants';
import {calculatorMethodInfo} from '../widgets/biochem-properties-widget';
import {convertNotationForColumn} from '../utils/convert-notation-utils';
import {cruxInchi} from '../crux/crux-inchi';
import {getIdentifiersSingle} from '../widgets/identifiers';
import {getStructuralAlerts} from '../widgets/structural-alerts';
import {getPharmacophoreFeatures} from '../widgets/pharmacophore-features';
import {mmpLinkFragmentsByEngine} from '../engine/analysis-jobs';
import {getInverseSubstructuresAndAlign} from '../analysis/molecular-matched-pairs/mmp-viewer/mmp-mol-rendering';
import {MAX_SMILES_LENGTH} from '../utils/chem-constants';
import {TOXICITY_ALERTS_CRUX_FAILED, getToxicityAlertsColumns} from '../widgets/toxicity';

/** The renderer's internals the tests read (the molecule's engine, a raster). */
interface RendererInternals {
  _fetchMol(molString: string, scaffolds: {molecule: string, color?: string}[], molRegenerateCoords: boolean,
    scaffoldRegenerateCoords: boolean, details: object, alignByFirstSubstructure: boolean): {engine: ChemicalEngine,
    molCtx: {mol: any}};
  _fetchRender(width: number, height: number, molString: string, scaffolds: {molecule: string, color?: string}[],
    molRegenerateCoords: boolean, scaffoldRegenerateCoords: boolean, alignByFirstSubstructure: boolean,
    details?: object): ImageData;
  cachedMoleculeCount: number;
  _molToSvg(molString: string, w: number, h: number): string | null;
}

function newRenderer(): RendererInternals {
  return new RDKitCellRenderer(chemCommonRdKit.getRdKitModule()) as unknown as RendererInternals;
}

/** Pixels with any ink (the renderer draws on a transparent background). */
function inked(image: ImageData): number {
  let n = 0;
  for (let i = 3; i < image.data.length; i += 4)
    n += image.data[i] > 0 ? 1 : 0;
  return n;
}

function sameImage(a: ImageData, b: ImageData): boolean {
  return a.width === b.width && a.height === b.height && a.data.every((v, i) => v === b.data[i]);
}

function molblockOf(smiles: string, v3000 = false): string {
  const mol = chemCommonRdKit.getRdKitModule().get_mol(smiles);
  try {
    mol.set_new_coords();
    return v3000 ? mol.get_v3Kmolblock() : mol.get_molblock();
  } finally {
    mol.delete();
  }
}

const ASPIRIN = 'CC(=O)Oc1ccccc1C(=O)O';
/** AZT monophosphate as the Crux sketcher wrote it (V3000, an MDLV30/STEABS collection): a reported bug, Crux's
 * descriptors read it as nothing (a row of nulls; RDKit and the container read it). */
const AZT_MP_V3000 = `
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 23 24 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 6.615700 -0.062500 0.000000 0
M  V30 2 C 6.615700 -0.887500 0.000000 0
M  V30 3 N 5.903700 -1.295800 0.000000 0
M  V30 4 C 5.191700 -0.887500 0.000000 0
M  V30 5 N 5.191700 -0.062500 0.000000 0
M  V30 6 C 5.903700 0.354200 0.000000 0
M  V30 7 O 5.911500 1.179200 0.000000 0
M  V30 8 O 4.482400 -1.300000 0.000000 0
M  V30 9 C 7.328200 0.354200 0.000000 0
M  V30 10 C 5.908300 -2.779200 0.000000 0
M  V30 11 C 4.791600 -3.495900 0.000000 0
M  V30 12 C 5.537500 -3.495900 0.000000 0
M  V30 13 O 5.200000 -2.212500 0.000000 0
M  V30 14 C 4.422400 -2.782600 0.000000 0
M  V30 15 C 4.004200 -2.062500 0.000000 0
M  V30 16 O 3.179200 -2.058300 0.000000 0
M  V30 17 P 2.354200 -2.058300 0.000000 0
M  V30 18 O 1.529200 -2.058300 0.000000 0
M  V30 19 O 2.350000 -1.233300 0.000000 0
M  V30 20 O 2.350000 -2.883300 0.000000 0
M  V30 21 N 4.200000 -4.079200 0.000000 0
M  V30 22 N 3.495800 -4.087500 0.000000 0 CHG=1
M  V30 23 N 2.779200 -4.087500 0.000000 0 CHG=-1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 4 5
M  V30 2 1 5 6
M  V30 3 1 12 10
M  V30 4 1 10 13
M  V30 5 1 13 14
M  V30 6 1 14 11
M  V30 7 2 6 7
M  V30 8 1 14 15 CFG=1
M  V30 9 2 1 2
M  V30 10 1 15 16
M  V30 11 2 4 8
M  V30 12 1 16 17
M  V30 13 1 1 6
M  V30 14 1 17 18
M  V30 15 1 1 9
M  V30 16 2 17 19
M  V30 17 1 2 3
M  V30 18 1 17 20
M  V30 19 1 10 3 CFG=1
M  V30 20 1 11 21 CFG=3
M  V30 21 1 11 12
M  V30 22 2 21 22
M  V30 23 1 3 4
M  V30 24 2 22 23
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(3 10 11 14)
M  V30 END COLLECTION
M  V30 END CTAB
M  END`;
/** 2D descriptors the descriptor tests compare (one of each kind: masses, counts, surface, logP, rings). */
const DESCRIPTORS_2D = ['MolWt', 'HeavyAtomCount', 'NumHDonors', 'NumHAcceptors', 'TPSA', 'MolLogP', 'NumRotatableBonds',
  'RingCount'];

/** Where `actual`'s values differ from `expected`'s (beyond 1e-6 relative; text values read as numbers), or are null. */
function descriptorDifferences(expected: DG.Column[], actual: DG.Column[]): string[] {
  const differences: string[] = [];
  const number = (c: DG.Column, i: number) => {
    const v = c.isNone(i) ? null : c.get(i);
    return v === null || v === '' ? null : Number(v);
  };
  for (const e of expected) {
    const a = actual.find((c) => c.name === e.name);
    for (let i = 0; i < e.length; i++) {
      const [x, y] = [number(e, i), a ? number(a, i) : null];
      if (y === null || x === null || Number.isNaN(y) || Math.abs(x - y) > 1e-6 * Math.max(1, Math.abs(x)))
        differences.push(`${e.name} row ${i}: ${y} (the container's ${x})`);
    }
  }
  return differences;
}
/** A ~600-atom peptide-like chain: a very large molecule for the grid. */
const LARGE = 'N' + 'C(C)C(=O)N'.repeat(120) + 'CC(=O)O';

category('chemical engine', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  after(async () => {
    setCruxLoadFailureForTests(null);
    setChemicalEngine(null);
    setSubstructureSearchEngine(null);
  });

  test('the settings reach the test bundle', async () => {
    let own: string;
    try {
      own = JSON.stringify(chemPackage.settings);
    } catch (e: any) {
      own = `throws ${e?.message ?? e}`;
    }
    // package.ts's _package in this bundle is not the published package: chemSetting asks the published one
    console.log(`chemical engine | in package-test.js, package.ts's _package "${chemPackage.name}" settings: ${own}; ` +
      `chemSetting: ChemicalEngine ${chemSetting('ChemicalEngine')}, Sketcher ${chemSetting('Sketcher')}`);
    expect(chemSetting('ChemicalEngine') !== undefined || chemSetting('Sketcher') !== undefined, true,
      'the settings of the published Chem are read');
  });

  test('the setting of the stand reaches the renderers', async () => {
    // a run with the stand's Chemical Engine set to Crux (misc/crux/engine-setting.sh) draws with Crux here too
    setChemicalEngine(null);
    const engine = requestedEngine();
    if (engine === ChemicalEngine.Crux)
      expect(await cruxReady(), true, 'the Crux module loads');
    const drawnBy = newRenderer()._fetchMol(ASPIRIN, [], false, false, {}, false).engine;
    console.log(`chemical engine | the stand's setting: ${engine}; aspirin drawn by ${drawnBy}`);
    expect(drawnBy, engine);
  });

  test('RDKit mode does not ask for the Crux file', async () => {
    await inEngine(ChemicalEngine.RDKit, async () => {
      const before = cruxLoadRequests().length;
      const renderer = newRenderer();
      for (const molecule of [ASPIRIN, molblockOf(ASPIRIN), 'c1ccccc1[N,O]', LARGE])
        renderer._fetchRender(200, 100, molecule, [{molecule: 'c1ccccc1'}], false, false, false);
      const canvas = document.createElement('canvas');
      chemCommonRdKit.drawMoleculeToCanvas(0, 0, 200, 100, canvas, ASPIRIN, 'c1ccccc1');
      PackageFunctions.drawMolecule(ASPIRIN, 200, 100, false);
      expect(isCruxRequested(), false, 'descriptors and the scaffold tree go to the container and the script');
      expect(engineInEffect(), ChemicalEngine.RDKit);
      expect(cruxLoadRequests().slice(before).join(', '), '');
    });
  });

  test('switching both ways', async () => {
    setChemicalEngine(ChemicalEngine.Crux);
    expect(await cruxReady(), true, 'the Crux module loads');
    const renderer = newRenderer();
    expect(renderer._fetchMol(ASPIRIN, [], false, false, {}, false).engine, ChemicalEngine.Crux);
    setChemicalEngine(ChemicalEngine.RDKit);
    expect(renderer._fetchMol(ASPIRIN, [], false, false, {}, false).engine, ChemicalEngine.RDKit);
    setChemicalEngine(ChemicalEngine.Crux);
    expect(renderer._fetchMol(ASPIRIN, [], false, false, {}, false).engine, ChemicalEngine.Crux);
  });

  test('a switch draws with RDKit until Crux is ready, then repaints', async () => {
    setChemicalEngine(ChemicalEngine.RDKit);
    const renderer = newRenderer();
    let repaints = 0;
    const stop = onEngineSwitch(() => repaints++);
    const release = holdCruxForTests();
    try {
      setChemicalEngine(ChemicalEngine.Crux);
      expect(isCruxRequested(), true);
      expect(renderer._fetchMol(ASPIRIN, [], false, false, {}, false).engine, ChemicalEngine.RDKit,
        'while the file loads');
      expect(repaints, 0);
      release();
      expect(await cruxReady(), true);
      expect(repaints > 0, true, 'the switch listeners (the grids\' repaint) ran once Crux was ready');
      expect(renderer._fetchMol(ASPIRIN, [], false, false, {}, false).engine, ChemicalEngine.Crux);
    } finally {
      release();
      stop();
    }
  });

  test('a load failure falls back to RDKit', async () => {
    setChemicalEngine(ChemicalEngine.Crux);
    setCruxLoadFailureForTests(new Error('simulated: chem_wasm_bg.wasm 404'));
    try {
      expect(isCruxRequested(), false, 'Crux mode with a failed load');
      expect(engineInEffect(), ChemicalEngine.RDKit);
      const renderer = newRenderer();
      expect(renderer._fetchMol(ASPIRIN, [], false, false, {}, false).engine, ChemicalEngine.RDKit);
      expect(inked(renderer._fetchRender(200, 100, ASPIRIN, [], false, false, false)) > 0, true);
      expect(await cruxReady(), false);
    } finally {
      setCruxLoadFailureForTests(null);
    }
    expect(await cruxReady(), true, 'Crux back once the simulated failure is cleared');
  });

  test('a missing Crux file falls back to RDKit', async () => {
    // the real loader (fetch, compile) on a URL with no file behind it
    setChemicalEngine(ChemicalEngine.Crux);
    resetCruxLoaderForTests(cruxWasmUrl().replace(CRUX_WASM_FILE, 'no_such_chem_wasm_bg.wasm'));
    try {
      expect(await cruxReady(), false, 'the load fails');
      expect(isCruxRequested(), false, 'Crux mode with a failed load');
      expect(engineInEffect(), ChemicalEngine.RDKit);
      const renderer = newRenderer();
      expect(renderer._fetchMol(ASPIRIN, [], false, false, {}, false).engine, ChemicalEngine.RDKit);
      expect(inked(renderer._fetchRender(200, 100, ASPIRIN, [], false, false, false)) > 0, true);
      const canvas = document.createElement('canvas');
      chemCommonRdKit.drawMoleculeToCanvas(0, 0, 200, 100, canvas, ASPIRIN, null);
      expect(inked(canvas.getContext('2d')!.getImageData(0, 0, canvas.width, canvas.height)) > 0, true,
        'canvasMol draws with RDKit');
    } finally {
      resetCruxLoaderForTests(null);
    }
    expect(await cruxReady(), true, 'the real file loads after the reset');
    expect(newRenderer()._fetchMol(ASPIRIN, [], false, false, {}, false).engine, ChemicalEngine.Crux);
  });

  test('a stalled Crux file falls back to RDKit and the container', async () => {
    // a fetch that never settles: past the load's limit (1 s here, 60 s in a session) the descriptors tree comes from
    // the container and the search runs on RDKit, with RDKit's hits
    await ensureContainerRunning('name = "chem-chem"', CONTAINER_TIMEOUT);
    const col = (await readDataframe('smiles.csv')).col('canonical_smiles')!;
    col.semType = DG.SEMTYPE.MOLECULE;
    setChemicalEngine(ChemicalEngine.Crux);
    setSubstructureSearchEngine(SubstructureSearchEngine.RDKit);
    const rdkitHits: DG.BitSet = (await PackageFunctions.searchSubstructure(col, 'c1ccccc1', '')).get(0);
    const containerTree = await container.getDescriptorsTree();
    await resetCruxServiceForTests();
    forgetCruxTreeForTests();
    resetCruxLoaderForTests(null);
    setCruxFetchForTests(() => new Promise<Response>(() => {}), 1000);
    const started = performance.now();
    try {
      const tree = await getDescriptorsTree();
      expect(cruxLoadFailed(), true, 'the stalled load has failed');
      expect(JSON.stringify(tree), JSON.stringify(containerTree), 'the descriptors tree: the container\'s');
      resetCruxLoaderForTests(null); // the search meets the stalled fetch too
      setSubstructureSearchEngine(SubstructureSearchEngine.Crux);
      const hits: DG.BitSet = (await PackageFunctions.searchSubstructure(col, 'c1ccccc1', '')).get(0);
      expect(hits.toString(), rdkitHits.toString(), 'the search: RDKit\'s hits');
      expect(performance.now() - started < 20000, true, 'within the limits, not hung');
    } finally {
      setCruxFetchForTests(null);
      resetCruxLoaderForTests(null);
      setSubstructureSearchEngine(null);
    }
    expect(await cruxReady(), true, 'Crux back after the reset');
  }, {timeout: 60000 + CONTAINER_TIMEOUT});

  test('a trap in a Crux call leaves RDKit for the session', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      try {
        noteCruxError(new Error('simulated: a molecule Crux refuses'));
        expect(isCruxRequested(), true, 'an error other than a trap leaves Crux on');
        noteCruxError(new WebAssembly.RuntimeError('simulated: unreachable'));
        expect(isCruxRequested(), false, 'after a trap: RDKit');
        expect(engineInEffect(), ChemicalEngine.RDKit);
        expect(newRenderer()._fetchMol(ASPIRIN, [], false, false, {}, false).engine, ChemicalEngine.RDKit);
      } finally {
        resetCruxLoaderForTests(null); // the next session's Crux
      }
      expect(await cruxReady(), true, 'Crux back after the reset');
    });
  });

  test('a trap while Crux parses a molecule leaves RDKit for the session', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      setCruxParseTrapForTests(true);
      try {
        const info = newRenderer()._fetchMol(ASPIRIN, [], false, false, {}, false);
        expect(isCruxRequested(), false, 'the trap in get_mol marked Crux failed');
        expect(info.engine, ChemicalEngine.RDKit, 'RDKit drew the molecule');
      } finally {
        setCruxParseTrapForTests(false);
        resetCruxLoaderForTests(null); // the next session's Crux
      }
      expect(await cruxReady(), true, 'Crux back after the reset');
    });
  });

  test('a repeatable Crux job error keeps the workers', async () => {
    // only a trap, a timeout or a lost worker restarts them (and drops their search indexes)
    await inEngine(ChemicalEngine.Crux, async () => {
      const service = await getCruxService();
      const generation = service.generation;
      let failure = '';
      try {
        await service.failJobForTests('simulated: no scaffold tree for these molecules');
      } catch (e: any) {
        failure = String(e?.message ?? e);
      }
      expect(failure.includes('simulated'), true, 'the job failed');
      expect(service.generation, generation, 'the workers were kept');
    });
  });

  test('Crux jobs: each chunk is timed from its own start', async () => {
    // three chunks per worker, each keeping its worker busy 400 ms, within 1,000 ms each: all done, none timed out
    await inEngine(ChemicalEngine.Crux, async () => {
      const service = await getCruxService();
      const generation = service.generation;
      const chunks = service.workerCount * 3;
      expect(await service.busyChunksForTests(chunks, 400, 1000), chunks);
      expect(service.generation, generation, 'no chunk timed out (no restart)');
    });
  }, {timeout: 60000});

  test('a failed Crux worker job falls back to RDKit and restarts the workers', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      const service = await getCruxService();
      const generation = service.generation;
      const fallbacks = cruxFallbackCount('test job');
      const answer = await cruxOrRdKitAsync('test job', () => service.failJobForTests(), async () => 'RDKit');
      expect(answer, 'RDKit', 'the RDKit path answered');
      expect(cruxFallbackCount('test job'), fallbacks + 1, 'the fallback is counted');
      expect(service.generation, generation + 1, 'the workers restarted');
      expect(isCruxRequested(), true, 'a worker\'s failure leaves Crux mode on');
      const tree = JSON.parse(await service.scaffoldTree([ASPIRIN], 10, false));
      expect(Array.isArray(tree) && tree.length > 0, true, 'the new workers compute');
    });
  });

  test('renders in both engines', async () => {
    const benzeneScaffold = molblockOf('c1ccccc1');
    type Case = {name: string, molecule: string, scaffolds?: {molecule: string, color?: string}[], align?: boolean,
      crux: boolean, ink: boolean};
    const cases: Case[] = [
      {name: 'empty', molecule: '', crux: true, ink: false},
      {name: 'invalid', molecule: 'C1CC(', crux: false, ink: true}, // RDKit's error cross
      {name: 'SMILES', molecule: ASPIRIN, crux: true, ink: true},
      {name: 'molblock V2000', molecule: molblockOf(ASPIRIN), crux: true, ink: true},
      {name: 'molblock V3000', molecule: molblockOf(ASPIRIN, true), crux: true, ink: true},
      {name: 'SMARTS', molecule: 'c1ccccc1[N,O;!H0]', crux: true, ink: true},
      {name: 'substructure highlight', molecule: ASPIRIN, scaffolds: [{molecule: 'c1ccccc1', color: '#ff0000'}],
        crux: true, ink: true},
      {name: 'aligned to a scaffold', molecule: ASPIRIN, scaffolds: [{molecule: benzeneScaffold}], align: true,
        crux: true, ink: true},
      {name: 'very large', molecule: LARGE, crux: true, ink: true},
      {name: 'longer than Chem reads', molecule: 'C'.repeat(6000), crux: false, ink: true},
    ];
    const results: string[] = [];
    for (const c of cases) {
      const ink = await inBothEngines(async (engine) => {
        const renderer = newRenderer();
        const info = renderer._fetchMol(c.molecule, c.scaffolds ?? [], false, false, {}, c.align ?? false);
        const drawnBy = engine === ChemicalEngine.Crux && c.crux ? ChemicalEngine.Crux : ChemicalEngine.RDKit;
        expect(info.engine, drawnBy, `${c.name}: drawn by`);
        const pixels = inked(renderer._fetchRender(200, 100, c.molecule, c.scaffolds ?? [], false, false, c.align ?? false));
        expect(pixels > 0, c.ink, `${c.name} in ${engine}: inked`);
        return pixels;
      });
      // the same picture to within a third of its ink (Crux draws as RDKit 2026.03, Chem's RDKit is 2024.09)
      if (c.ink)
        expect(Math.abs(ink.crux - ink.rdkit) <= ink.rdkit / 3, true, `${c.name}: ink ${ink.rdkit} RDKit, ${ink.crux} Crux`);
      results.push(`${c.name}: ${ink.rdkit}/${ink.crux}`);
    }
    console.log(`chemical engine | inked pixels RDKit/Crux: ${results.join('; ')}`);
  });

  test('canvasMol, drawMolecule and the SVG export draw with Crux', async () => {
    await inEngine(ChemicalEngine.Crux, () => expectNoFallback(['canvasMol', 'rendering', 'drawing'], async () => {
      const canvas = document.createElement('canvas');
      chemCommonRdKit.drawMoleculeToCanvas(0, 0, 200, 100, canvas, ASPIRIN, 'c1ccccc1');
      expect(inked(canvas.getContext('2d')!.getImageData(0, 0, canvas.width, canvas.height)) > 0, true,
        'canvasMol inked');
      const host = PackageFunctions.drawMolecule(ASPIRIN, 200, 100, false);
      const drawn = (host instanceof HTMLCanvasElement ? host : host.querySelector('canvas'))!;
      await awaitCheck(() => inked(drawn.getContext('2d')!.getImageData(0, 0, drawn.width, drawn.height)) > 0,
        'drawMolecule did not ink its canvas', 10000);
      const renderer = newRenderer();
      expect(renderer._fetchMol(ASPIRIN, [], false, false, {}, false).engine, ChemicalEngine.Crux);
      const svg = renderer._molToSvg(ASPIRIN, 200, 100) ?? '';
      expect(svg.includes('<svg') && svg.includes('</svg>'), true, 'an SVG');
    }));
  });

  test('a scaffold on the column realigns the Crux pictures, and they come back without it', async () => {
    await inEngine(ChemicalEngine.Crux, () => expectNoFallback(['rendering', 'drawing'], async () => {
      const df = await readDataframe('tests/smiles-50.csv');
      const col = df.col('canonical_smiles')!;
      col.semType = DG.SEMTYPE.MOLECULE;
      const grid = DG.Viewer.grid(df);
      const renderer = new RDKitCellRenderer(chemCommonRdKit.getRdKitModule());
      const r = window.devicePixelRatio;
      const canvas = document.createElement('canvas');
      canvas.width = 200 * r;
      canvas.height = 100 * r;
      const g = canvas.getContext('2d', {willReadFrequently: true})!;
      const picture = (row: number) => {
        g.clearRect(0, 0, canvas.width, canvas.height);
        const cell = grid.cell('canonical_smiles', row);
        renderer.render(g, 0, 0, 200, 100, cell, cell.style);
        return g.getImageData(0, 0, canvas.width, canvas.height);
      };
      /** Benzene as a sketcher draws it, its first atom at `deg` degrees. */
      const benzene = (deg: number) => {
        const f = (v: number) => v.toFixed(4).padStart(10);
        const atoms = Array.from({length: 6}, (_, i) => (deg + 60 * i) * Math.PI / 180)
          .map((a) => `${f(1.5 * Math.cos(a))}${f(1.5 * Math.sin(a))}    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0`);
        const bonds = Array.from({length: 6}, (_, i) =>
          `${String(i + 1).padStart(3)}${String((i + 1) % 6 + 1).padStart(3)}  ${i % 2 ? 1 : 2}  0`);
        return ['', '  test', '', '  6  6  0  0  0  0  0  0  0  0999 V2000', ...atoms, ...bonds, 'M  END'].join('\n');
      };
      const rows = [0, 1, 4, 5, 6];
      const before = rows.map(picture);
      col.setTag(ALIGN_BY_SCAFFOLD_LAYOUT_PERSISTED_TAG, benzene(0));
      const aligned = rows.map(picture);
      expect(rows.filter((_, i) => sameImage(aligned[i], before[i])).join(','), '', 'rows a flat benzene leaves as they were');
      // a structure whose benzene already lies as the scaffold's is not turned (RDKit 2026.03's alignment; HOST-063's
      // first cell with a pointed benzene)
      col.setTag(ALIGN_BY_SCAFFOLD_LAYOUT_PERSISTED_TAG, benzene(90));
      expect(sameImage(picture(0), before[0]), true, 'row 1 with a pointed benzene');
      expect(sameImage(picture(1), before[1]), false, 'row 2 with a pointed benzene');
      delete col.tags[ALIGN_BY_SCAFFOLD_LAYOUT_PERSISTED_TAG];
      const cleared = rows.map(picture);
      expect(rows.filter((_, i) => !sameImage(cleared[i], before[i])).join(','), '', 'rows not back without the scaffold');
    }));
  });

  test('caches are emptied on a switch', async () => {
    setChemicalEngine(ChemicalEngine.RDKit);
    const renderer = newRenderer();
    const smiles = (await readDataframe('smiles.csv')).col('canonical_smiles')!.toList().slice(0, 20) as string[];
    const rdkitInfos = smiles.map((s) => renderer._fetchMol(s, [], false, false, {}, false));
    expect(renderer.cachedMoleculeCount, 20);
    setChemicalEngine(ChemicalEngine.Crux);
    expect(await cruxReady(), true);
    const info = renderer._fetchMol(smiles[0], [], false, false, {}, false);
    expect(info.engine, ChemicalEngine.Crux);
    expect(renderer.cachedMoleculeCount, 1);
    expect(rdkitInfos.every((i) => i.molCtx.mol === null), true, 'the RDKit molecules are freed');
  });

  test('evicted molecules are freed and drawn again', async () => {
    // the cache holds 100 molecules: the first ones are evicted (and freed) by the time they are drawn again
    const smiles = (await readDataframe('smiles.csv')).col('canonical_smiles')!.toList().slice(0, 150) as string[];
    await inBothEngines(async (engine) => {
      const renderer = newRenderer();
      const first = smiles.slice(0, 20).map((s) => renderer._fetchMol(s, [], false, false, {}, false));
      for (const s of smiles)
        renderer._fetchRender(200, 100, s, [], false, false, false);
      expect(first.every((i) => i.molCtx.mol === null), true, `${engine}: evicted molecules freed`);
      expect(renderer.cachedMoleculeCount <= 100, true);
      const highlight = [{molecule: 'c1ccccc1', color: '#00ff00'}];
      const fresh = newRenderer();
      for (const s of smiles.slice(0, 20)) {
        const again = renderer._fetchRender(300, 150, s, highlight, false, false, false);
        expect(sameImage(again, fresh._fetchRender(300, 150, s, highlight, false, false, false)), true,
          `${engine}: ${s} drawn again as a fresh renderer draws it`);
      }
    });
  });

  test('the Crux heap does not grow over thousands of renders', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      const smiles = (await readDataframe('smiles.csv')).col('canonical_smiles')!.toList() as string[];
      const renderer = newRenderer();
      const pass = (w: number) => {
        for (const s of smiles)
          renderer._fetchRender(w, w / 2, s, [{molecule: 'c1ccccc1'}], false, false, false);
      };
      pass(200);
      const warm = cruxHeapBytes();
      pass(240);
      pass(280);
      pass(320);
      const end = cruxHeapBytes();
      console.log(`chemical engine | Crux heap ${warm} B after 1,000 renders, ${end} B after 4,000`);
      expect(end - warm <= 8 * 1024 * 1024, true, `the heap grew ${end - warm} B over 3,000 renders`);
    });
  }, {timeout: 300000});

  test('descriptors: Crux against the container', async () => {
    await ensureContainerRunning('name = "chem-chem"', CONTAINER_TIMEOUT);
    await inEngine(ChemicalEngine.Crux, () => expectNoFallback(['descriptorTree', 'calcDescriptors'], async () => {
      const containerTree = await container.getDescriptorsTree();
      const cruxTree = await getDescriptorsTree();
      // what the dialog builds from: the groups (key, name, version) and their descriptors' names, in order; the texts
      // (descriptions, tags) are reported: the stand's container image may run another chem.py than the package's
      const outline = (tree: any) => JSON.stringify(Object.entries(tree).map(([key, g]: [string, any]) =>
        [key, g.name, g.version, g.descriptors.map((d: any) => d.name)]));
      expect(outline(cruxTree), outline(containerTree), 'the tree: groups and names');
      const texts = (tree: any) => Object.values(tree).flatMap((g: any) => g.descriptors.map((d: any) =>
        `${d.name}: ${d.description} [${d.tags}]`));
      const [cruxTexts, containerTexts] = [texts(cruxTree), texts(containerTree)];
      console.log(`chemical engine | descriptor tree texts: ${cruxTexts.filter((t, i) => t !== containerTexts[i]).length} ` +
        `of ${cruxTexts.length} differ from the container's${JSON.stringify(cruxTree) === JSON.stringify(containerTree) ?
          ' (the JSON is identical)' : ''}`);
      const names: string[] = [];
      const names3D = new Set<string>();
      for (const tree of [containerTree, cruxTree]) {
        for (const group of Object.values(tree) as any[]) {
          for (const d of group.descriptors) {
            if (tree === containerTree)
              names.push(d.name);
            if (d.tags?.includes('3D'))
              names3D.add(d.name);
          }
        }
      }
      // the 2D names on 30 molecules (fewer than the 64 the container splits over processes; the 100,000-molecule
      // comparison, misc/crux/compare100k.mjs, is the large-scale check), the 3D ones (the container embeds each
      // molecule: slow; values reported only) on 2; an empty and an unreadable row with each
      const smiles = (await readDataframe('smiles.csv')).col('canonical_smiles')!.toList().slice(0, 20);
      const molblocks = (await readDataframe('tests/spgi-100.csv')).col('Structure')!.toList().slice(0, 10);
      const molecules = DG.Column.fromStrings('molecules', [...smiles, ...molblocks, '', 'C1CC(']);
      const few = DG.Column.fromStrings('molecules', [smiles[0], molblocks[0], '', 'C1CC(']);
      const [names2D, only3D] = [names.filter((n) => !names3D.has(n)), names.filter((n) => names3D.has(n))];
      const expected = [...await container.calculateDescriptors(molecules, names2D),
        ...only3D.length ? await container.calculateDescriptors(few, only3D) : []];
      const actual = [...await calculateDescriptors(molecules, names2D),
        ...only3D.length ? await calculateDescriptors(few, only3D) : []];
      expect(actual.map((c) => c.name).sort().join(','), expected.map((c) => c.name).sort().join(','), 'names');
      // 2D: within 1e-12 of the container's (relative, at least absolute 1e-12); 3D (the container embeds unseeded, so
      // its values depend on what its process embedded before; Crux embeds as a fresh process does) and Ipc / AvgIpc
      // (RDKit's own value is numpy summation-order noise on big molecules): reported only
      const TOLERANCE = 1e-12;
      const differences: string[] = [];
      const reported: string[] = [];
      // a container image from before GROK-20505 (commit 30215919ce) gives some descriptors as strings (a missing value
      // 'None', numbers printed to some 8 digits): their values are compared as numbers within 1e-6, the type reported
      const stringTyped: string[] = [];
      for (const e of expected) {
        const a = actual.find((c) => c.name === e.name)!;
        if (a.type !== e.type) {
          if (e.type === DG.COLUMN_TYPE.STRING)
            stringTyped.push(`${e.name} (Crux ${a.type})`);
          else
            differences.push(`${e.name}: type ${a.type}, the container's ${e.type}`);
        }
        const number = (v: any) => typeof v !== 'string' ? v : ['', 'None', 'nan', 'NaN'].includes(v) ? null : Number(v);
        const tolerance = e.type === DG.COLUMN_TYPE.STRING ? 1e-6 : TOLERANCE;
        let bad = 0;
        let maxRel = 0;
        let first = '';
        for (let i = 0; i < e.length; i++) {
          // a missing value: the int / double column's null, or the string column's '' (or 'None', 'nan')
          const x = e.isNone(i) ? null : number(e.get(i));
          const y = a.isNone(i) ? null : number(a.get(i));
          const same = x === null || y === null ? x === y :
            Math.abs(x - y) / Math.max(1, Math.abs(x)) <= tolerance;
          if (x !== null && y !== null)
            maxRel = Math.max(maxRel, Math.abs(x - y) / Math.max(1, Math.abs(x)));
          if (!same && bad++ === 0)
            first = `, row ${i} ${JSON.stringify(e.get(i))} / ${JSON.stringify(a.get(i))}`;
        }
        if (names3D.has(e.name) || e.name === 'Ipc' || e.name === 'AvgIpc')
          reported.push(`${e.name} ${bad}/${e.length} beyond ${TOLERANCE} (max ${maxRel.toPrecision(2)})`);
        else if (bad > 0)
          differences.push(`${e.name}: ${bad} of ${e.length} (max relative ${maxRel.toPrecision(2)}${first})`);
      }
      console.log(`chemical engine | ${names2D.length} 2D descriptors on ${molecules.length} molecules, ${only3D.length} ` +
        `3D ones on ${few.length}; reported only: ` +
        `${reported.join('; ')}; typed string by the container: ${stringTyped.join(', ') || 'none'}`);
      expect(differences.join('; '), '', `2D descriptors beyond ${TOLERANCE}`);
    }));
  }, {timeout: 600000 + CONTAINER_TIMEOUT});

  test('descriptors: Crux reads the Crux sketcher\'s V3000 molblocks', async () => {
    // a reported bug: the V3000 molblocks the Crux sketcher writes (an MDLV30/STEABS collection) came back as nulls
    await ensureContainerRunning('name = "chem-chem"', CONTAINER_TIMEOUT);
    await inEngine(ChemicalEngine.Crux, async () => {
      // a molecule with a stereocenter drawn in Crux and read back as its V3000 molblock
      const previous = grok.chem.currentSketcherType;
      grok.chem.currentSketcherType = 'Crux';
      const host = new DG.chem.Sketcher();
      host.setSmiles('N[C@@H](C)C(=O)O');
      const dialog = ui.dialog().add(host).show();
      let drawn = '';
      try {
        await host.sketcherReady();
        await delay(300);
        drawn = host.sketcher!.molV3000;
      } finally {
        dialog.close();
        grok.chem.currentSketcherType = previous;
      }
      expect(drawn.includes('V3000') && drawn.includes('CFG='), true, 'the sketcher\'s V3000 molblock with its stereocenter');
      const molecules = DG.Column.fromStrings('molecules', [AZT_MP_V3000, drawn]);
      // Crux read them: no row from the container
      const actual = await expectNoFallback(['calcDescriptors'], () => calculateDescriptors(molecules, DESCRIPTORS_2D));
      const molWt = actual.find((c) => c.name === 'MolWt')!.get(0);
      expect(Math.abs(molWt - 347.224) < 0.01, true, `AZT monophosphate's MolWt: ${molWt}`);
      const expected = await container.calculateDescriptors(molecules, DESCRIPTORS_2D);
      expect(descriptorDifferences(expected, actual).join('; '), '', 'the values, against the container\'s');
    });
  }, {timeout: 120000 + CONTAINER_TIMEOUT});

  test('descriptors: a row Crux cannot read comes from the container', async () => {
    // a gap in Crux's reader (a row of nulls for a molecule given) is recomputed on the container, counted and logged
    await ensureContainerRunning('name = "chem-chem"', CONTAINER_TIMEOUT);
    await inEngine(ChemicalEngine.Crux, async () => {
      const molecules = DG.Column.fromStrings('molecules', [ASPIRIN, 'CCO', 'c1ccccc1O', '', 'C1CC(']);
      const expected = await container.calculateDescriptors(molecules, DESCRIPTORS_2D);
      setCruxDescriptorGapForTests(['CCO']);
      try {
        const fallbacks = cruxFallbackCount('calcDescriptors');
        const actual = await calculateDescriptors(molecules, DESCRIPTORS_2D);
        expect(cruxFallbackCount('calcDescriptors'), fallbacks + 1, 'one fallback: the row Crux could not read');
        const rowsWithValues = (columns: DG.Column[]) => columns.length ? Array.from({length: 3}, (_, i) => i) : [];
        const pick = (columns: DG.Column[]) => columns.map((c) => DG.Column.fromList(c.type as DG.ColumnType, c.name,
          rowsWithValues(columns).map((i) => c.get(i))));
        expect(descriptorDifferences(pick(expected), pick(actual)).join('; '), '', 'the three molecules\' values');
        expect(actual.every((c) => c.isNone(4)), true, 'a molecule neither reads stays empty');
      } finally {
        setCruxDescriptorGapForTests(null);
      }
    });
  }, {timeout: 60000 + CONTAINER_TIMEOUT});

  test('descriptors: the typed columns are the list-built ones', async () => {
    // batch 22: the columns made from typed arrays, each distinct molecule calculated once, hold what
    // DG.Column.fromList(type, name, values with nulls) held: types, storage, values and missing values
    await inEngine(ChemicalEngine.Crux, async () => {
      const molecules = DG.Column.fromStrings('molecules', [ASPIRIN, '', 'C1CC(', 'CCO', ASPIRIN, 'c1ccccc1O', 'CCO']);
      const names = ['MolWt', 'NumHDonors', 'TPSA', 'RingCount', 'MolLogP', 'FractionCSP3'];
      const actual = await calculateDescriptors(molecules, names);
      const packed = await (await getCruxService()).calcDescriptors(molecules.toList(), names);
      const differences: string[] = [];
      for (const [name, values] of Object.entries(packed)) {
        const expected = DG.Column.fromList(values.type, name,
          Array.from(values.value, (v, i) => isNoneValue(values, i) ? null : v));
        const column = actual.find((c) => c.name === name);
        if (!column || column.type !== expected.type ||
          column.getRawData().constructor !== expected.getRawData().constructor) {
          differences.push(`${name}: ${column?.type} ${column?.getRawData().constructor.name}, list-built ` +
            `${expected.type} ${expected.getRawData().constructor.name}`);
          continue;
        }
        for (let i = 0; i < expected.length; i++) {
          if (column.isNone(i) !== expected.isNone(i) || (!expected.isNone(i) && column.get(i) !== expected.get(i)))
            differences.push(`${name} row ${i}: ${column.get(i)}, list-built ${expected.get(i)}`);
        }
      }
      expect(differences.join('; '), '');
    });
  });

  test('descriptors: 1,000 molecules use all the workers, with the same values', async () => {
    // batch 22: fewer molecules than 1,000 per worker are split evenly over the workers; chunking changes no value
    await inEngine(ChemicalEngine.Crux, async () => {
      const smiles = (await readDataframe('smiles.csv')).col('canonical_smiles')!.toList() as string[];
      const molecules = DG.Column.fromStrings('molecules', smiles);
      const names = ['MolWt', 'NumHDonors', 'TPSA', 'MolLogP', 'qed', 'FractionCSP3'];
      setDescriptorTimingForTests(true);
      let together: DG.Column[];
      let timing: NonNullable<ReturnType<typeof lastDescriptorTiming>>;
      try {
        together = await calculateDescriptors(molecules, names);
        timing = lastDescriptorTiming()!; // (read before the switch goes off, which forgets it)
      } finally {
        setDescriptorTimingForTests(false);
      }
      const workers = new Set(timing.chunks.map((c) => c.worker)).size;
      const distinct = molecules.categories.length;
      // the chunks are listed as they finish, so the last, smaller one may come first
      const largest = Math.max(...timing.chunks.map((c) => c.size));
      const what = `${timing.chunks.length} chunks of up to ${largest} on ${workers} of ${timing.workers} workers`;
      expect(largest, Math.ceil(distinct / timing.workers), what);
      expect(workers, timing.chunks.length, `one chunk per worker: ${what}`);
      expect(workers > timing.workers / 2, true, what);
      // the same molecules in tens: each call its own chunking
      const differences: string[] = [];
      for (let start = 0; start < smiles.length; start += 10) {
        const part = await calculateDescriptors(DG.Column.fromStrings('molecules', smiles.slice(start, start + 10)),
          names);
        for (const column of part) {
          const whole = together.find((c) => c.name === column.name)!;
          for (let i = 0; i < column.length; i++) {
            if (column.isNone(i) !== whole.isNone(start + i) || (!column.isNone(i) && column.get(i) !== whole.get(start + i)))
              differences.push(`${column.name} row ${start + i}`);
          }
        }
      }
      expect(differences.slice(0, 10).join(', '), '', 'values that depend on the chunking');
    });
  }, {timeout: 300000});

  test('descriptors: 50,000 rows of 1,000 molecules within 2 s', async () => {
    // batch 22's budget: each distinct molecule once and typed arrays (before it: 4.8 s for these 20 descriptors). A
    // timing: with --benchmark only
    if (!DG.Test.isInBenchmark)
      return;
    await inEngine(ChemicalEngine.Crux, async () => {
      const smiles = (await readDataframe('smiles.csv')).col('canonical_smiles')!.toList() as string[];
      const molecules = DG.Column.fromStrings('molecules', Array.from({length: 50000}, (_, i) => smiles[i % smiles.length]));
      // the tree's first 20 2D descriptors, as the timing test's middle set
      const names = (Object.values(await getDescriptorsTree()) as any[]).flatMap((g) => g.descriptors)
        .filter((d: any) => !d.tags?.includes('3D')).map((d: any) => d.name as string).slice(0, 20);
      await calculateDescriptors(DG.Column.fromStrings('warm', [ASPIRIN]), ['MolWt']); // the pool up
      const started = performance.now();
      const columns = await calculateDescriptors(molecules, names);
      const ms = performance.now() - started;
      console.log(`descriptor timing | budget: 50,000 rows of 1,000 molecules, ${names.length} names: ${ms.toFixed(0)} ms`);
      expect(columns.length, names.length);
      expect(ms < 2000, true, `${ms.toFixed(0)} ms`);
    });
  }, {timeout: 120000, benchmark: true});

  test('descriptors: Crux timing on 1,000 and 5,000 molecules', async () => {
    // the small inputs batch 22 splits over all the workers: smiles.csv, and 5,000 unique molecules (each with one of
    // five small extra fragments), with 3 and 20 descriptors; timings logged. A measurement, no check: with
    // --benchmark only
    if (!DG.Test.isInBenchmark)
      return;
    await inEngine(ChemicalEngine.Crux, async () => {
      const smiles = (await readDataframe('smiles.csv')).col('canonical_smiles')!.toList() as string[];
      const names20 = (Object.values(await getDescriptorsTree()) as any[]).flatMap((g) => g.descriptors)
        .filter((d: any) => !d.tags?.includes('3D')).map((d: any) => d.name as string).slice(0, 20);
      const extras = ['C', 'CC', 'CO', 'CN', 'CCl'];
      const tables: [string, DG.Column][] = [['1,000', DG.Column.fromStrings('molecules', smiles)],
        ['5,000', DG.Column.fromStrings('molecules', Array.from({length: 5000},
          (_, i) => `${smiles[i % smiles.length]}.${extras[Math.floor(i / smiles.length)]}`))]];
      await calculateDescriptors(DG.Column.fromStrings('warm', [ASPIRIN]), ['MolWt']); // the pool up
      setDescriptorTimingForTests(true);
      try {
        for (const [label, molecules] of tables) {
          for (const names of [['MolWt', 'TPSA', 'MolLogP'], names20]) {
            const started = performance.now();
            await calculateDescriptors(molecules, names);
            const ms = performance.now() - started;
            const timing = lastDescriptorTiming()!;
            console.log(`descriptor timing | ${label} molecules (${molecules.categories.length} distinct), ` +
              `${names.length} names: ${ms.toFixed(0)} ms; ${timing.chunks.length} chunks of ${timing.chunks[0].size} on ` +
              `${new Set(timing.chunks.map((c) => c.worker)).size} of ${timing.workers} workers`);
          }
        }
      } finally {
        setDescriptorTimingForTests(false);
      }
    });
  }, {timeout: 300000, benchmark: true});

  test('descriptors: Crux timing on a large table', async () => {
    // smiles.csv's molecules over and over, 500,000 rows with 3, 20 and all the 2D descriptors, repeated and then
    // unique (each with one of 500 small extra fragments: no two rows alike); each call's phases, pool and chunks
    // logged. A measurement: with --benchmark only (the budget test above checks 50,000 rows)
    if (!DG.Test.isInBenchmark)
      return;
    await inEngine(ChemicalEngine.Crux, async () => {
      const rows = 500000;
      const smiles = (await readDataframe('smiles.csv')).col('canonical_smiles')!.toList() as string[];
      const suffixes = ['', 'O', 'N', 'Cl', 'F', 'Br', 'S', 'C=O', 'C#N', 'OC', 'NC', 'CO', 'CN', 'C(=O)O', 'C(=O)N',
        'OCC', 'NCC', 'SC', 'C(C)C', 'C(C)(C)C', 'CF', 'CCl', 'C(F)(F)F', 'OO', 'NN'];
      const fragment = (k: number) => 'C'.repeat(1 + k % 20) + suffixes[Math.floor(k / 20) % suffixes.length];
      const tables: [string, DG.Column][] = [
        ['repeated', DG.Column.fromStrings('molecules', Array.from({length: rows}, (_, i) => smiles[i % smiles.length]))],
        ['unique', DG.Column.fromStrings('molecules', Array.from({length: rows},
          (_, i) => `${smiles[i % smiles.length]}.${fragment(Math.floor(i / smiles.length))}`))]];
      const tree = await getDescriptorsTree();
      const names2D = (Object.values(tree) as any[]).flatMap((g) => g.descriptors)
        .filter((d: any) => !d.tags?.includes('3D')).map((d: any) => d.name as string);
      // all the 2D descriptors on the first 100,000 rows only (some 200 columns of 500,000 would not fit a tab)
      const sets: [string[], number][] =
        [[['MolWt', 'TPSA', 'MolLogP'], rows], [names2D.slice(0, 20), rows], [names2D, 100000]];
      setDescriptorTimingForTests(true);
      const runs = tables.flatMap(([kind, column]) => sets.map(([names, n]) => ({kind, names, molecules: n === rows ?
        column : DG.Column.fromStrings('molecules', (column.toList() as string[]).slice(0, n))})));
      try {
        for (const {kind, molecules, names} of runs) {
          const started = performance.now();
          const columns = await calculateDescriptors(molecules, names);
          const calculated = performance.now();
          const table = DG.DataFrame.create(molecules.length);
          for (const col of columns)
            table.columns.add(col);
          const added = performance.now();
          const timing = lastDescriptorTiming()!;
          const chunkMs = timing.chunks.map((c) => c.ms).sort((a, b) => a - b);
          const sum = (f: (c: typeof timing.chunks[0]) => number) => timing.chunks.reduce((s, c) => s + f(c), 0);
          console.log(`descriptor timing | ${molecules.length} ${kind} rows (${molecules.categories.length} distinct), ` +
            `${names.length} names: ${(calculated - started).toFixed(0)} ms ` +
            `(adding the columns ${(added - calculated).toFixed(0)} ms); ${timing.workers} workers, ` +
            `${timing.chunks.length} chunks of ${timing.chunks[0]?.size}; phases ` +
            Object.entries(timing.phases).map(([k, v]) => `${k} ${v.toFixed(0)}`).join(', ') + '; per chunk median ' +
            `${chunkMs[Math.floor(chunkMs.length / 2)]?.toFixed(0)} ms, Crux ${(sum((c) => c.wasmMs) /
              timing.chunks.length).toFixed(0)} ms, parse ${(sum((c) => c.parseMs) / timing.chunks.length).toFixed(0)} ms; ` +
            `throughput ${(timing.rows / ((timing.phases['chunks'] ?? 1) / 1000)).toFixed(0)} molecules/s; storage ` +
            [...new Set(columns.map((c) => `${c.type} ${c.getRawData().constructor.name}`))].join(', '));
          expect(columns.length, names.length);
        }
      } finally {
        setDescriptorTimingForTests(false);
      }
    });
  }, {timeout: 600000, benchmark: true, benchmarkTimeout: 3600000});

  test('scaffold tree: Crux against the Python script', async () => {
    // smiles.csv's first 250 rows (the 100,000-molecule comparison, misc/crux/compare100k.mjs, is the large-scale check)
    const datasets: [string, string, number][] = [
      ['tests/sar-small_test.csv', 'smiles', 1000], ['smiles.csv', 'canonical_smiles', 250],
      ['tests/spgi-100.csv', 'Structure', 1000]];
    const differences: string[] = [];
    const layouts: string[] = [];
    for (const [file, column, rows] of datasets) {
      const all = await readDataframe(file);
      const df = DG.DataFrame.fromColumns([all.col(column)!]);
      if (df.rowCount > rows)
        df.rows.removeAt(rows, df.rowCount - rows);
      await grok.data.detectSemanticTypes(df);
      const python = JSON.parse(await inEngine(ChemicalEngine.RDKit, () => PackageFunctions.getScaffoldTree(df, 10, false)));
      // Crux built it: no fallback to the script
      const crux = JSON.parse(await inEngine(ChemicalEngine.Crux, () =>
        expectNoFallback(['scaffoldTree'], () => PackageFunctions.getScaffoldTree(df, 10, false))));
      // the structure (nodes, SMILES, children, their order) must be the script's; the nodes' MOL blocks are compared
      // apart: wasm's libm moves a few coordinates in their last digits
      const [pythonShape, cruxShape] = [python, crux].map((t) => JSON.stringify(t, (k, v) => k === 'scaffold' ? undefined : v));
      if (cruxShape !== pythonShape)
        differences.push(`${file}: ${firstDifference(JSON.parse(pythonShape), JSON.parse(cruxShape), '')}`);
      const [pythonBlocks, cruxBlocks] = [python, crux].map((t) => molblocksOf(t));
      layouts.push(`${file}: ${pythonBlocks.filter((b, i) => b !== cruxBlocks[i]).length} of ${pythonBlocks.length} MOL blocks differ`);
    }
    console.log(`chemical engine | scaffold tree MOL blocks, Crux against the script: ${layouts.join('; ')}`);
    expect(differences.join('; '), '', 'the trees');
  }, {timeout: 600000});

  test('the substructure search follows the Chemical Engine', async () => {
    // no setting of its own since 2026-10-22 (the product owner): RDKit mode searches with RDKit, Crux mode with Crux,
    // whatever a server saved for the retired SubstructureSearchEngine (this stand's maps hold Crux)
    const col = (await readDataframe('smiles.csv')).col('canonical_smiles')!;
    col.semType = DG.SEMTYPE.MOLECULE;
    setSubstructureSearchEngine(null);
    const search = async () => (await PackageFunctions.searchSubstructure(col, 'c1ccncc1', '')).get(0) as DG.BitSet;
    const rdkit = await inEngine(ChemicalEngine.RDKit, async () => {
      expect(getSubstructureSearchEngine(), SubstructureSearchEngine.RDKit, 'RDKit mode: the search engine');
      expect(isCruxSearch(SubstructureSearchType.CONTAINS), false, 'RDKit mode: Contains on RDKit');
      return search();
    });
    const crux = await inEngine(ChemicalEngine.Crux, async () => {
      expect(getSubstructureSearchEngine(), SubstructureSearchEngine.Crux, 'Crux mode: the search engine');
      expect(isCruxSearch(SubstructureSearchType.NOT_CONTAINS), true, 'Crux mode: Not contains on Crux');
      return search();
    });
    expect(crux.toString(), rdkit.toString(), 'the same hits in both modes');
  });

  test('the property calculator is Chemical Properties, its method the engine\'s', async () => {
    // the product owner, 2026-10-22: "(OCL)" leaves the calculator's shown name (its function name stays, for saved
    // HitTriage campaigns and flows); in Crux mode the dialog names RDKit's definitions as its method
    const funcs = DG.Func.find({package: 'Chem', meta: {function_family: 'biochem-calculator'}});
    const func = funcs.find((f) => f.name === 'addChemPropertiesColumns')!;
    expect(!!func, true, `the OpenChemLib calculator among ${funcs.map((f) => f.name).join(', ')}`);
    expect(func?.friendlyName, 'Chemical Properties', 'the calculator\'s name');
    await inEngine(ChemicalEngine.RDKit, async () =>
      expect(calculatorMethodInfo(func).authors, 'Open Chem Lib Team', 'RDKit mode: OpenChemLib'));
    await inEngine(ChemicalEngine.Crux, async () => {
      expect(calculatorMethodInfo(func).authors, 'RDKit\'s definitions, calculated by Crux', 'Crux mode: the method');
      expect(calculatorMethodInfo(func).github, 'https://github.com/rdkit/rdkit', 'Crux mode: the source');
    });
  });

  // batch 19b: fingerprints and structural alerts in Crux mode

  /** Chem's molecules for the worker jobs: SMILES (smiles.csv), MOL blocks (spgi-100, mol1K), '' and an unreadable one. */
  const jobMolecules = async () => [
    ...(await readDataframe('smiles.csv')).col('canonical_smiles')!.toList(),
    ...(await readDataframe('tests/spgi-100.csv')).col('Structure')!.toList(),
    ...(await readDataframe('mol1K.csv')).col('molecule')!.toList().slice(0, 300),
    '', 'C1CC(',
  ] as string[];

  test('fingerprints: Crux against the RDKit pool', async () => {
    const molecules = await jobMolecules();
    const rdkitPool = await chemCommonRdKit.getRdKitService();
    const cruxPool = await getCruxService();
    const report: string[] = [];
    const tooMany: string[] = [];
    for (const type of Object.values(Fingerprint)) {
      const rdkit = await rdkitPool.getFingerprints(type, molecules, true);
      const crux = await cruxPool.getFingerprints(type, molecules, true);
      const differ: number[] = [];
      for (let i = 0; i < molecules.length; i++) {
        const [a, b] = [rdkit.fps[i], crux.fps[i]];
        if ((a === null) !== (b === null) || (a && b && (a.length !== b.length || a.some((v, k) => v !== b[k]))))
          differ.push(i);
      }
      // every differing molecule listed: RDKit 2024.09 (Chem's) against 2026.03 (Crux's)
      report.push(`${type} ${molecules.length - differ.length}/${molecules.length} identical` +
        (differ.length ? ` (differ: ${differ.map((i) => `#${i} ${rdkit.smiles?.[i] ?? crux.smiles?.[i]}`).join(', ')})` : ''));
      if (differ.length > molecules.length / 50)
        tooMany.push(`${type}: ${differ.length}`);
    }
    console.log(`chemical engine | fingerprints, Crux against the RDKit pool on ${molecules.length} molecules: ` +
      report.join('; '));
    expect(tooMany.join(', '), '', 'fingerprint types with more than 2 % of the molecules differing');
  }, {timeout: 300000});

  test('fingerprints and alerts: RDKit does the molecules Crux cannot read', async () => {
    const PHENOL = 'Oc1ccccc1';
    const molecules = [ASPIRIN, PHENOL, 'CCN', '', 'C1CC('];
    await inEngine(ChemicalEngine.Crux, () => withCruxUnreadable([PHENOL], async () => {
      const fallbacks = cruxFallbackCount('fingerprints');
      const crux = await fingerprintsByEngine(Fingerprint.Morgan, molecules, false);
      const rdkit = await (await chemCommonRdKit.getRdKitService()).getFingerprints(Fingerprint.Morgan, molecules, false);
      expect(crux.engine, ChemicalEngine.Crux);
      expect(crux.fps[1]?.toString(), rdkit.fps[1]?.toString(), 'the molecule Crux could not read: RDKit\'s fingerprint');
      expect(rdkit.fps[1] !== null, true);
      expect(crux.fps[4], null, 'a molecule neither engine reads: none');
      expect(cruxFallbackCount('fingerprints'), fallbacks + 1, 'one fallback: the molecule RDKit read');
      const alerts = await structuralAlertsByEngine({'PAINS': ['[OX2H]c1ccccc1']}, molecules);
      expect(alerts[0][1].join(','), 'false,true,false,false,false', 'RDKit matched the molecule Crux could not read');
    }));
  });

  test('fingerprints and alerts: RDKit does a CXSMILES query and an alert SMARTS Crux cannot read', async () => {
    const QUERY = '*c1ccccc1 |$Q_e;;;;;;$|'; // a query atom's label, which Crux's reader skips (its pattern bits differ)
    const BENZENE = 'c1ccccc1';
    const molecules = [ASPIRIN, QUERY, 'CCN'];
    await inEngine(ChemicalEngine.Crux, async () => {
      const crux = await fingerprintsByEngine(Fingerprint.Pattern, molecules, false);
      const rdkit = await (await chemCommonRdKit.getRdKitService()).getFingerprints(Fingerprint.Pattern, molecules, false);
      expect(crux.fps.map((fp) => fp?.toString()).join('|'), rdkit.fps.map((fp) => fp?.toString()).join('|'),
        'the CXSMILES query: RDKit\'s pattern fingerprint');
      // the SMARTS: one Crux reads nothing from
      const alerts = await withCruxUnreadable([BENZENE], () =>
        structuralAlertsByEngine({'PAINS': ['[Se]', BENZENE], 'LINT': ['[Se]']}, molecules));
      expect(alerts.map(([rule, flags]) => `${rule} ${flags.join(',')}`).join('; '),
        'PAINS true,true,false; LINT false,false,false', 'the SMARTS Crux cannot read: RDKit\'s matches');
    });
  });

  test('fingerprints: one engine\'s cache at a time, edits reach both engines', async () => {
    const df = await readDataframe('smiles.csv');
    const col = df.col('canonical_smiles')!;
    col.semType = DG.SEMTYPE.MOLECULE;
    const cached = () => ['~canonical_smiles.Morgan', '~canonical_smiles.Morgan.Crux'].filter((n) => !!df.col(n)).join(', ');
    const asText = (fps: any[]) => fps.map((fp) => fp?.toString() ?? '');
    const morgan = async () => asText(await chemGetFingerprints(col, Fingerprint.Morgan, false));
    const fresh = async () => { // the fingerprints of the column's current molecules, without its caches
      const copy = DG.Column.fromStrings('smiles', col.toList());
      DG.DataFrame.fromColumns([copy]).name = 'fresh';
      copy.semType = DG.SEMTYPE.MOLECULE;
      return asText(await chemGetFingerprints(copy, Fingerprint.Morgan, false));
    };
    const rdkit = await inEngine(ChemicalEngine.RDKit, morgan);
    expect(cached(), '~canonical_smiles.Morgan');
    // Crux computed them, then an edit read first in Crux mode (it consumes the column's change set)
    const crux = await inEngine(ChemicalEngine.Crux, () => expectNoFallback(['fingerprints'], async () => {
      const fps = await morgan();
      col.set(5, 'c1ccccc1O');
      const edited = await morgan();
      expect(edited.every((fp, i) => fp === (i === 5 ? edited[5] : fps[i])), true);
      return fps;
    }));
    const same = rdkit.filter((fp, i) => fp === crux[i]).length;
    console.log(`chemical engine | Morgan of smiles.csv through chemGetFingerprints: ${same}/${col.length} identical`);
    expect(same >= col.length * 0.99, true);
    // then in RDKit mode, and back
    await inEngine(ChemicalEngine.RDKit, async () => {
      const rdkitEdited = await morgan();
      const rdkitFresh = await fresh();
      expect(rdkitEdited[5], rdkitFresh[5], 'RDKit mode sees the edit Crux mode read');
      expect(rdkitEdited.every((fp, i) => fp === rdkitFresh[i]), true, 'RDKit mode: the fingerprints of the molecules');
      expect(cached(), '~canonical_smiles.Morgan');
      col.set(6, 'CCO');
      expect((await morgan())[6], (await fresh())[6], 'RDKit mode reads its edit');
    });
    await inEngine(ChemicalEngine.Crux, async () => {
      const cruxAgain = await morgan();
      const cruxFresh = await fresh();
      expect(cruxAgain[6], cruxFresh[6], 'Crux mode sees the edit RDKit mode read');
      expect(cruxAgain.every((fp, i) => fp === cruxFresh[i]), true, 'Crux mode: the fingerprints of the molecules');
      expect(cached(), '~canonical_smiles.Morgan.Crux', 'one engine\'s cache at a time');
    });
  }, {timeout: 120000});

  test('Crux mode leaves the RDKit search its own fingerprints after edits', async () => {
    // the RDKit substructure search refreshes its cached pattern fingerprints for edited rows: RDKit's, in Crux mode too
    const rdkitSearch = (col: DG.Column, query: string) => {
      setSubstructureSearchEngine(SubstructureSearchEngine.RDKit);
      return chemSubstructureSearchLibrary(col, query, '', FILTER_TYPES.substructure, false, true);
    };
    await inEngine(ChemicalEngine.Crux, async () => {
      const col = (await readDataframe('smiles.csv')).col('canonical_smiles')!;
      await rdkitSearch(col, 'C1CC1');
      col.set(0, 'C1CC1CCN');
      col.set(500, 'c1ccccc1');
      col.set(998, '');
      await rdkitSearch(col, 'C1CC1'); // refreshes the edited rows' cached fingerprints
      const edited = await rdkitSearch(col, 'c1ccccc1'); // reads them
      const fresh = DG.Column.fromStrings('smiles', col.toList());
      DG.DataFrame.fromColumns([fresh]).name = 'fresh';
      fresh.semType = DG.SEMTYPE.MOLECULE;
      const expected = await rdkitSearch(fresh, 'c1ccccc1');
      expect(edited.equals(expected), true, `after edits ${edited.trueCount()} hits, without caches ${expected.trueCount()}`);
      setSubstructureSearchEngine(null);
    });
  }, {timeout: 120000});

  test('similarity and diversity in both engines', async () => {
    const col = (await readDataframe('smiles.csv')).col('canonical_smiles')!;
    col.semType = DG.SEMTYPE.MOLECULE;
    const query = 'O=C(O)c1ccccc1O';
    const {rdkit, crux} = await inBothEngines(async () => ({
      similarities: (await chemGetSimilarities(col, query))!.toList() as number[],
      diverse: (await chemGetDiversities(col, 10))!.toList() as string[],
    }));
    const differ = rdkit.similarities.filter((v, i) => Math.abs(v - crux.similarities[i]) > 1e-12).length;
    // the diverse subset starts from a random molecule (getDiverseSubset): its members are not compared, its size is
    console.log(`chemical engine | similarity to ${query}: ${col.length - differ}/${col.length} identical`);
    expect(differ <= col.length / 100, true, `${differ} similarities differ`);
    expect(crux.diverse.length, rdkit.diverse.length, 'the diverse subsets');
  }, {timeout: 120000});

  test('structural alerts: Crux against the RDKit pool', async () => {
    const alertsDf = DG.DataFrame.fromCsv(await _package.files.readAsText('alert-collection.csv'));
    const molecules = DG.Column.fromStrings('molecules', await jobMolecules());
    const ruleSet: RuleSet = {'BMS': true, 'Dundee': true, 'Glaxo': true, 'Inpharmatica': true, 'LINT': true,
      'MLSMR': true, 'PAINS': true, 'SureChEMBL': true};
    const rdkit = await runStructuralAlertsDetection(molecules, ruleSet, alertsDf,
      await chemCommonRdKit.getRdKitService());
    // Crux matched them
    const crux = await inEngine(ChemicalEngine.Crux, () => expectNoFallback(['structural alerts'], () =>
      runStructuralAlertsDetection(molecules, ruleSet, alertsDf, alertsByEngine)));
    const report: string[] = [];
    const differences: string[] = [];
    for (const column of rdkit.columns.toList()) {
      const other = crux.col(column.name)!;
      const rows = Array.from({length: column.length}, (_, i) => i).filter((i) => column.get(i) !== other.get(i));
      report.push(`${column.name} ${column.stats.sum} hits`);
      if (rows.length)
        differences.push(`${column.name}: rows ${rows.slice(0, 10).join(', ')}${rows.length > 10 ? '...' : ''}`);
    }
    console.log(`chemical engine | structural alerts on ${molecules.length} molecules (RDKit): ${report.join('; ')}`);
    expect(differences.join('; '), '', 'the alerts per rule');
  }, {timeout: 300000});

  test('pharmacophore features: Crux against the RDKit pool', async () => {
    const featuresDf = await readDataframe('pharmacophore-features.csv');
    const molecules = DG.Column.fromStrings('molecules', await jobMolecules());
    const families = Object.fromEntries(PHARM_FAMILY_NAMES.map((f) => [f, true])) as PharmFamilySet;
    // Crux matched them
    const {rdkit, crux} = await inBothEngines(() =>
      runPharmacophoreDetection(molecules, families, featuresDf, alertsByEngine), {noFallback: ['structural alerts']});
    const report: string[] = [];
    const differences: string[] = [];
    for (const column of rdkit.columns.toList()) {
      const other = crux.col(column.name)!;
      const rows = Array.from({length: column.length}, (_, i) => i).filter((i) => column.get(i) !== other.get(i));
      report.push(`${column.name} ${column.stats.sum}`);
      if (rows.length)
        differences.push(`${column.name}: rows ${rows.slice(0, 10).join(', ')}${rows.length > 10 ? '...' : ''}`);
    }
    console.log(`chemical engine | pharmacophore features on ${molecules.length} molecules (RDKit): ${report.join('; ')}`);
    expect(rdkit.columns.length > 0, true, 'the families');
    expect(differences.join('; '), '', 'the features per family');
  }, {timeout: 300000});

  // batch 23 (2026-10-23): conversions, validation, the one-molecule panels and the MMP viewer's joins on Crux

  /** RDKit's canonical SMILES of a molecule string (a query's: of its SMARTS), null where RDKit reads none. */
  const rdkitCanonical = (s: string, query = false): string | null => {
    const module = chemCommonRdKit.getRdKitModule();
    let mol = null;
    try {
      mol = query ? module.get_qmol(s) : module.get_mol(s);
      return mol ? (query ? mol.get_smarts() : mol.get_smiles()) : null;
    } catch {
      return null;
    } finally {
      mol?.delete();
    }
  };

  test('conversions: Crux writes the same molecules, RDKit those Crux cannot read', async () => {
    // Crux writes RDKit.js 2026.03.6's strings (the product owner, 2026-10-23); RDKit mode's are 2024.09's: the same
    // molecules, read back by RDKit (a SMARTS: the molecule matches it)
    const all = (await jobMolecules()).filter((m) => m);
    const molecules = [...all.slice(0, 300), ...all.slice(999, 1199)]; // SMILES, then spgi-100's and mol1K's MOL blocks
    const targets = [DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock, DG.chem.Notation.V3KMolBlock,
      DG.chem.Notation.Smarts, 'cxsmiles' as DG.chem.Notation];
    const convertAll = async (): Promise<[string[][], string[]]> => [
      targets.map((t) => molecules.map((m) => PackageFunctions.convertMolNotation(m, DG.chem.Notation.Unknown, t))),
      await convertNotationForColumn(DG.Column.fromStrings('m', molecules), DG.chem.Notation.Smiles)];
    // Crux converted every molecule
    const {rdkit: [rdkit, rdkitColumn], crux: [crux, cruxColumn]} = await inBothEngines(convertAll,
      {noFallback: ['convertMolNotation']});
    const differences: string[] = [];
    let identical = 0;
    targets.forEach((t, k) => molecules.forEach((m, i) => {
      const [r, c] = [rdkit[k][i], crux[k][i]];
      identical += r === c ? 1 : 0;
      const same = r === c || (t === DG.chem.Notation.Smarts ? rdkitCanonical(c, true) !== null :
        rdkitCanonical(r) === rdkitCanonical(c));
      if (!same && differences.length < 5)
        differences.push(`${t} of ${m.slice(0, 60)}: ${c.slice(0, 80)} vs RDKit's ${r.slice(0, 80)}`);
    }));
    molecules.forEach((m, i) => {
      if (rdkitCanonical(rdkitColumn[i]) !== rdkitCanonical(cruxColumn[i]) && differences.length < 10)
        differences.push(`the column's SMILES of ${m.slice(0, 60)}: ${cruxColumn[i]} vs RDKit's ${rdkitColumn[i]}`);
    });
    console.log(`chemical engine | conversions: ${identical} of ${targets.length * molecules.length} strings ` +
      `identical to RDKit 2024.09's`);
    expect(differences.join('; '), '', 'the same molecules');
    // a molecule Crux cannot read: the RDKit pool converts it in the column
    await inEngine(ChemicalEngine.Crux, () => withCruxUnreadable([molecules[0]], async () => {
      const fallbacks = cruxFallbackCount('convertMolNotation');
      const column = await convertNotationForColumn(DG.Column.fromStrings('m', molecules.slice(0, 3)),
        DG.chem.Notation.Smiles);
      expect(column[0], rdkitColumn[0], 'RDKit\'s string for the molecule Crux cannot read');
      expect(cruxFallbackCount('convertMolNotation') > fallbacks, true, 'the fallback is counted');
    }));
  }, {timeout: 300000});

  test('InChI: Crux mode writes RDKit\'s InChIs and keys, the panel\'s and the column\'s', async () => {
    // the split (crux-core decision chem-wasm-inchi): Crux fills InChI's input, crux-sketch's InChI module runs the
    // library; RDKit.js 2026.03.6's strings, which Chem's RDKit 2024.09 writes too on these molecules (200 SMILES, 50
    // MOL blocks of SPGI and 50 of mol1K; compare100k.mjs compares 100,000)
    const all = (await jobMolecules()).filter((m) => m);
    const molecules = [...all.slice(0, 200), ...all.slice(1000, 1050), ...all.slice(1100, 1150), '', 'C1CC(', '*CC',
      'C[Na]', 'c1cccc1', '[2H]C([2H])O', 'C/C=C/C', 'C[C@H](N)O |&1:1|'];
    // Crux wrote every InChI (RDKit none)
    const {rdkit: [rdkit, rdkitKeys], crux: [crux, cruxKeys]} = await inBothEngines(async () =>
      [await inchisByEngine(molecules, false), await inchisByEngine(molecules, true)], {noFallback: ['InChI']});
    const differences = molecules.flatMap((m, i) => crux[i] === rdkit[i] && cruxKeys[i] === rdkitKeys[i] ? [] :
      [`${m.slice(0, 50)}: ${crux[i]} ${cruxKeys[i]} vs RDKit's ${rdkit[i]} ${rdkitKeys[i]}`]);
    expect(differences.slice(0, 5).join('; '), '', `the same InChIs and keys (${differences.length} differ)`);
    const panel = (engine: ChemicalEngine) => inEngine(engine, async () => {
      const identifiers = await getIdentifiersSingle(molecules[0]) as {[k: string]: string};
      return `${identifiers['Inchi']} ${identifiers['Inchi key']}`;
    });
    expect(await panel(ChemicalEngine.Crux), await panel(ChemicalEngine.RDKit), 'the Identifiers panel');
    await inEngine(ChemicalEngine.Crux, async () => {
      expect((await cruxInchi(molecules[0]))?.inchi, rdkit[0], 'the panel\'s InChI is Crux\'s');
      expect(await withCruxUnreadable([molecules[1]], async () => (await inchisByEngine(molecules.slice(0, 3), false))[1]),
        rdkit[1], 'RDKit\'s for an unreadable molecule');
    });
  }, {timeout: 300000});

  test('isSmiles and validateMolecule: Crux mode answers as RDKit mode', async () => {
    const strings = ['CCO', 'c1ccccc1O', 'C[C@H](N)C(=O)O', '[#6]~[#7]', 'c1cccc1', 'C1CC(', 'not a molecule', '',
      'C(C)(C)(C)(C)C', ...(await jobMolecules()).slice(-150, -100)];
    const answers = async () => strings.map((s) => // (RDKit's log lines start with their time)
      `${PackageFunctions.isSmiles(s)}|${PackageFunctions.validateMolecule(s)?.replace(/\[\d\d:\d\d:\d\d\] /g, '')}`);
    const {rdkit, crux} = await inBothEngines(answers);
    const differ = strings.flatMap((s, i) => rdkit[i] === crux[i] ? [] : [`${s.slice(0, 40)}: ${crux[i]} vs ${rdkit[i]}`]);
    expect(differ.join('; '), '');
  });

  test('structural alerts and pharmacophores panels: Crux matches as RDKit', async () => {
    const all = (await jobMolecules()).filter((m) => m && m !== 'C1CC(');
    const molecules = [...all.slice(0, 40), ...all.slice(999, 1039)]; // SMILES and MOL blocks
    const panels = async () => {
      const out: string[] = [];
      for (const m of molecules) {
        const features = (await getPharmacophoreFeatures(m)).map((f) => `${f.index}:${f.atoms.join(',')}`);
        out.push(`${(await getStructuralAlerts(m)).join(',')}|${features.join(';')}`);
      }
      return out;
    };
    // Crux matched every molecule
    const {rdkit, crux} = await inBothEngines(panels, {noFallback: ['structuralAlerts', 'pharmacophoreFeatures']});
    const differ = molecules.flatMap((m, i) => rdkit[i] === crux[i] ? [] : [m.slice(0, 60)]);
    expect(differ.join('; '), '', 'the same alerts and features');
    expect(rdkit.some((r) => r.split('|')[0] !== ''), true, 'some molecule has an alert');
  }, {timeout: 300000});

  /** Molecule strings meant to break a Crux path (2026-10-23 review): explicit H, a query block, a dative bond, a core,
   * enhanced stereo, CXSMILES query features, salts, isotopes, radicals, SMARTS, unreadable and over-long strings. */
  const breakers = (): {[what: string]: string} => {
    const rdkit = chemCommonRdKit.getRdKitModule();
    const block = (s: string, query = false, hs = false) => {
      const mol = query ? rdkit.get_qmol(s) : rdkit.get_mol(s)!;
      try {
        if (hs)
          mol.add_hs_in_place();
        return mol.get_molblock();
      } finally {
        mol.delete();
      }
    };
    return {
      'explicit H block': block('C[C@H](N)C(=O)O', false, true), 'query block': block('[#7,#8]-c1ccccc1', true),
      'dative block': block('[NH3]->[Pt+2]<-[NH3]'), 'core': 'c1ccc([*:1])cc1', 'two-point core': '[*:1]C(=O)N[*:2]',
      'AND group': 'C[C@H](N)C(=O)O |&1:1|', 'OR group': 'C[C@H](O)[C@@H](C)O |o1:1,3|',
      'query label': '*c1ccccc1 |$Q_e;;;;;;$|', 'ring bond count': 'C1CCCCC1 |rb:0:2|', 'N-oxide': '[O-][n+]1ccccc1',
      'dative SMILES': '[NH3]->[Pt+2]<-[NH3]', 'salt': 'CC(=O)[O-].[Na+]', 'mixture': 'c1ccccc1.O',
      'isotopes': '[2H]C([2H])([2H])Oc1ccccc1', 'carbon-13': '[13CH3]O', 'radical': '[CH3]', 'oxygen radical': 'C[O]',
      'SMARTS': '[#6]~[#7]', 'Daylight SMARTS': '[NX3][CX3]=[OX1]', 'not kekulizable': 'c1cccc1', 'unreadable': 'C1CC(',
      'not a molecule': 'MALFORMED_INPUT_VALUE', 'over-long SMILES': 'C'.repeat(MAX_SMILES_LENGTH + 1),
    };
  };

  test('break it: conversions of hard molecules, Crux mode against RDKit mode', async () => {
    const molecules = breakers();
    const targets = [DG.chem.Notation.Smiles, 'cxsmiles' as DG.chem.Notation, DG.chem.Notation.MolBlock,
      DG.chem.Notation.V3KMolBlock, DG.chem.Notation.Smarts, 'cxsmarts' as DG.chem.Notation];
    const all = () => Object.fromEntries(Object.entries(molecules).map(([what, m]) => [what, [
      ...targets.map((t) => PackageFunctions.convertMolNotation(m, DG.chem.Notation.Unknown, t)),
      `${PackageFunctions.isSmiles(m)}|${PackageFunctions.validateMolecule(m)?.replace(/\[\d\d:\d\d:\d\d\] /g, '')}`]]));
    /** Whether RDKit reads `smarts` as a query its molecule `m` matches (a SMARTS written by RDKit 2026.03 can differ from
     * 2024.09's, [CH3] as [#6] |^1:0| where 2024.09 wrote [#6H3], the same query on that molecule). */
    const matchesItself = (smarts: string, m: string) => {
      const module = chemCommonRdKit.getRdKitModule();
      let mol = null;
      let query = null;
      try {
        mol = module.get_mol(m);
        query = module.get_qmol(smarts.replace(/ \|.*\|$/, ''));
        return !!mol && !!query && mol.get_substruct_match(query) !== '{}';
      } catch {
        return false;
      } finally {
        mol?.delete();
        query?.delete();
      }
    };
    // and a SMILES with a centre to a MOL block
    const run = async () => ({strings: all(), centre: PackageFunctions.convertMolNotation('C[C@H](N)C(=O)O',
      DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock)});
    const {rdkit: {strings: rdkit, centre: rdkitCentre}, crux: {strings: crux, centre}} = await inBothEngines(run);
    const differences: string[] = [];
    for (const what of Object.keys(molecules)) {
      crux[what].forEach((c: string, k: number) => {
        const r = rdkit[what][k];
        const t = targets[k] ?? 'isSmiles|validateMolecule';
        const query = t === DG.chem.Notation.Smarts || t === 'cxsmarts';
        const plain = (x: string) => x.replace(/ \|.*\|$/, ''); // a CXSMILES' molecule (its extension: below)
        const same = c === r || (k < targets.length && !query && rdkitCanonical(plain(c)) !== null &&
          rdkitCanonical(plain(c)) === rdkitCanonical(plain(r))) || (query && matchesItself(c, molecules[what]));
        if (!same)
          differences.push(`${what} -> ${t}: ${c.slice(0, 90)} vs RDKit's ${r.slice(0, 90)}`);
      });
    }
    expect(differences.join('\n'), '', 'the same molecules and answers');
    // what RDKit reads and Crux skips goes to RDKit: the same strings
    for (const what of ['query label', 'ring bond count'])
      expect(JSON.stringify(crux[what]), JSON.stringify(rdkit[what]), what);
    expect(crux['AND group'][1].includes('&1:'), true, 'the AND group kept in CXSMILES');
    expect(crux['AND group'][3].includes('STERAC'), true, 'the AND group kept in V3000');
    expect(centre.split('\n')[3].slice(12, 15), rdkitCentre.split('\n')[3].slice(12, 15),
      'a SMILES with a centre: the chiral flag of RDKit mode');
    expect(rdkitCanonical(centre), rdkitCanonical(rdkitCentre), 'the same molecule');
  });

  test('break it: a column\'s conversion across chunks, with rows Crux cannot read', async () => {
    const base = (await jobMolecules()).filter((m) => m && m !== 'C1CC(');
    const molecules = Array.from({length: 2500}, (_, i) => i === 1001 ? '' : base[i % base.length]);
    const unreadable = [0, 999, 1000, 2499];
    const convert = () => convertNotationForColumn(DG.Column.fromStrings('m', molecules), DG.chem.Notation.MolBlock);
    const rdkit = await inEngine(ChemicalEngine.RDKit, convert);
    const crux = await inEngine(ChemicalEngine.Crux, () => withCruxUnreadable(unreadable.map((i) => molecules[i]), convert));
    expect(crux.length, 2500);
    expect(crux[1001], '', 'an empty cell stays empty');
    // a row whose molecule Crux "cannot read" gets RDKit's block (each of them is also elsewhere in the column)
    for (const i of unreadable)
      expect(crux[i], rdkit[i], `row ${i}: RDKit's block`);
    const other = crux.filter((b, i) => !unreadable.includes(i) && i !== 1001 && !b.includes('M  END')).length;
    expect(other, 0, 'every other row a MOL block');
  }, {timeout: 300000});

  test('break it: a trap in the panels and conversions leaves RDKit for the session', async () => {
    const molecule = 'CC(=O)Oc1ccccc1C(=O)O';
    const rdkit = await inEngine(ChemicalEngine.RDKit, async () => [JSON.stringify(await getStructuralAlerts(molecule)),
      PackageFunctions.convertMolNotation(molecule, DG.chem.Notation.Unknown, DG.chem.Notation.MolBlock)]);
    await inEngine(ChemicalEngine.Crux, async () => {
      setCruxParseTrapForTests(true);
      try {
        const alerts = JSON.stringify(await getStructuralAlerts(molecule));
        expect(isCruxRequested(), false, 'the trap marked Crux failed');
        expect(alerts, rdkit[0], 'RDKit\'s alerts');
        expect(PackageFunctions.convertMolNotation(molecule, DG.chem.Notation.Unknown, DG.chem.Notation.MolBlock), rdkit[1],
          'RDKit\'s block');
        expect(PackageFunctions.isSmiles(molecule), true);
        expect(PackageFunctions.validateMolecule(molecule), '');
        let error = '';
        try {
          await getToxicityAlertsColumns(DG.Column.fromStrings('m', [molecule]), ['Brenk']);
        } catch (e: any) {
          error = e.message;
        }
        expect(error, TOXICITY_ALERTS_CRUX_FAILED, 'toxicity alerts name the failure');
      } finally {
        setCruxParseTrapForTests(false);
        resetCruxLoaderForTests(null); // the next session's Crux
      }
      expect(await cruxReady(), true, 'Crux back after the reset');
    });
  });

  test('break it: a Kekulé conversion leaves the RDKit search its molecules', async () => {
    // the RDKit pool's convertMolNotation kekulized the molecules its searches had cached (2026-10-23 review)
    await inEngine(ChemicalEngine.RDKit, async () => {
      const molecules = (await readDataframe('smiles.csv')).col('canonical_smiles')!;
      const before = (await chemSubstructureSearchLibrary(molecules, 'c1ccccc1', '')).trueCount;
      await convertNotationForColumn(molecules, DG.chem.Notation.Smiles, true);
      const after = (await chemSubstructureSearchLibrary(molecules, 'c1ccccc1', '')).trueCount;
      expect(after, before, 'the same hits after the Kekulé conversion');
    });
  });

  test('MMP: generated molecules and a pair\'s alignment on Crux', async () => {
    const cores = ['c1ccc(cc1)[*:1]', 'O=C(N[*:1])c1ccccc1', 'C1CCN(CC1)[*:1]', ''];
    const fragments = ['[*:1]C(F)(F)F', '[*:1]CCO', '[*:1]c1ccncc1', '[*:1]C'];
    const run = async () => ({linked: await mmpLinkFragmentsByEngine(cores, fragments),
      pair: await getInverseSubstructuresAndAlign(['c1ccc(cc1)[*:1]'], ['Cc1ccccc1'], ['Oc1ccccc1'],
        chemCommonRdKit.getRdKitModule())});
    // Crux did them
    const {rdkit, crux} = await inBothEngines(run, {noFallback: ['mmpLinkFragments', 'mmpAlign']});
    expect(JSON.stringify(crux.linked.map((s) => s && rdkitCanonical(s))),
      JSON.stringify(rdkit.linked.map((s) => s && rdkitCanonical(s))), 'the same molecules');
    expect(crux.linked[3], '', 'none for an empty core');
    expect(crux.pair.fromAligned[0].includes('M  END') && crux.pair.toAligned[0].includes('M  END'), true, 'aligned');
    expect(JSON.stringify([crux.pair.inverse1[0]?.atoms, crux.pair.inverse2[0]?.atoms]),
      JSON.stringify([rdkit.pair.inverse1[0]?.atoms, rdkit.pair.inverse2[0]?.atoms]), 'the atoms outside the core');
  });
});

/** The `scaffold` MOL blocks of a scaffold tree, depth first. */
function molblocksOf(tree: any): string[] {
  const blocks: string[] = [];
  const walk = (node: any) => {
    if (Array.isArray(node))
      node.forEach(walk);
    else if (node && typeof node === 'object') {
      if (typeof node.scaffold === 'string')
        blocks.push(node.scaffold);
      walk(node.child_nodes);
    }
  };
  walk(tree);
  return blocks;
}

/** Where two JSON values first differ, as a path. */
function firstDifference(a: any, b: any, path: string): string {
  if (typeof a !== typeof b || Array.isArray(a) !== Array.isArray(b) || a === null || b === null || typeof a !== 'object')
    return a === b ? '' : `${path || '/'}: ${JSON.stringify(a)?.slice(0, 80)} vs ${JSON.stringify(b)?.slice(0, 80)}`;
  const keys = Array.isArray(a) ? Array.from({length: Math.max(a.length, b.length)}, (_, i) => String(i)) :
    Array.from(new Set([...Object.keys(a), ...Object.keys(b)]));
  for (const k of keys) {
    const d = firstDifference(a[k], b[k], `${path}/${k}`);
    if (d)
      return d;
  }
  return '';
}
