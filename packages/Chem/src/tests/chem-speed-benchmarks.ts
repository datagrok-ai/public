/**
 * Benchmarks of Chem's functions end to end where users wait (package start, the grid's drawing, column detection,
 * the substructure and similarity searches, the column jobs), in both chemical engines, run interleaved in one page:
 * `grok test --benchmark --category "chem speed"`. Each test returns its medians (ms) as its result. Nothing of them
 * runs without `--benchmark`; `chem speed: same results` (every run) pins that the replaced code's results stay.
 */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';
import {PackageFunctions, getMoleculeRenderer} from '../package';
import {getMCS} from '../utils/most-common-subs';
import {getMapIdentifiers} from '../widgets/identifiers';
import {ensureContainerRunning} from '@datagrok-libraries/test/src/test-container-utils';
import {runTransformationReaction, setReactionYieldEvery50ForTests} from '../utils/reactions/reactions';
import {NamedReaction} from '../utils/reactions/types';
import {addDeprotectedColumn, defaultFragment, setDeprotectOnWorkersForTests} from '../analysis/deprotect';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {ChemicalEngine} from '../engine/chem-engine';
import {compileCrux, cruxModuleSync, cruxWasmUrl} from '../engine/crux-module';
import {getCruxService, resetCruxServiceForTests} from '../crux/crux-searches';
import {chemFindSimilar, chemGetFingerprints, chemGetDiversities, chemSubstructureSearchLibrary}
  from '../chem-searches';
import {chemDiversitySearch} from '../analysis/chem-diversity-viewer';
import {getDiverseSubset} from '@datagrok-libraries/utils/src/similarity-metrics';
import {tanimotoSimilarity} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import {diverseSubset} from '../utils/diverse-subset';
import {getMmpFrags, getMmpRulesCPU, mmpRulesOfPairs}
  from '../analysis/molecular-matched-pairs/mmp-analysis/mmpa-fragments';
import {getAllRulesOcasions} from '../analysis/molecular-matched-pairs/mmp-analysis/mmpa-differences';
import {MmpFragments, MmpRules, MolecularPair} from '../analysis/molecular-matched-pairs/mmp-analysis/mmpa-misc';
import {chemSimilarityBitset} from '../analysis/chem-similarity-viewer';
import {BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {Fingerprint} from '../utils/chem-common';
import {FILTER_SCAFFOLD_TAG} from '../constants';
import {getPropertiesAsColumns} from '../widgets/properties';
import {convertNotationForColumn} from '../utils/convert-notation-utils';
import {getStructuralAlertsByRules, RuleSet} from '../panels/structural-alerts';
import {getInchiKeysImpl, getInchisImpl} from '../panels/inchi';
import {convertMolNotationByEngine, fingerprintsByEngine, inchisByEngine} from '../engine/worker-jobs';
import {beautifyMolsByEngine, coordGenCoordsByEngine, flattenMoleculesByEngine} from '../engine/molecule-jobs';
import {rGroupsByEngine} from '../engine/analysis-jobs';
import {cruxInchi} from '../crux/crux-inchi';
import {loadCruxSketch} from '../crux/crux-sketcher';
import {inchiToInchiKey, inchiToSmiles, smilesToInchi, smilesToInchiKey} from '../python-free';
import {inEngine} from './engine-utils';
import {readDataframe} from './utils';

const ENGINES = [ChemicalEngine.Crux, ChemicalEngine.RDKit];
const REPS = 3;
/** ChEMBL's 100k molecules on the stand (the demo files); 10x smi10K where they are not there. */
const LARGE = 'System:DemoFiles/chem/chembl/chembl-100k.csv';

const median = (xs: number[]): number => {
  const s = [...xs].sort((a, b) => a - b);
  return s.length % 2 ? s[(s.length - 1) / 2] : (s[s.length / 2 - 1] + s[s.length / 2]) / 2;
};

async function ms(run: () => unknown): Promise<number> {
  const start = performance.now();
  await run();
  return performance.now() - start;
}

/** The medians (ms) of each leg, the legs run in turn `reps` times (interleaved); `prepare` before each leg. */
async function interleaved(legs: {[label: string]: () => unknown}, reps = REPS,
  prepare?: (label: string) => unknown): Promise<{[label: string]: number}> {
  const times: {[label: string]: number[]} = {};
  for (let rep = 0; rep < reps; rep++) {
    for (const [label, leg] of Object.entries(legs)) {
      await prepare?.(label);
      (times[label] ??= []).push(await ms(leg));
    }
  }
  const medians = Object.fromEntries(Object.entries(times).map(([label, t]) => [label, +median(t).toFixed(2)]));
  console.log(`chem speed: ${JSON.stringify(medians)} (all: ${JSON.stringify(times)})`);
  return medians;
}

/** `legs` per engine (`<label> <engine>`), each leg run in its engine. */
function byEngine(legs: {[label: string]: () => unknown}): {[label: string]: () => unknown} {
  const out: {[label: string]: () => unknown} = {};
  for (const engine of ENGINES) {
    for (const [label, leg] of Object.entries(legs))
      out[`${label} ${engine}`] = () => inEngine(engine, async () => leg());
  }
  return out;
}

let large: string | null = null;
async function largeCsv(): Promise<string> {
  if (large === null) {
    try {
      large = await grok.dapi.files.readAsText(LARGE);
    } catch {
      const small = (await readDataframe('tests/smi10K.csv')).toCsv();
      const rows = small.split('\n').slice(1).filter((l) => l !== '');
      large = ['smiles', ...Array.from({length: 10}, () => rows).flat()].join('\n');
    }
  }
  return large;
}

async function molTable(csv: string): Promise<DG.DataFrame> {
  const df = DG.DataFrame.fromCsv(csv);
  df.col('smiles')!.semType = DG.SEMTYPE.MOLECULE;
  return df;
}

category('chem speed', () => {
  let smi10k: string;

  before(async () => {
    if (!DG.Test.isInBenchmark)
      return;
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    smi10k = (await readDataframe('tests/smi10K.csv')).toCsv();
  });

  test('package start: wasm compile, Crux pool', async () => {
    if (!DG.Test.isInBenchmark)
      return;
    const rdkitUrl = `${chemCommonRdKit.getRdKitWebRoot()}/dist/RDKit_minimal_1.2.23.wasm`;
    const compile = (url: string) => WebAssembly.compileStreaming(fetch(url, {cache: 'no-store'}));
    const compiled = await interleaved({
      'compile Crux wasm': () => compile(cruxWasmUrl()),
      'compile RDKit wasm': () => compile(rdkitUrl),
    });
    await compileCrux();
    const pool = await interleaved({'Crux pool start': async () => {
      await resetCruxServiceForTests();
      await getCruxService();
    }});
    return JSON.stringify({...compiled, ...pool});
  }, {benchmark: true, timeout: 300000});

  test('grid drawing: new molecules, cache hits, filter highlight', async () => {
    if (!DG.Test.isInBenchmark)
      return;
    const df = await molTable(smi10k);
    const grid = df.plot.grid();
    const renderer = getMoleculeRenderer();
    const g = ui.canvas(200, 100).getContext('2d')!;
    const empty = () => {
      (renderer as any).cacheGeneration = -1;
      renderer._syncEngine();
    };
    const draw = (from: number, count: number) => {
      for (let i = from; i < from + count; i++) {
        const cell = grid.cell('smiles', i);
        renderer.render(g, 0, 0, 200, 100, cell, cell.style);
      }
    };
    const fresh = await interleaved(byEngine({'draw 300 new': () => draw(0, 300)}), REPS, () => empty());
    // cache hits: 40 cells drawn again, 25 times (scrolling back and forth over them)
    const hits = await interleaved(byEngine({'redraw 40 x 25': () => {
      draw(1000, 40);
      for (let k = 0; k < 25; k++)
        draw(1000, 40);
    }}));
    const col = df.col('smiles')!;
    col.temp[FILTER_SCAFFOLD_TAG] = JSON.stringify([{molecule: 'c1ccccc1', isSuperstructure: false, align: true}]);
    const highlighted = await interleaved(byEngine({'draw 150 new, filter highlight': () => draw(2000, 150)}), REPS,
      () => empty());
    delete col.temp[FILTER_SCAFFOLD_TAG];
    return JSON.stringify({...fresh, ...hits, ...highlighted});
  }, {benchmark: true, timeout: 600000});

  test('column detection on open', async () => {
    if (!DG.Test.isInBenchmark)
      return;
    const csv100k = await largeCsv();
    return JSON.stringify(await interleaved({
      'detect 10k': async () => grok.data.detectSemanticTypes(DG.DataFrame.fromCsv(smi10k)),
      'detect 100k': async () => grok.data.detectSemanticTypes(DG.DataFrame.fromCsv(csv100k)),
    }));
  }, {benchmark: true, timeout: 300000});

  test('similarity 100k: find similar, the viewer, diversity', async () => {
    if (!DG.Test.isInBenchmark)
      return;
    const csv = await largeCsv();
    const query = 'CC(=O)Oc1ccccc1C(=O)O';
    let df: DG.DataFrame;
    const first = await interleaved(byEngine({'find similar (new table)': () =>
      chemFindSimilar(df.col('smiles')!, query, {limit: 10, minScore: 0.3})}), REPS,
    async () => df = await molTable(csv));
    // the fingerprints of a new table: the engine's call alone, and Chem's function (the cache columns, BitArrays)
    const glue = await interleaved(byEngine({
      'fingerprints 100k, the workers\' call': () =>
        fingerprintsByEngine(Fingerprint.Morgan, df.col('smiles')!.toList(), false),
      'fingerprints 100k, chemGetFingerprints (new table)': () => chemGetFingerprints(df.col('smiles')!,
        Fingerprint.Morgan, false),
    }), REPS, async () => df = await molTable(csv));
    df = await molTable(csv);
    const col = df.col('smiles')!;
    // the cached legs per engine (each engine keeps its own fingerprint cache and drops the other's)
    let later: {[label: string]: number} = {};
    for (const engine of ENGINES) {
      later = {...later, ...await inEngine(engine, async () => {
        await chemGetFingerprints(col, Fingerprint.Morgan, false);
        const all = DG.BitSet.create(col.length).setAll(true);
        const medians = await interleaved({
          'find similar (fingerprints cached)': () => chemFindSimilar(col, query, {limit: 10, minScore: 0.3}),
          'similarity viewer: a row clicked': () =>
            chemSimilarityBitset(col, query, BitArrayMetricsNames.Tanimoto, 0.3, Fingerprint.Morgan),
          'cached fingerprints read': () => chemGetFingerprints(col, Fingerprint.Morgan, false),
          'diversity 10 (cached)': () => chemGetDiversities(col, 10),
          'diversity viewer 10 (cached)': () =>
            chemDiversitySearch(col, tanimotoSimilarity, 10, Fingerprint.Morgan, all),
        });
        return Object.fromEntries(Object.entries(medians).map(([k, v]) => [`${k} ${engine}`, v]));
      })};
    }
    // the diverse subset before (the utils library's getDiverseSubset) and now (utils/diverse-subset.ts): the same
    // picks for the same random draws, on the column's fingerprints
    const fps = await chemGetFingerprints(col, Fingerprint.Morgan, false);
    const indexes = fps.flatMap((f, i) => f && !f.allFalse ? [i] : []);
    const dist = (i1: number, i2: number) => 1 - tanimotoSimilarity(fps[indexes[i1]]!, fps[indexes[i2]]!);
    const random = Math.random;
    const picks: {[seed: number]: string[]} = {};
    let seed = 0;
    const seeded = () => {
      let state = ++seed * 7919;
      Math.random = () => (state = (state * 16807) % 2147483647) / 2147483647;
    };
    let diversity: {[label: string]: number};
    try {
      diversity = await interleaved({
        'diverse subset 10 of 100k, before': () => {
          seeded();
          (picks[seed] ??= []).push(getDiverseSubset(indexes.length, 10, dist).join(','));
        },
        'diverse subset 10 of 100k, now': () => {
          seed--; // the same draws as the leg before
          seeded();
          (picks[seed] ??= []).push(diverseSubset(indexes.length, 10, dist).join(','));
        },
      });
    } finally {
      Math.random = random;
    }
    const differing = Object.values(picks).filter((p) => p.length !== 2 || p[0] !== p[1]).length;
    return JSON.stringify({...first, ...glue, ...later, ...diversity,
      'diverse subsets compared': Object.keys(picks).length, 'diverse subsets differing': differing});
  }, {benchmark: true, timeout: 900000});

  test('substructure search 100k: first and later', async () => {
    if (!DG.Test.isInBenchmark)
      return;
    const csv = await largeCsv();
    const queries = ['c1ccncc1', 'C(=O)N', '[#6]1:[#6]:[#7]:[#6]:[#6]:1'];
    let df: DG.DataFrame;
    const first = await interleaved(byEngine({'first search (new table)': () =>
      chemSubstructureSearchLibrary(df.col('smiles')!, 'c1ccccc1', '')}), REPS,
    async () => df = await molTable(csv));
    df = await molTable(csv);
    for (const engine of ENGINES) // the index (Crux) or the pattern fingerprints (RDKit) built
      await inEngine(engine, () => chemSubstructureSearchLibrary(df.col('smiles')!, 'c1ccccc1', ''));
    const later = await interleaved(byEngine({'3 later searches': async () => {
      for (const q of queries)
        await chemSubstructureSearchLibrary(df.col('smiles')!, q, '');
    }}));
    return JSON.stringify({...first, ...later});
  }, {benchmark: true, timeout: 900000});

  test('column jobs 10k: properties, notation, alerts, InChI', async () => {
    if (!DG.Test.isInBenchmark)
      return;
    const df = await molTable(smi10k);
    const col = df.col('smiles')!;
    const pains = {PAINS: true} as unknown as RuleSet;
    await getCruxService();
    return JSON.stringify(await interleaved(byEngine({
      'properties (all 10)': () => getPropertiesAsColumns(col, ['MW', 'HBA', 'HBD', 'LogP', 'LogS', 'PSA',
        'Rotatable bonds', 'Stereo centers', 'Molecule charge', 'Molecular formula']),
      'to molblock': () => convertNotationForColumn(col, DG.chem.Notation.MolBlock),
      'to canonical smiles': () => convertNotationForColumn(col, DG.chem.Notation.Smiles),
      'PAINS alerts': () => getStructuralAlertsByRules(col, pains),
      'InChI': () => getInchisImpl(col),
      'InChI keys': () => getInchiKeysImpl(col),
    })));
  }, {benchmark: true, timeout: 900000});

  test('InChI 10k: column vs molecule by molecule, first call', async () => {
    if (!DG.Test.isInBenchmark)
      return;
    const molecules = Array.from((await molTable(smi10k)).col('smiles')!.toList() as (string | null)[],
      (m) => m ?? '');
    return inEngine(ChemicalEngine.Crux, async () => {
      // a fresh pool: the first InChI job of each worker loads the InChI module
      await resetCruxServiceForTests();
      await getCruxService();
      const first1k = await ms(() => inchisByEngine(molecules.slice(0, 1000), false));
      const [inchis, keys] = [await inchisByEngine(molecules, false), await inchisByEngine(molecules, true)];
      // byte identity: each molecule's InChI and key by the main thread's module, as the column's were written before
      // (the input by Crux, the library by crux-sketch's InChI module)
      let differ = 0;
      const lib = await loadCruxSketch();
      const crux = cruxModuleSync() as any;
      const inputs: (Uint8Array | null)[] = crux.chemInchiInputs(molecules);
      const ids = {inchi: await smilesToInchi(molecules), key: await smilesToInchiKey(molecules)};
      let differIds = 0;
      for (let i = 0; i < molecules.length; i++) {
        const one = molecules[i] ? await cruxInchi(molecules[i]) : {inchi: '', inchiKey: ''};
        if ((one?.inchi ?? null) !== inchis[i] || (one?.inchiKey ?? null) !== keys[i])
          differ++;
        const inchi = inputs[i] ? await lib.inchiFromInput(inputs[i]!) : '';
        const key = inchi ? await lib.inchiKeyFromInchi(inchi) : '';
        if (inchi !== ids.inchi[i] || key !== ids.key[i])
          differIds++;
      }
      // Map Identifiers from InChI: on the workers now, crux-sketch's module on the main thread before
      const fromInchi = {smiles: await inchiToSmiles(ids.inchi), key: await inchiToInchiKey(ids.inchi)};
      let differFromInchi = 0;
      const mainThreadFromInchi = async () => {
        let differs = 0;
        for (let i = 0; i < ids.inchi.length; i++) {
          const inchi = ids.inchi[i];
          const smiles = inchi ? await lib.inchiToMol(inchi).catch(() => '') : '';
          const key = inchi ? (await lib.inchiKeyFromInchi(inchi)) || null : null;
          if (smiles !== fromInchi.smiles[i] || key !== fromInchi.key[i])
            differs++;
        }
        return differs;
      };
      differFromInchi = await mainThreadFromInchi();
      const times = await interleaved({
        'InChI 10k (column)': () => inchisByEngine(molecules, false),
        'InChI keys 10k (column)': () => inchisByEngine(molecules, true),
        'Map Identifiers smiles to inchi 10k': () => smilesToInchi(molecules),
        'Map Identifiers inchi to smiles 10k, workers (now)': () => inchiToSmiles(ids.inchi),
        'Map Identifiers inchi to smiles and key 10k, main thread (before)': () => mainThreadFromInchi(),
        'Map Identifiers inchi to smiles and key 10k, workers (now)': async () => {
          await inchiToSmiles(ids.inchi);
          await inchiToInchiKey(ids.inchi);
        },
      });
      return JSON.stringify({'first InChI 1k after pool start': +first1k.toFixed(1), ...times,
        'molecules': molecules.length, 'InChI or key differing from the molecule path': differ,
        'Map Identifiers InChI or key differing from the main-thread library': differIds,
        'Map Identifiers SMILES or key from InChI differing from the main-thread library': differFromInchi});
    });
  }, {benchmark: true, timeout: 900000});

  test('Crux workers without pauses: conversions, coordinates, R-groups', async () => {
    if (!DG.Test.isInBenchmark)
      return;
    const molecules = Array.from((await molTable(smi10k)).col('smiles')!.toList() as (string | null)[],
      (m) => m ?? '');
    const sar = (await readDataframe('sar-small.csv')).col('smiles')!.toList() as string[];
    const series = Array.from({length: 5}, () => sar).flat(); // 1,000 molecules of one core
    const core = 'O=C1CN=C(c2ccccc2N1)C1CCCCC1';
    const options = JSON.stringify({matchingStrategy: 'Greedy', includeTargetMolInResults: true,
      onlyMatchAtRGroups: false});
    const service = await inEngine(ChemicalEngine.Crux, () => getCruxService());
    const outputs: {[leg: string]: string} = {};
    const legs: {[label: string]: () => Promise<unknown>} = {
      'to molblock 10k': () => convertMolNotationByEngine(molecules, DG.chem.Notation.MolBlock, false),
      'CoordGen 10k': () => coordGenCoordsByEngine(molecules),
      'flatten 10k': () => flattenMoleculesByEngine(molecules),
      'beautifyMols 2k': () => beautifyMolsByEngine(molecules.slice(0, 2000)),
      'R-groups 1k of one core': () => rGroupsByEngine(series, core, false, options),
    };
    const ab: {[label: string]: () => Promise<void>} = {};
    for (const [label, run] of Object.entries(legs)) {
      for (const when of ['pauses every 50 (before)', 'no pauses (now)']) {
        ab[`${label}, ${when}`] = () => inEngine(ChemicalEngine.Crux, async () => {
          outputs[`${label}, ${when}`] = JSON.stringify(await run());
        });
      }
    }
    let times: {[label: string]: number};
    try {
      times = await interleaved(ab, REPS, (label) => service.pausesForTests(label.endsWith('(before)') ? 50 : null));
    } finally {
      await service.pausesForTests(null);
    }
    const differing = Object.keys(legs).filter((label) =>
      outputs[`${label}, pauses every 50 (before)`] !== outputs[`${label}, no pauses (now)`]);
    return JSON.stringify({...times, 'outputs differing': differing});
  }, {benchmark: true, timeout: 900000});

  test('reactions 10k: transformation, deprotection', async () => {
    if (!DG.Test.isInBenchmark)
      return;
    const molecules = Array.from((await molTable(smi10k)).col('smiles')!.toList() as (string | null)[],
      (m) => m ?? '');
    const reactions: NamedReaction[] = JSON.parse(await _package.files.readAsText('reactions/default-reactions.json'));
    const hydrolysis = reactions.find((r) => r.name === 'Ester Hydrolysis (Saponification)')!.reactionSmarts;
    const rdkit = PackageFunctions.getRdKitModule();
    const products: {[leg: string]: string} = {};
    const legs: {[label: string]: () => Promise<void>} = {};
    for (const engine of ENGINES) {
      for (const every50 of [true, false]) {
        const label = `ester hydrolysis 10k, ${every50 ? 'yield every 50 (before)' : 'yield by time (now)'} ${engine}`;
        legs[label] = () => inEngine(engine, async () => {
          setReactionYieldEvery50ForTests(every50);
          try {
            products[label] = JSON.stringify((await runTransformationReaction(rdkit, hydrolysis, molecules)).products);
          } finally {
            setReactionYieldEvery50ForTests(false);
          }
        });
      }
    }
    const times = await interleaved(legs);
    const differing = ENGINES.filter((e) => products[`ester hydrolysis 10k, yield every 50 (before) ${e}`] !==
      products[`ester hydrolysis 10k, yield by time (now) ${e}`]);
    // Deprotect: the column on the Crux workers (now) or on the main thread (before), and RDKit mode's
    const df = await molTable(smi10k);
    const deprotected: {[leg: string]: string} = {};
    const deprotectLeg = (engine: ChemicalEngine, workers: boolean, label: string) => () => inEngine(engine,
      async () => {
        setDeprotectOnWorkersForTests(workers);
        try {
          const col = await addDeprotectedColumn(df, df.col('smiles')!, defaultFragment, rdkit);
          deprotected[label] = JSON.stringify(col!.toList());
          df.columns.remove(col!.name);
        } finally {
          setDeprotectOnWorkersForTests(true);
        }
      });
    const deprotect = await interleaved({
      'deprotect 10k, main thread (before) Crux': deprotectLeg(ChemicalEngine.Crux, false, 'before'),
      'deprotect 10k, Crux workers (now) Crux': deprotectLeg(ChemicalEngine.Crux, true, 'now'),
      'deprotect 10k RDKit': deprotectLeg(ChemicalEngine.RDKit, true, 'rdkit'),
    });
    return JSON.stringify({...times, ...deprotect, 'engines whose products differ': differing,
      'deprotected columns identical (before, now)': deprotected.before === deprotected.now});
  }, {benchmark: true, timeout: 900000});

  test('Map Identifiers step by step', async () => {
    // calculate's "map identifiers" test, step by step, each within its own 60 s: the chem-chem container (a publish
    // leaves it stopped), Crux's InChIs and keys on the workers, the ChEMBL and UniChem lookups
    if (!DG.Test.isInBenchmark)
      return;
    const log: string[] = [];
    const step = async (name: string, run: () => Promise<unknown>) => {
      const start = performance.now();
      let timer: ReturnType<typeof setTimeout> | undefined;
      try {
        const late = new Promise((_, reject) =>
          timer = setTimeout(() => reject(new Error('no answer in 60 s')), 60000));
        await Promise.race([run(), late]);
        log.push(`${name}: ${(performance.now() - start).toFixed(0)} ms`);
      } catch (e) {
        const message = e instanceof Error ? e.message : String(e);
        log.push(`${name}: ${(performance.now() - start).toFixed(0)} ms, ${message}`);
      } finally {
        clearTimeout(timer);
      }
      console.log(`chem speed: ${log[log.length - 1]}`);
    };
    return inEngine(ChemicalEngine.Crux, async () => {
      // the test's order: the container first (a publish leaves it stopped)
      await step('ensureContainerRunning chem-chem (first)', () => ensureContainerRunning('name = "chem-chem"', 90000));
      const df = grok.data.demo.molecules(20);
      await grok.data.detectSemanticTypes(df);
      const mols = df.col('smiles')!.toList();
      await step('smilesToInchi (Crux workers)', () => smilesToInchi(mols));
      await step('smilesToInchiKey (Crux workers)', () => smilesToInchiKey(mols));
      await step('getMapIdentifiers to inchi', () => getMapIdentifiers(df, df.col('smiles')!, 'smiles', 'inchi'));
      for (const to of ['chembl', 'actor', 'pubchem'])
        await step(`getMapIdentifiers to ${to}`, () => getMapIdentifiers(df, df.col('smiles')!, 'smiles', to));
      await step('ensureContainerRunning chem-chem', () => ensureContainerRunning('name = "chem-chem"', 90000));
      await step('getMapIdentifiers to chembl, again',
        () => getMapIdentifiers(df, df.col('smiles')!, 'smiles', 'chembl'));
      return log.join('; ');
    });
  }, {benchmark: true, timeout: 900000});

  test('analyses: scaffold tree, MCS', async () => {
    if (!DG.Test.isInBenchmark)
      return;
    const all = Array.from((await molTable(smi10k)).col('smiles')!.toList() as (string | null)[], (m) => m ?? '')
      .filter((m) => m !== '');
    const table = (n: number) => {
      const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('smiles', all.slice(0, n))]);
      df.col('smiles')!.semType = DG.SEMTYPE.MOLECULE;
      df.col('smiles')!.meta.units = DG.UNITS.Molecule.SMILES;
      return df;
    };
    const [t1k, t5k] = [table(1000), table(5000)];
    // the MCS of a series (the R-group dialog's, Cluster MCS'): the MMP demo's first 500 molecules
    const demo = await grok.dapi.files.readAsText('System:DemoFiles/chem/mmp_demo.csv');
    const series = DG.Column.fromStrings('smiles',
      demo.split('\n').slice(1).map((l) => l.split('\t')[0]).filter((m) => m !== '').slice(0, 500));
    const crux = await inEngine(ChemicalEngine.Crux, () => interleaved({
      'scaffold tree 1k (Crux)': () => PackageFunctions.getScaffoldTree(t1k, 10, false),
      'scaffold tree 5k (Crux)': () => PackageFunctions.getScaffoldTree(t5k, 10, false),
    }, 2));
    const mcs = await interleaved(byEngine({'MCS of 500': () => getMCS(series, true, true)}), 2);
    return JSON.stringify({...crux, ...mcs});
  }, {benchmark: true, timeout: 900000});

  test('MMP: fragments, rules, rule occasions', async () => {
    if (!DG.Test.isInBenchmark)
      return;
    // the MMP demo's molecules (ChEMBL series: many pairs), the first 3,000
    const demo = await grok.dapi.files.readAsText('System:DemoFiles/chem/mmp_demo.csv');
    const molecules = demo.split('\n').slice(1).map((l) => l.split('\t')[0]).filter((m) => m !== '').slice(0, 3000);
    const fragments = await interleaved(byEngine({'MMP fragments 3k': () => getMmpFrags(molecules)}));
    const [frags] = await inEngine(ChemicalEngine.Crux, () => getMmpFrags(molecules));
    const rules: {[leg: string]: string} = {};
    const occasions: {[leg: string]: string} = {};
    const times = await interleaved({
      'MMP rules (CPU) 3k, before': () => {
        rules.before = JSON.stringify(mmpRulesCpuBefore(frags, 0.4));
      },
      'MMP rules (CPU) 3k, now': () => {
        rules.now = JSON.stringify(getMmpRulesCPU(frags, 0.4));
      },
    });
    const [mmpRules] = getMmpRulesCPU(frags, 0.4);
    const occ = await interleaved({
      'MMP rule occasions, before': () => {
        occasions.before = JSON.stringify(allRulesOcasionsBefore(mmpRules, frags, {fragmentIdxs: [], frequencies: []}));
      },
      'MMP rule occasions, now': () => {
        const sorting = {fragmentIdxs: [] as number[], frequencies: [] as number[]};
        occasions.now = JSON.stringify([getAllRulesOcasions(mmpRules, frags, sorting as any), sorting]);
      },
    });
    // the rules of the pairs alone (the GPU path's work once the GPU has found the pairs), on the whole MMP demo
    const allDemo = demo.split('\n').slice(1).map((l) => l.split('\t')[0]).filter((m) => m !== '');
    const [demoFrags] = await inEngine(ChemicalEngine.Crux, () => getMmpFrags(allDemo));
    const pairs = mmpPairsBefore(demoFrags, 0.4);
    const arrays = ['firstR', 'first', 'secondR', 'second', 'core'].map((k) =>
      Int32Array.from(pairs, (p) => (p as any)[k] as number));
    const fill: {[leg: string]: string} = {};
    const fillTimes = await interleaved({
      'MMP rules of the demo\'s pairs (GPU path), before': () => {
        const rules: MmpRules = {rules: [], smilesFrags: []};
        let [ruleCounter, cases] = [0, 0];
        for (let i = 0; i < pairs.length; i++) {
          [ruleCounter, cases] = fillRulesBefore(rules, arrays[0][i], arrays[1][i], arrays[2][i], arrays[3][i],
            arrays[4][i], ruleCounter, cases);
        }
        fill.before = JSON.stringify([rules, cases]);
      },
      'MMP rules of the demo\'s pairs (GPU path), now': () => {
        fill.now = JSON.stringify(mmpRulesOfPairs(arrays[0], arrays[1], arrays[2], arrays[3], arrays[4]));
      },
    });
    return JSON.stringify({...fragments, ...times, ...occ, ...fillTimes, 'demo molecules': allDemo.length,
      'demo pairs': pairs.length, 'demo rules identical': fill.before === fill.now, 'rules': mmpRules.rules.length,
      'pairs': mmpRules.rules.reduce((n, r) => n + r.pairs.length, 0), 'rules identical': rules.before === rules.now,
      'occasions identical': occasions.before === occasions.now});
  }, {benchmark: true, timeout: 900000});
});

// The speed work's equivalences, cheap enough for every run: the new code gives what the code it replaced gave.
category('chem speed: same results', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('diverse subset: the picks of the library\'s getDiverseSubset', async () => {
    const df = await readDataframe('smiles.csv');
    const fps = await chemGetFingerprints(df.col('canonical_smiles')!, Fingerprint.Morgan, false);
    // 200 molecules: the library's search is cubic in the subset's size (a whole subset of 1,000 takes a minute)
    const indexes = fps.flatMap((f, i) => f && !f.allFalse ? [i] : []).slice(0, 200);
    const dist = (i1: number, i2: number) => 1 - tanimotoSimilarity(fps[indexes[i1]]!, fps[indexes[i2]]!);
    const random = Math.random;
    try {
      for (const [seed, n] of [[1, 10], [2, 25], [3, 1], [4, indexes.length]]) {
        const picks = [getDiverseSubset, diverseSubset].map((pick) => {
          let state = seed * 7919;
          Math.random = () => (state = (state * 16807) % 2147483647) / 2147483647;
          return pick(indexes.length, n, dist).join(',');
        });
        expect(picks[1], picks[0], `seed ${seed}, ${n} of ${indexes.length}`);
      }
    } finally {
      Math.random = random;
    }
  });

  test('MMP rules and occasions: the arrays\' searches\' results', async () => {
    const demo = await grok.dapi.files.readAsText('System:DemoFiles/chem/mmp_demo.csv');
    const molecules = demo.split('\n').slice(1).map((l) => l.split('\t')[0]).filter((m) => m !== '').slice(0, 400);
    const [frags] = await getMmpFrags(molecules);
    const before = mmpRulesCpuBefore(frags, 0.4);
    const now = getMmpRulesCPU(frags, 0.4);
    expect(now[0].rules.length > 0, true, 'some rules');
    expect(JSON.stringify(now), JSON.stringify(before), 'the CPU rules');
    const pairs = mmpPairsBefore(frags, 0.4);
    const arrays = ['firstR', 'first', 'secondR', 'second', 'core'].map((k) =>
      Int32Array.from(pairs, (p) => (p as any)[k] as number));
    expect(JSON.stringify(mmpRulesOfPairs(arrays[0], arrays[1], arrays[2], arrays[3], arrays[4])),
      JSON.stringify([before[0], before[1]]), 'the rules of the pairs');
    const sorting = {fragmentIdxs: [] as number[], frequencies: [] as number[]};
    expect(JSON.stringify([getAllRulesOcasions(now[0], frags, sorting as any), sorting]),
      JSON.stringify(allRulesOcasionsBefore(before[0], frags, {fragmentIdxs: [], frequencies: []})), 'the occasions');
  });
});

// The MMP rules as they were built before (mmpa-fragments.ts' getMmpRulesCPU and fillRules searching the arrays per
// pair), the reference the new build is compared with.
function mmpRulesCpuBefore(fragsOut: MmpFragments, fragmentCutoff: number): [MmpRules, number, boolean] {
  const mmpRules: MmpRules = {rules: [], smilesFrags: []};
  const pairs = mmpPairsBefore(fragsOut, fragmentCutoff);
  let ruleCounter = 0;
  let allCasesCounter = 0;
  for (let i = 0; i < pairs.length; i++) {
    [ruleCounter, allCasesCounter] = fillRulesBefore(mmpRules, pairs[i].firstR, pairs[i].first, pairs[i].secondR,
      pairs[i].second, pairs[i].core, ruleCounter, allCasesCounter);
  }
  return [mmpRules, allCasesCounter, false];
}

/** getMmpRulesCPU's pairs (unchanged): each pair of molecules' best common core within the cutoff. */
function mmpPairsBefore(fragsOut: MmpFragments, fragmentCutoff: number): MolecularPair[] {
  const dim = fragsOut.fragCodes.length;
  const pairs: MolecularPair[] = [];
  for (let i = 0; i < dim; i++) {
    const dimFirstMolecule = fragsOut.fragCodes[i].length;
    for (let j = i + 1; j < dim; j++) {
      const dimSecondMolecule = fragsOut.fragCodes[j].length;
      let core: number | null = null;
      let r1: number | null = null;
      let r2: number | null = null;
      for (let p1 = 0; p1 < dimFirstMolecule; p1++) {
        for (let p2 = 0; p2 < dimSecondMolecule; p2++) {
          if (fragsOut.fragCodes[i][p1][0] === fragsOut.fragCodes[j][p2][0]) {
            const newCore = fragsOut.fragCodes[i][p1][0];
            if (core === null || fragsOut.sizes[newCore] > fragsOut.sizes[core]) {
              core = newCore;
              r1 = fragsOut.fragCodes[i][p1][1];
              r2 = fragsOut.fragCodes[j][p2][1];
            }
          }
        }
      }
      if (core === null)
        continue;
      if (fragsOut.idToName[core] === '' || fragsOut.sizes[r1!] / fragsOut.sizes[core] > fragmentCutoff ||
        fragsOut.sizes[r2!] / fragsOut.sizes[core] > fragmentCutoff)
        continue;
      pairs.push({first: i, second: j, core, firstR: r1!, secondR: r2!});
    }
  }
  return pairs;
}

function fillRulesBefore(mmpRules: MmpRules, rFirst: number, idxFirst: number, rSecond: number, idxSecond: number,
  core: number, ruleCounter: number, allCasesCounter: number): [number, number] {
  let ruleSmiles1 = mmpRules.smilesFrags.indexOf(rFirst);
  let ruleSmiles2 = mmpRules.smilesFrags.indexOf(rSecond);
  let ruleIndexStraight = -1;
  let ruleIndexInverse = -1;
  for (let ind = 0; ind < mmpRules.rules.length; ind++) {
    if (mmpRules.rules[ind].sr1 === ruleSmiles1 && mmpRules.rules[ind].sr2 === ruleSmiles2)
      ruleIndexStraight = ind;
    if (mmpRules.rules[ind].sr1 === ruleSmiles2 && mmpRules.rules[ind].sr2 === ruleSmiles1)
      ruleIndexInverse = ind;
  }
  if (ruleSmiles1 === -1) {
    mmpRules.smilesFrags.push(rFirst);
    ruleSmiles1 = mmpRules.smilesFrags.length -1;
  }
  if (ruleSmiles2 === -1) {
    mmpRules.smilesFrags.push(rSecond);
    ruleSmiles2 = mmpRules.smilesFrags.length -1;
  }
  const indxFirst = ruleSmiles1 < ruleSmiles2;
  if (ruleIndexStraight === -1) {
    mmpRules.rules.push({sr1: indxFirst ? ruleSmiles1: ruleSmiles2, sr2: indxFirst ? ruleSmiles2: ruleSmiles1,
      pairs: []});
    mmpRules.rules[ruleCounter].pairs.push({fs: indxFirst ? idxFirst : idxSecond,
      ss: indxFirst ? idxSecond : idxFirst, core});
    ruleCounter++;
    mmpRules.rules.push({sr1: indxFirst ? ruleSmiles2: ruleSmiles1, sr2: indxFirst ? ruleSmiles1: ruleSmiles2,
      pairs: []});
    mmpRules.rules[ruleCounter].pairs.push({fs: indxFirst ? idxSecond : idxFirst,
      ss: indxFirst ? idxFirst : idxSecond, core});
    ruleCounter++;
    allCasesCounter += 2;
  } else {
    mmpRules.rules[ruleIndexStraight].pairs.push({fs: idxFirst, ss: idxSecond, core});
    mmpRules.rules[ruleIndexInverse].pairs.push({fs: idxSecond, ss: idxFirst, core});
    allCasesCounter += 2;
  }
  return [ruleCounter, allCasesCounter];
}

/** mmpa-differences.ts' getAllRulesOcasions as it was (an indexOf per rule), and the sorting info it filled. */
function allRulesOcasionsBefore(mmpr: MmpRules, frags: MmpFragments,
  fragSortingInfo: {fragmentIdxs: number[], frequencies: number[]}): unknown {
  const allSize = mmpr.rules.length;
  const fromFrag = new Array<string>(allSize);
  const toFrag = new Array<string>(allSize);
  const occasions = new Int32Array(allSize);
  for (let i = 0; i < allSize; i++) {
    const fromFragment = mmpr.smilesFrags[mmpr.rules[i].sr1];
    fromFrag[i] = frags.idToName[fromFragment];
    toFrag[i] = frags.idToName[mmpr.smilesFrags[mmpr.rules[i].sr2]];
    occasions[i] = mmpr.rules[i].pairs.length;
    const fragIdxInSortingInfo = fragSortingInfo.fragmentIdxs.indexOf(fromFragment);
    if (fragIdxInSortingInfo === -1) {
      fragSortingInfo.fragmentIdxs.push(fromFragment);
      fragSortingInfo.frequencies.push(mmpr.rules[i].pairs.length);
    } else
      fragSortingInfo.frequencies[fragIdxInSortingInfo] += mmpr.rules[i].pairs.length;
  }
  return [[fromFrag, toFrag, occasions], fragSortingInfo];
}
