import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';
import {PackageFunctions} from '../package';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {readDataframe} from './utils';
import {ChemicalEngine, cruxReady, setChemicalEngine} from '../engine/chem-engine';
import {cruxFallbackCount} from '../engine/fallback';
import {setCruxUnreadableForTests} from '../engine/worker-jobs';
import {getCruxService} from '../crux/crux-searches';
import {getMmpFrags} from '../analysis/molecular-matched-pairs/mmp-analysis/mmpa-fragments';
import {MmpFragments} from '../analysis/molecular-matched-pairs/mmp-analysis/mmpa-misc';
import {rGroupsMinilib} from '../analysis/r-group-analysis';
import {getMCS} from '../utils/most-common-subs';
import {decomposeClusters} from '../analysis/sar-matrix/sar-matrix-decompose';
import {linkStaged} from '../analysis/sar-matrix/sar-matrix-assemble';
import {LinkStages} from '../analysis/sar-matrix/sar-matrix-link';

/**
 * MMPA fragmentation, R-group decomposition and FMCS in Crux mode ("Chemical engine"): engine/analysis-jobs.ts. Each
 * feature's result comes from the Crux workers (no RDKit worker asked, no fallback) and equals RDKit mode's on Chem's
 * test data, or RDKit 2026.03.6's where Chem's RDKit 2024.09 differs (said where).
 * Fragment linking and molzip: in Crux mode the SAR matrix's joins by molzip, as the string assembly in RDKit mode;
 * molzip's expected values are the pinned RDKit 2026.03.6's (Python `Chem.molzip`, the SMILES read again).
 */

/** How many times the RDKit pool was asked for (each `getRdKitService()` counts). */
const rdkitPoolAsks = () => chemCommonRdKit._rdKitService?.timesInitialized ?? 0;

/** `run()` in RDKit mode, then in Crux mode: both results; in Crux mode no RDKit worker asked and no fallback to
 * RDKit for `calls`. */
async function inBothEngines<T>(run: () => Promise<T>, calls: string[]): Promise<{rdkit: T, crux: T}> {
  try {
    setChemicalEngine(ChemicalEngine.RDKit);
    const rdkit = await run();
    setChemicalEngine(ChemicalEngine.Crux);
    expect(await cruxReady(), true, 'the Crux module loads');
    await getCruxService();
    const [asks, fallbacks] = [rdkitPoolAsks(), calls.map((c) => cruxFallbackCount(c))];
    const crux = await run();
    expect(calls.map((c) => cruxFallbackCount(c)).join(','), fallbacks.join(','), `RDKit did ${calls} for Crux`);
    expect(rdkitPoolAsks(), asks, 'the RDKit workers were asked for');
    return {rdkit, crux};
  } finally {
    setChemicalEngine(null);
  }
}

/** Each molecule's fragment pairs by name. */
function fragmentNames(frags: MmpFragments): string[][][] {
  return frags.fragCodes.map((pairs) => pairs.map(([a, b]) => [frags.idToName[a], frags.idToName[b]]));
}

const column = (name: string, values: string[]) => DG.Column.fromStrings(name, values);

/** The R-group test molecules (top menu r-groups' sample table, the core benzene). */
const R_GROUP_SAMPLE = ['O=C1CN=C(c2ccccc2N1)C3CCCCC3', 'CN1C(=O)CN=C(c2ccccc12)C3CCCCC3',
  'CCCCN1C(=O)CN=C(c2ccccc12)C3CCCCC3', 'CC(C)CCN1C(=O)CN=C(c2ccccc12)C3CCCCC3',
  'O=C1CN=C(c2ccccc2N1C)C3CCCCC3', 'O=C1CN=C(c2cc(Cl)ccc2N1)C3CCCCC3', 'CN1C(=O)CN=C(c2cc(Cl)ccc12)C3CCCCC3'];

/** Chem's four MCS comparisons: [exact atoms, exact bonds]. */
const MCS_COMPARISONS: [boolean, boolean][] = [[true, true], [false, true], [true, false], [false, false]];

/** The pinned RDKit 2026.03.6's FindMCS SMARTS of the SAR table's 17 molecules and SPGI-100's first 25, for
 * {@link MCS_COMPARISONS} (scratch mcs_expected.py). */
const MCS_R26 = {
  sar: new Array<string>(4).fill('[#8]=[#6]1-[#6]-[#7]=[#6](-[#6]2:[#6](-[#7]-1):[#6]:[#6]:[#6]:[#6]:2)-[#6]1-[#6]-' +
    '[#6]-[#6]-[#6]-[#6]-1'),
  spgi: ['[#6]-[#6]-[#7](-[#6])-[#6]', '[#6,#7,#9]-[#6,#7]-[#7,#6](-[#6,#8,#9])-[#6,#7,#9]',
    '[#6]-,:[#6]-,:[#7](-[#6]-[#6])-,:[#6]',
    '[#6]-,:[#6,#7]-,:[#7,#6](-,:[#6,#7,#16]-,:[#6,#7]-,:[#6,#7]-,:[#7,#6]-,:[#6,#7,#8]-,:[#6,#7])-,:[#6,#7]-,:' +
      '[#6,#7]-,:[#6,#7]-,:[#6,#7,#9]'],
};

category('crux analyses', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  after(async () => {
    setCruxUnreadableForTests(null);
    setChemicalEngine(null);
  });

  test('MMPA fragments: Crux mode on the Crux workers, as RDKit mode', async () => {
    // the MMP test table: the same fragments and SMILES in both engines; and two molecules fragmentMol does not cut
    // (no cuttable bond; 24 of them, over its 20), which RDKit's code gives no SMILES on both engines: no RDKit asked
    const molecules = [...(await readDataframe('demo_files/matched_molecular_pairs.csv')).col('smiles')!.toList(), '',
      'c1ccccc1', 'C'.repeat(25)].map((m) => m ?? '') as string[];
    const {rdkit, crux} = await inBothEngines(() => getMmpFrags(molecules), ['MMPA fragments']);
    expect(JSON.stringify(crux[1]), JSON.stringify(rdkit[1]), 'the canonical SMILES');
    const [a, b] = [fragmentNames(rdkit[0]), fragmentNames(crux[0])];
    const differ = a.flatMap((pairs, i) => JSON.stringify(pairs) === JSON.stringify(b[i]) ? [] : [molecules[i]]);
    expect(differ.join(' '), '', 'molecules whose fragments differ');
    expect(a.reduce((n, pairs) => n + pairs.length, 0), 170, 'the cuts of the 38 molecules');
  });

  test('MMPA fragments: Crux mode keeps the stereo at a cut, as RDKit 2026.03', async () => {
    // RDKit.js 2026.03.6's fragments (scratch mmp-expected.mjs). Chem's RDKit 2024.09 inverts the centre next to a cut
    // of the first (COC(=O)[C@@H](C)[*:1]) and drops the E/Z of the second's cut double bond (FC=C[*:1]): RDKit's own
    // fixes between the releases, as on 5,617 of the MMP demo's 20,267 molecules
    const expected = [
      [['COC(=O)[C@@H](N)[*:1]', 'C[*:1]'], ['COC(=O)[C@H](C)[*:1]', 'N[*:1]'], ['C[C@H](N)[*:1]', 'COC(=O)[*:1]'],
        ['C[C@H](N)C(=O)O[*:1]', 'C[*:1]']],
      [['C(=C\\[*:1])/c1ccccc1', 'F[*:1]'], ['c1ccc([*:1])cc1', 'F/C=C/[*:1]']]];
    setChemicalEngine(ChemicalEngine.Crux);
    try {
      expect(await cruxReady(), true);
      const [frags] = await getMmpFrags(['C[C@H](N)C(=O)OC', 'F/C=C/c1ccccc1']);
      expect(JSON.stringify(fragmentNames(frags)), JSON.stringify(expected));
    } finally {
      setChemicalEngine(null);
    }
  });

  test('MMPA fragments: RDKit fragments a molecule Crux cannot read', async () => {
    const molecules = ['c12ccc(Cl)cc2nccc1CC', 'CC(=O)Oc1ccccc1C(=O)O', 'CCN(CC)C(=O)c1ccccc1'];
    setChemicalEngine(ChemicalEngine.RDKit);
    const [rdkit] = await getMmpFrags(molecules);
    setChemicalEngine(ChemicalEngine.Crux);
    setCruxUnreadableForTests([molecules[1]]);
    try {
      expect(await cruxReady(), true);
      const fallbacks = cruxFallbackCount('MMPA fragments');
      const [crux] = await getMmpFrags(molecules);
      expect(JSON.stringify(fragmentNames(crux)), JSON.stringify(fragmentNames(rdkit)));
      expect(cruxFallbackCount('MMPA fragments'), fallbacks + 1, 'one fallback: the molecule RDKit fragmented');
    } finally {
      setCruxUnreadableForTests(null);
      setChemicalEngine(null);
    }
  });

  test('R-group decomposition: Crux mode on a Crux worker, as RDKit mode', async () => {
    // the R-group tests' sample on benzene; the SAR table on its MCS as a query core, every matching strategy
    const sar = (await readDataframe('tests/sar-small_test.csv')).col('smiles')!.toList() as string[];
    const runs: [string, string[], string, boolean, string][] = [
      ['sample', R_GROUP_SAMPLE, 'c1ccccc1', false, 'Greedy'],
      // ten R groups: RDKit orders the columns as strings (R1, R10, R2, ...); Crux ordered them by number (the 100K
      // comparison, 2026-10-23)
      ['ten R groups', ['Fc1c(F)c(F)c(-c2c(F)c(F)c(F)c(F)c2F)c(F)c1F', 'Clc1ccc(-c2ccccc2)cc1'], 'c1ccc(-c2ccccc2)cc1',
        false, 'Greedy'],
      ...['Greedy', 'GreedyChunks', 'Exhaustive', 'NoSymmetrization'].map((strategy) =>
        ['sar-small', sar, MCS_R26.sar[0], true, strategy] as [string, string[], string, boolean, string])];
    for (const [name, molecules, core, coreIsQMol, matchingStrategy] of runs) {
      const options = {matchingStrategy, includeTargetMolInResults: true, onlyMatchAtRGroups: false};
      const decompose = async () => {
        const table = DG.DataFrame.fromColumns([column('smiles', molecules)]);
        const {rGroups, highlightCol} = await rGroupsMinilib(table.col('smiles')!, core, coreIsQMol, 0, options);
        return JSON.stringify([rGroups.map((c) => [c.name, c.toList()]), highlightCol?.toList()]);
      };
      const {rdkit, crux} = await inBothEngines(decompose, ['R-group decomposition']);
      expect(crux, rdkit, `${name}, ${matchingStrategy}: columns, R groups and highlights`);
      expect(rdkit.includes('R2'), true, `${name}: decomposed`);
    }
  }, {timeout: 60000});

  test('MCS: Crux mode on a Crux worker, RDKit 2026.03\'s SMARTS', async () => {
    // Chem's RDKit 2024.09 writes some of these MCSs in another atom order (the same query graph): RDKit's own change
    const spgi = (await readDataframe('tests/spgi-100.csv')).col('Structure')!.toList().slice(0, 25);
    const clusters: [keyof typeof MCS_R26, DG.Column][] = [
      ['sar', column('smiles', (await readDataframe('tests/sar-small_test.csv')).col('smiles')!.toList())],
      ['spgi', column('molblocks', spgi)]];
    const rdkit = chemCommonRdKit.getRdKitModule();
    const graph = (smarts: string) => {
      const query = rdkit.get_qmol(smarts);
      const size = `${query.get_num_atoms()} atoms, ${query.get_num_bonds()} bonds`;
      query.delete();
      return size;
    };
    for (const [name, cluster] of clusters) {
      for (const [k, [exactAtoms, exactBonds]] of MCS_COMPARISONS.entries()) {
        const mcs = await inBothEngines(() => getMCS(cluster, exactAtoms, exactBonds), ['MCS']);
        const what = `${name}, atoms ${exactAtoms ? 'exact' : 'any'}, bonds ${exactBonds ? 'exact' : 'any'}`;
        expect(mcs.crux, MCS_R26[name][k], what);
        expect(graph(mcs.rdkit), graph(mcs.crux), `${what}: RDKit mode's MCS`);
      }
    }
  }, {timeout: 60000});

  test('cluster MCS: on the Crux workers, RDKit for a cluster Crux cannot read', async () => {
    // three series of SPGI's MOL blocks (rows 0-9, 10-19, 20-29), the SAR table, a single molecule, an empty cluster
    const spgi = (await readDataframe('tests/spgi-100.csv')).col('Structure')!.toList() as string[];
    const sar = (await readDataframe('tests/sar-small_test.csv')).col('smiles')!.toList() as string[];
    const clustered: [string, string][] = [
      ...spgi.slice(0, 30).map((m, i) => [m, `spgi ${Math.floor(i / 10)}`] as [string, string]),
      ...sar.map((m) => [m, 'sar'] as [string, string]), ['CCOc1ccccc1', 'one'], ['', 'empty'],
      ['CC(=O)Nc1ccc(O)cc1', 'another one']];
    const molCol = column('molecules', clustered.map(([m]) => m));
    molCol.semType = DG.SEMTYPE.MOLECULE;
    const clusterCol = column('cluster', clustered.map(([, c]) => c));
    const mcsOf = async () => (await PackageFunctions.performClusterMCS(molCol, clusterCol)).toList();
    const {rdkit, crux} = await inBothEngines(mcsOf, ['cluster MCS']);
    expect(JSON.stringify(crux), JSON.stringify(rdkit), 'each row\'s cluster MCS');
    const clusterMcs = ['CCNC(C)=O', 'CCN(C)C', 'CCN(C)C1:N:N:C(C):N:1CC', 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1',
      'CCOc1ccccc1', null, 'CC(=O)Nc1ccc(O)cc1'];
    expect(JSON.stringify([...new Set(rdkit)]), JSON.stringify(clusterMcs),
      'the clusters\' MCS (as RDKit.js 2026.03.6 writes them too)');
    setChemicalEngine(ChemicalEngine.Crux);
    // a molecule of the third SPGI series (the RDKit pool's MCS of it) and the one-molecule cluster's (its SMILES
    // written by RDKit's module)
    setCruxUnreadableForTests([spgi[25], 'CC(=O)Nc1ccc(O)cc1']);
    try {
      const fallbacks = cruxFallbackCount('cluster MCS');
      expect(JSON.stringify(await mcsOf()), JSON.stringify(rdkit), 'RDKit\'s MCS for the clusters Crux cannot read');
      expect(cruxFallbackCount('cluster MCS'), fallbacks + 2, 'two fallbacks: the SPGI series, the one molecule');
    } finally {
      setCruxUnreadableForTests(null);
      setChemicalEngine(null);
    }
  }, {timeout: 60000});

  test('SAR matrix: clusters anchored and decomposed on Crux, as on RDKit', async () => {
    // the SAR table and SPGI's molecules, two clusters each: MCS anchors, then R-group decomposition per cluster
    const sar = (await readDataframe('tests/sar-small_test.csv')).col('smiles')!.toList() as string[];
    const spgi = (await readDataframe('tests/spgi-100.csv')).col('Structure')!.toList() as string[];
    const smiles = [...sar, ...spgi.slice(0, 24).map((m) => chemCommonRdKit.getRdKitModule().get_mol(m))
      .map((mol) => {
        const s = mol.get_smiles();
        mol.delete();
        return s;
      })];
    const members = [[0, 1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12, 13, 14, 15, 16],
      [17, 18, 19, 20, 21, 22, 23, 24], [25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]];
    const decompose = async () => JSON.stringify(await decomposeClusters(members, smiles, true));
    const {rdkit, crux} = await inBothEngines(decompose, ['cluster MCS', 'cluster R-groups']);
    expect(crux, rdkit, 'every cluster\'s records and positions');
    expect(rdkit.includes('records'), true, 'some cluster decomposed');
  }, {timeout: 60000});

  test('fragment linking: molzip in Crux mode, as the string assembly in RDKit mode', async () => {
    // the SAR matrix tests' joins: a single point, H, a point the core lacks, stereo by the fragment, a bridge
    const cases: [string, LinkStages, {[position: string]: string}][] = [
      ['c1ccc([*:1])cc1', [{R1: [1]}], {R1: 'C[*:1]'}],
      ['c1ccc([*:1])cc1[*:2]', [{R1: [1], R2: [2]}], {R1: 'C[*:1]', R2: ''}],
      ['c1ccc([*:1])cc1', [{R1: [1], R2: [2]}], {R1: 'C[*:1]', R2: 'Cl[*:2]'}],
      ['c1ccc([*:1])cc1[*:3]', [{R1: [1]}], {R1: 'CC[*:1]'}],
      ...['[*:1][C@@H](C)O', '[*:1][C@H](C)O', '[*:1]/C=C/C', '[*:1][C@](C)(F)Cl'].map((fragment) =>
        ['c%11ccccc%11[*:1]', [{R1: [1]}], {R1: fragment}] as [string, LinkStages, {[position: string]: string}]),
      ['C1CC(C[*:100])(N[*:1])C([*:101])N1[*:2]', [{R1: [1], R2: [2], R100: [100, 101]}],
        {R1: 'CC(=O)[*:1]', R2: 'C[*:2]', R100: 'C(C[*:101])[*:100]'}],
      ['O=C(N[*:1])c1cc([*:2])ccc1', [{R1: [1]}, {R2: [2]}], {R1: '[*:1]C1CC1', R2: 'Br[*:2]'}],
    ];
    const link = () => linkStaged(cases.map(([core]) => core), cases.map(([, plan]) => plan),
      (i, position) => cases[i][2][position] ?? '');
    const {rdkit, crux} = await inBothEngines(link, ['fragment linking']);
    expect(JSON.stringify(crux), JSON.stringify(rdkit), 'each joined molecule');
    expect(rdkit.every((s) => !!s), true, 'every case joined');
  });

  test('molzip: RDKit\'s molzip in both engines', async () => {
    // [molecule, params, the pinned RDKit 2026.03.6's result] (scratch molzip_expected.py)
    const cases: [string, object, string][] = [
      ['C[*:1].N[*:1]', {}, 'CN'],
      ['F[C@H]([*:1])Cl.[*:1]O', {}, 'O[C@@H](F)Cl'],
      ['[*:1]c1ccc([*:2])cc1.[*:1]C(=O)O.[*:2]N', {}, 'Nc1ccc(C(=O)O)cc1'],
      ['C1CC([*:1])CC1[*:2].[*:1]OC[*:2]', {}, 'C1CC2CC1CO2'],
      ['CC([*:1])[*:2].O[*:1]', {}, 'CC(O)[*:2]'],
      ['C/C=C/[*:1].[*:1]Br', {}, 'C/C=C/Br'],
      // a stereo double bond between two attachment points: Crux panicked here (the 100K comparison, 2026-10-23)
      ['C(=C\\[*:2])/[*:1].C[*:1].C[*:2]', {}, 'C/C=C/C'],
      ['C(=N\\[*:2])/[*:1].C[*:1].C[*:2]', {}, 'C/C=N/C'],
      ['C[1*].N[1*]', {Label: 'Isotope'}, 'CN'],
      ['C[Xe].N[Xe]', {Label: 'AtomType', AtomSymbols: ['Xe']}, 'CN'],
      ['CCO', {}, 'CCO'],
      ['C[*:1].N[*:1].O[*:1]', {}, ''], // a label used twice: RDKit throws
      ['[*:1]C(F)(F)(F).[*:1]=O', {}, ''], // over carbon's valence: RDKit throws
      ['', {}, ''],
      ['C1CC(', {}, ''],
    ];
    const molblock = (() => {
      const mol = chemCommonRdKit.getRdKitModule().get_mol('C[*:1].N[*:1]');
      const block = mol.get_molblock();
      mol.delete();
      return block;
    })();
    const molzip = async (molecules: string[], params?: object): Promise<string> =>
      JSON.stringify(await grok.functions.call('Chem:molzip', {molecules, params}));
    for (const engine of [ChemicalEngine.RDKit, ChemicalEngine.Crux]) {
      setChemicalEngine(engine);
      try {
        for (const [molecule, params, expected] of cases)
          expect(await molzip([molecule], params), JSON.stringify([expected]), `${engine}: ${molecule}`);
        expect(await molzip([molblock, 'C[*:1].Cl[*:1]']), JSON.stringify(['CN', 'CCl']),
          `${engine}: a MOL block of two pieces, and a second molecule`);
      } finally {
        setChemicalEngine(null);
      }
    }
  });

  test('a replaced worker: the index a search was using is freed on the other workers', async () => {
    // CruxService.restartWorker drops the column indexes with a segment on the worker it replaces; one a search is
    // using is released (its segments on the other workers freed) when the search ends
    setChemicalEngine(ChemicalEngine.Crux);
    try {
      expect(await cruxReady(), true);
      const service = await getCruxService();
      const smiles = (await readDataframe('smiles.csv')).col('canonical_smiles')!.toList() as string[];
      const index = await service.getIndex(column('molecules', [...smiles, ...smiles, ...smiles]));
      await service.build(index, () => false);
      expect(new Set(index.segments.map((segment) => segment.worker)).size > 1, true, 'indexed on several workers');
      const search = service.search(index, 'c1ccccc1', '', () => {}, () => false);
      service.restartWorkerForTests(index.segments[0].worker);
      await search.catch(() => {});
      expect(index.lost, true, 'the index is dropped');
      expect(index.released.length, 0, 'its segments on the other workers are freed');
    } finally {
      setChemicalEngine(null);
    }
  });

  test('Cluster MCS: a cluster waits for an idle worker, not behind a busy one', async () => {
    // a long call on worker 0 (an FMCS, say): the clusters go to the idle workers, none past its time limit (1.5 s
    // here), and the pool is not restarted
    setChemicalEngine(ChemicalEngine.Crux);
    expect(await cruxReady(), true);
    const service = await getCruxService();
    const sar = (await readDataframe('tests/sar-small_test.csv')).col('smiles')!.toList() as string[];
    service.clusterLimitForTests = 1500;
    try {
      const generation = service.generation;
      const busy = service.busyWorkerForTests(0, 4000);
      await new Promise((resolve) => setTimeout(resolve, 200)); // the busy call posted
      const {results, timedOut} = await service.clusterMCS([sar.slice(0, 6), sar.slice(6, 12), sar.slice(12)],
        true, true);
      expect(JSON.stringify(timedOut), '[]', 'clusters past their time limit');
      expect(results.every((mcs) => !!mcs), true, 'every cluster has an MCS');
      await busy;
      expect(service.generation, generation, 'the pool was restarted');
    } finally {
      service.clusterLimitForTests = null;
      setChemicalEngine(null);
    }
  }, {timeout: 60000});

  test('a replaced worker\'s calls fail without restarting the pool', async () => {
    // the jobs of a worker CruxService.restartWorker replaces fail (their callers fall back to RDKit); the other
    // workers and their indexes stay
    setChemicalEngine(ChemicalEngine.Crux);
    try {
      expect(await cruxReady(), true);
      const service = await getCruxService();
      const generation = service.generation;
      const busy = service.busyWorkerForTests(1, 3000);
      await new Promise((resolve) => setTimeout(resolve, 200)); // the busy call posted
      service.restartWorkerForTests(1);
      const failure = await busy.then(() => 'no failure', (e: any) => String(e?.message ?? e));
      expect(failure.includes('was replaced'), true, failure);
      expect(service.generation, generation, 'the pool was restarted');
      // the new worker answers (a new client: the old one's terminate reason is not carried over)
      expect(typeof await service.busyWorkerForTests(1, 10), 'number', 'the replacement answers');
    } finally {
      setChemicalEngine(null);
    }
  }, {timeout: 60000});

  test('timing: RDKit and Crux', async () => {
    // a measurement, no check: with --benchmark only (it returns at once otherwise). Median of 3, each engine: MMPA
    // fragmentation of the MMP demo (20,267 molecules), R-group decomposition of the R-group tests' benchmark table
    // (smiles.csv's 1,000 on benzene), FMCS of the MCS tests' cluster (SPGI-100) by element and by any atom
    if (!DG.Test.isInBenchmark)
      return;
    const mmp = (await readDataframe('demo_files/mmp_demo.csv')).col('smiles')!.toList() as string[];
    const smiles = (await readDataframe('smiles.csv')).col('canonical_smiles')!.toList() as string[];
    const rgroups = DG.DataFrame.fromColumns([column('smiles', smiles)]).col('smiles')!;
    const mcs = column('molblocks', (await readDataframe('tests/spgi-100.csv')).col('Structure')!.toList() as string[]);
    const options = {matchingStrategy: 'Greedy', includeTargetMolInResults: true, onlyMatchAtRGroups: false};
    const jobs: [string, () => Promise<unknown>][] = [
      [`MMPA fragments, ${mmp.length} molecules`, () => getMmpFrags(mmp)],
      [`R-groups, ${rgroups.length} molecules`, () => rGroupsMinilib(rgroups, 'c1ccccc1', false, 0, options)],
      [`FMCS by element, ${mcs.length} molecules`, () => getMCS(mcs, true, true)],
      [`FMCS any atom, ${mcs.length} molecules`, () => getMCS(mcs, false, true)],
    ];
    const report: string[] = [];
    for (const [name, job] of jobs) {
      const medians: string[] = [];
      for (const engine of [ChemicalEngine.RDKit, ChemicalEngine.Crux]) {
        setChemicalEngine(engine);
        try {
          await cruxReady();
          await job(); // warm: the pools up, the molecules read once
          const times: number[] = [];
          for (let k = 0; k < 3; k++) {
            const started = performance.now();
            await job();
            times.push(performance.now() - started);
          }
          medians.push(`${engine} ${times.sort((a, b) => a - b)[1].toFixed(0)} ms`);
        } finally {
          setChemicalEngine(null);
        }
      }
      report.push(`${name}: ${medians.join(', ')}`);
    }
    console.log(`crux analyses | timing (median of 3): ${report.join('; ')}`);
  }, {timeout: 600000, benchmark: true, benchmarkTimeout: 3600000});
});
