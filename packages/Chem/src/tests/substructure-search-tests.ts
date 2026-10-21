import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {readDataframe, _testSearchSubstructure, _testSearchSubstructureAllParameters,
  _testSearchSubstructureSARSmall,
  checkBitSetIndices} from './utils';
import {after, before, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {_package} from '../package-test';
import {PackageFunctions} from '../package';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {FILTER_TYPES, chemSubstructureSearchLibrary} from '../chem-searches';
import { SubstructureSearchType } from '../constants';
import { Fingerprint } from '../utils/chem-common';
import {SubstructureSearchEngine, setSubstructureSearchEngine} from '../crux/crux-searches';
import {getMolSafe, getQueryMolSafe} from '../utils/mol-creation_rdkit';
import {queryNeedsSanitizedMolecules} from '../rdkit-service/query-features';

export const testSubstructure = 'C1CC1';
//csv with C1CC1 as substructure in pos 0 and 2
export const testCsv = `smiles
CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4
COc1ccc2cc(ccc2c1)C(C)C(=O)Oc3ccc(C)cc3OC
Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5
CC(C(=O)NCCS)c1cccc(c1)C(=O)c2ccccc2
COc1ccc2c(c1)c(CC(=O)N3CCCC3C(=O)Oc4ccc(C)cc4OC)c(C)n2C(=O)c5ccc(Cl)cc5`;

const included_in_molblock = `
RDKit          2D

 11 12  0  0  0  0  0  0  0  0999 V2000
    3.0000   -2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000   -2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000   -2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -3.8971    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -3.8971    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  2  0
  3  8  1  0
  8  9  2  0
  9 10  1  0
 10 11  2  0
 11  2  1  0
  8  7  1  0
M  END
`

category('substructure search', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('searchSubstructure.sar-small_awaitAll', async () => {
    await _testSearchSubstructureAllParameters(
      _testSearchSubstructureSARSmall);
  });

  test('searchSubstructureColWithoutDf_awaitAll', async () => {
    const col = DG.Column.fromStrings('smiles', [
      'CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4',
      'Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5',
    ]);
    await grok.chem.searchSubstructure(col, 'C1CC1');
  });

  //Number of molecules is smaller than a number of threads
  test('searchSubstructure.5_rows_awaitAll', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
  });

  test('searchSubstructureWithMalformedMolString_awaitAll', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    df.columns.byName('smiles').set(4, 'qq');
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
  });

  test('searchSubstructureWithNull_awaitAll', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    df.columns.byName('smiles').set(4, null);
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
  });

  test('searchSubstructurePerformance_awaitAll', async () => {
    const df = grok.data.demo.molecules(50000);
    const startTime = performance.now();
    await grok.chem.searchSubstructure(df.columns.byName('smiles'), 'c1ccccc1')!;
    const midTime1 = performance.now();
    for (let i = 0; i < 20; i++)
      await grok.chem.searchSubstructure(df.columns.byName('smiles'), 'c1ccccc1')!;
    const midTime2 = performance.now();
    await grok.chem.searchSubstructure(df.columns.byName('smiles'), 'C1CC1')!;
    const endTime = performance.now();

    console.log(`first Call to searchSubstructure took ${midTime1 - startTime} milliseconds`);
    console.log(`loop Call to searchSubstructure took ${midTime2 - midTime1} milliseconds`);
    console.log(`last Call to searchSubstructure took ${endTime - midTime2} milliseconds`);
  }, {timeout: 180000});

  test('searchSubstructureAfterChanges_awaitAll', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
    df.columns.byName('smiles').set(0, 'COc1ccc2cc(ccc2c1)C(C)C(=O)Oc3ccc(C)cc3OC');
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [2]);
  });

  test('searchSubstructureWith2Columns_awaitAll', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    const col = DG.Column.fromStrings('smiles1', [
      'COc1ccc2cc(ccc2c1)C(C)C(=O)Oc3ccc(C)cc3OC',
      'CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4',
      'COc1ccc2c(c1)c(CC(=O)N3CCCC3C(=O)Oc4ccc(C)cc4OC)c(C)n2C(=O)c5ccc(Cl)cc5',
      'COc1ccc2c(c1)c(CC(=O)N3CCCC3C(=O)Oc4ccc(C)cc4OC)c(C)n2C(=O)c5ccc(Cl)cc5',
      'Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5',
    ]);
    df.columns.add(col);
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
    await _testSearchSubstructure(df, 'smiles1', testSubstructure, [1, 4]);
  });

  test('searchSubstructure2Dataframes_awaitAll', async () => {
    const df1 = grok.data.demo.molecules(50);
    const df2 = grok.data.demo.molecules(50);

    await _testSearchSubstructure(df1, 'smiles', 'c1ccncc1', [6, 26, 46]);
    await _testSearchSubstructure(df2, 'smiles', 'c1ccccc1',
      [0, 4, 5, 7, 8, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 24,
        25, 27, 28, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 47, 48]);
  });

  test('searchSubstructureExplicitHydrogen_awaitAll', async () => {
    const df = await readDataframe('tests/explicit_h_test.csv');
    await _testSearchSubstructure(df, 'smiles', `
  MJ201900                      

  5  5  0  0  0  0  0  0  0  0999 V2000
    -0.1562   -0.7750    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    -0.8707   -0.3625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -0.8707    0.4625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5582    0.4625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5582   -0.3625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  5  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
M  END`, [0, 1, 2, 3]);

    await _testSearchSubstructure(df, 'smiles', `
  MJ201900                      

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.5357    0.1401    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    -0.1787    0.5526    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -0.1787    1.3776    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2501    1.3776    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2501    0.5526    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5357   -0.6848    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  5  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  1  6  1  0  0  0  0
M  END
  `, [2, 3]);
  });


  test('searchSubstructureAromaticBond_awaitAll', async () => {
    const df = await readDataframe('tests/aromatic_bond_test.csv');
    await _testSearchSubstructure(df, 'smiles', `
  MJ201900                      

  6  6  0  0  0  0  0  0  0  0999 V2000
    -1.3616    0.4008    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -2.0761   -0.0116    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -2.0761   -0.8367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -1.3616   -1.2491    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -0.6471   -0.8367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -0.6471   -0.0116    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  6  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  4  0  0  0  0
M  END
`, [2, 3]);
  });

  test('search_benzene_awaitAll', async () => {
    const df = DG.Test.isInBenchmark ? await readDataframe('tests/smi10K.csv') :
      grok.data.demo.molecules(500);
    await performanceTestWithConsoleLog(df.col('smiles')!, 'c1ccccc1');
  }, {benchmark: true});

  test('search_2_benzene_rings_awaitAll', async () => {
    const df = DG.Test.isInBenchmark ? await readDataframe('tests/smi10K.csv') :
      grok.data.demo.molecules(500);
    await performanceTestWithConsoleLog(df.col('smiles')!, 'c1ccc2ccccc2c1');
  }, {benchmark: true});

  test('search_complex_structure_awaitAll', async () => {
    const df = DG.Test.isInBenchmark ? await readDataframe('tests/smi10K.csv') :
      grok.data.demo.molecules(500);
    await performanceTestWithConsoleLog(df.col('smiles')!, 'c1nn2cnnc2s1');
  }, {benchmark: true});

  test('search_non_ring_structure_awaitAll', async () => {
    const df = DG.Test.isInBenchmark ? await readDataframe('tests/smi10K.csv') :
      grok.data.demo.molecules(500);
    await performanceTestWithConsoleLog(df.col('smiles')!, 'CNC(C)=O');
  }, {benchmark: true});

  test('rdkitReinitKeepsWorkerPool', async () => {
    // a second init (another category's before) used to start a new pool and leave the old one's workers running, each
    // holding an RDKit WebAssembly instance, until the page had room for no more ("Out of memory")
    const service = await chemCommonRdKit.getRdKitService();
    await chemCommonRdKit.initRdKitModuleLocal();
    expect(await chemCommonRdKit.getRdKitService() === service, true, 'the RDKit worker pool after a second init');
  });

  test('rdkitWorkerPoolIsBounded', async () => {
    // each worker holds an RDKit WebAssembly instance, and a page has room for only so many WebAssembly memories
    const service = await chemCommonRdKit.getRdKitService();
    expect(service.workerCount, Math.max(1, Math.min(navigator.hardwareConcurrency - 2, 16)), 'RDKit workers');
  });

});

category('substructure search: search types', () => { 

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('included_in_search_awaitAll', async () => {
    await testSearchType(included_in_molblock, SubstructureSearchType.INCLUDED_IN, [0, 1, 2]);
  });

  test('not_included_in_search_awaitAll', async () => {
    await testSearchType(included_in_molblock, SubstructureSearchType.NOT_INCLUDED_IN,
      [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]);
  });

  test('contains_search_awaitAll', async () => {
    await testSearchType('c1ccccc1', SubstructureSearchType.CONTAINS,
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 16, 17, 18, 19]);
  });

  test('not_contains_search_awaitAll', async () => {
    await testSearchType('c1ccccc1', SubstructureSearchType.NOT_CONTAINS, [13, 14]);
  });

  test('exact_search_awaitAll', async () => {
    await testSearchType('FC(F)(F)c1ccc(OC2CCNCC2)cc1', SubstructureSearchType.EXACT_MATCH, [4]);
  });

  test('is_similar_search_awaitAll', async () => {
    await testSearchType('c1ccccc1', SubstructureSearchType.IS_SIMILAR, [0, 1, 2, 3, 5, 15, 16], Fingerprint.Morgan, 0.05);
  });

});

category('substructure search: stereo agnostic', () => {

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('tetrahedral_stereo', async () => {
    // Query R-fluoroethanol -> should match R, S, and flat (indices 0, 1, 2)
    await testStereoAgnosticSearch('C[C@@H](O)F', [0, 1, 2]);
  });

  test('cis_trans_stereo', async () => {
    // Query trans-chlorofluoroethene -> should match trans, cis, and flat (indices 3, 4, 5)
    await testStereoAgnosticSearch('F/C=C/Cl', [3, 4, 5]);
  });

  test('multiple_stereocenters', async () => {
    // Query trans-cyclohexane-1,4-diol -> should match trans, cis, and flat (indices 6, 7, 8)
    await testStereoAgnosticSearch('O[C@@H]1CC[C@H](O)CC1', [6, 7, 8]);
  });

  test('complex_stereo_glucose', async () => {
    // Query D-glucose -> should match D-glucose, L-glucose, and flat (indices 9, 10, 11)
    await testStereoAgnosticSearch('C([C@H]1[C@@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O', [9, 10, 11]);
  });

  test('flat_query', async () => {
    // Query without stereo -> should still match all stereoisomers (indices 0, 1, 2)
    await testStereoAgnosticSearch('CC(O)F', [0, 1, 2]);
  });
});

category('substructure search: filter operators', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('descriptor', async () => {
    const funcs = DG.Func.find({package: 'Chem', meta: {role: 'filterOperators'}});
    expect(funcs.length, 1);
    expect(funcs[0].options['semType'], DG.SEMTYPE.MOLECULE);
    const set = await funcs[0].apply();
    expect(set.semType, DG.SEMTYPE.MOLECULE);
    expect(set.exclusive, true);
    expectArray(set.operators.map((o: any) => o.id),
      [SubstructureSearchType.CONTAINS, SubstructureSearchType.INCLUDED_IN, SubstructureSearchType.EXACT_MATCH,
        SubstructureSearchType.STEREO_AGNOSTIC, SubstructureSearchType.IS_SIMILAR]);
  });

  test('contains_matches_searchSubstructure', async () => {
    const col = grok.data.demo.molecules(200).col('smiles')!;
    col.semType = DG.SEMTYPE.MOLECULE;
    const contains = PackageFunctions.moleculeFilterOperators().operators
      .find((o) => o.id === SubstructureSearchType.CONTAINS)!;
    const mask = await contains.bitset(col, {property: col.name, operator: contains.id, value: 'c1ccccc1'},
      new AbortController().signal);
    const expected: DG.BitSet = (await PackageFunctions.searchSubstructure(col, 'c1ccccc1', '')).get(0);
    expect(mask.length, col.length);
    const actual = DG.BitSet.fromBytes(mask.bits.buffer as ArrayBuffer, mask.length);
    expect(actual.trueCount > 0 && actual.trueCount < col.length, true);
    expect(actual.toBinaryString(), expected.toBinaryString());
  });
});

const BMS_304 = '[#6;r16,r17,r18]~[#6]1~[#6]~[#6]~[#6](~[#6])~[#7]1';
/** Molecules the RDKit search answered ring, hybridization and radical queries for otherwise than RDKit does when it
 * read them unsanitized (their rings from a depth-first search, no hybridization, no radicals). */
const RING_CASES = [
  'CCOc1ccc2ccccc2c1', 'NCCc1cccc2ccccc12', // naphthalenes: a 10-ring and a 6-ring to the search, two 6-rings to RDKit
  'OC1CC2CCC1C2', // norbornane: a 6-ring and a 5-ring to the search, two 5-rings to RDKit
  'OCC12C3C4C1C5C2C3C45', // cubane: RDKit's six 4-rings (symmetrized SSSR), each atom in three
  'CCCCCCOC(=O)c1ccncc1', // sp3 carbons
  '[CH2]CCOC(=O)c1ccccc1', // a radical
  // BMS 304 (a structural alert): a true hit (ChEMBL-10k's row 8751) and a molecule only the search's rings hit (5418)
  'C=Cc1c(C)c2cc3nc(cc4[nH]c(cc5nc(cc1[nH]2)C(C)=C5C(C)OCCCCCC)c(C)c4CCC(=O)O)C(CCC(=O)O)=C3C',
  'c1ccc(-c2c3nc4ccccc4c-3[nH]c3ccccc23)cc1',
];
/** A query of a radical: the MOL block of an ethyl radical. */
function radicalQuery(): string {
  const mol = chemCommonRdKit.getRdKitModule().get_mol('[CH2]C');
  try {
    return mol.get_molblock();
  } finally {
    mol.delete();
  }
}
/** The rows of RING_CASES each query hits in RDKit (the pinned RDKit 2026.03 and Chem's RDKit.js agree; molecules
 * sanitized); the search answered [0,1], [0,1,2,3,4,5,7], [0,1,2,3,6,7], [3,6,7], [], [6,7] and [] for the first seven. */
const ringCasesHits = (): [string, number[]][] => [
  ['[c;r10]', []], ['[#6;r6]', [0, 1, 4, 5, 7]], ['[#6;R2]', [0, 1, 2, 6, 7]], ['[#6;R3]', [3]],
  ['[C^3]', [0, 1, 2, 3, 4, 5, 6]], [BMS_304, [6]], [radicalQuery(), [5]],
  // the search's answer was RDKit's for these: ring membership, ring bonds, no ring property at all
  ['[c;R2]', [0, 1, 6, 7]], ['[#6;R]', [0, 1, 2, 3, 4, 5, 6, 7]], ['[#6;x3]', [0, 1, 2, 3, 6, 7]], ['C1CCCC1', [2]],
  ['c1ccccc1', [0, 1, 5, 7]],
];

const hitRows = (hits: BitArray): number[] => Array.from({length: hits.length}, (_, i) => i).filter((i) => hits.getBit(i));
const queryLabel = (query: string): string => query.includes('M  END') ? 'a radical (MOL block)' : query;

category('substructure search: ring queries', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    setSubstructureSearchEngine(SubstructureSearchEngine.RDKit);
  });

  after(async () => setSubstructureSearchEngine(null));

  const rdkitSearch = (col: DG.Column, query: string) =>
    chemSubstructureSearchLibrary(col, query, '', FILTER_TYPES.substructure, false, true);

  test('queryNeedsSanitizedMolecules', async () => {
    // decided from the query molecule (its RDKit JSON, recursive SMARTS included): ring counts and sizes, hybridization
    // and radicals need the molecules sanitized; ring membership (R, R0, r, x, @) does not
    const cases: [string, boolean][] = [['[c;r10]', true], ['[#6;R2]', true], ['[#6;R{2-}]', true], ['[#6;!R1]', true],
      ['[C^3]', true], ['[$([#6;r6])]N', true], [BMS_304, true], [radicalQuery(), true], ['[#6;R]', false],
      ['[#6;R0]', false], ['[#6;!R]', false], ['[#6;r]', false], ['[#6;x2]', false], ['C@C', false], ['c1ccccc1', false],
      ['[CH2]', false], ['[N+](=O)[O-]', false], ['[$([#7;x2])]C', false]];
    for (const [query, needs] of cases) {
      const queryMol = getQueryMolSafe(query, '', chemCommonRdKit.getRdKitModule());
      try {
        expect(queryNeedsSanitizedMolecules(queryMol!), needs, queryLabel(query));
      } finally {
        queryMol?.delete();
      }
    }
  });

  test('editedRowsHitAsInRdkit', async () => {
    // an edited row's fingerprint comes from one RDKit worker (the last one, for one row) and the search reads the row on
    // another (the first one, for a small table), which does not hold the row's molecule: it read it unsanitized
    const col = DG.Column.fromStrings('smiles', RING_CASES.map(() => 'CC'));
    col.semType = DG.SEMTYPE.MOLECULE;
    DG.DataFrame.fromColumns([col]);
    await rdkitSearch(col, '[#6]');
    for (let i = 0; i < RING_CASES.length; i++) {
      col.set(i, RING_CASES[i]);
      await rdkitSearch(col, '[#6]');
    }
    for (const [query, rows] of ringCasesHits())
      expect(JSON.stringify(hitRows(await rdkitSearch(col, query))), JSON.stringify(rows), queryLabel(query));
  });

  test('workerHitsAsInRdkit', async () => {
    // the RDKit worker's search, on molecules its cache does not hold (canonical SMILES, as the search passes them)
    const rdkit = chemCommonRdKit.getRdKitModule();
    const molecules = ['COc1ccc2cc(CCC(C)=O)ccc2c1', 'Cc1ccc2cc(Br)ccc2c1', 'C12C3C4C1C5C2C3C45', 'CCCOC(=O)c1ccccc1']
      .map((smiles) => {
        const mol = getMolSafe(smiles, {}, rdkit).mol!;
        try {
          return mol.get_smiles();
        } finally {
          mol.delete();
        }
      });
    // RDKit's rows; the worker answered [0,1], [0,1,2] and [] for the first three
    const expected: [string, number[]][] = [['[c;r10]', []], ['[#6;R2]', [0, 1]], ['[C^3]', [0, 1, 2, 3]],
      ['[#6;R3]', [2]], ['c1ccccc1', [0, 1, 3]]];
    const worker = (await chemCommonRdKit.getRdKitService()).parallelWorkers[0];
    for (const [query, rows] of expected) {
      const hits = BitArray.fromUint32Array(molecules.length,
        await worker.searchSubstructure(query, '', molecules, SubstructureSearchType.CONTAINS));
      expect(JSON.stringify(hitRows(hits)), JSON.stringify(rows), query);
    }
  });

  test('ringQuerySpeed', async () => {
    // a plain query and a ring query on a fresh table (fingerprints, then the search), then again (fingerprints saved)
    const rows = DG.Test.isInBenchmark ? 100000 : 5000;
    const times: string[] = [];
    for (const query of ['c1ccccc1', '[#6;R2]', BMS_304]) {
      const col = grok.data.demo.molecules(rows).col('smiles')!;
      col.semType = DG.SEMTYPE.MOLECULE;
      const t0 = performance.now();
      const hits = await rdkitSearch(col, query);
      const t1 = performance.now();
      await rdkitSearch(col, query);
      const t2 = performance.now();
      times.push(`${query}: ${hits.trueCount()} hits, ${Math.round(t1 - t0)} + ${Math.round(t2 - t1)} ms`);
    }
    return `${rows} rows; ${times.join('; ')}`;
  }, {benchmark: true, timeout: 120000, benchmarkTimeout: 600000});
});


async function performanceTestWithConsoleLog(molCol: DG.Column, query: string) {
  const startTime = performance.now();
  await chemSubstructureSearchLibrary(molCol, query, '');
  const midTime1 = performance.now();
  await chemSubstructureSearchLibrary(molCol, query, '');
  const midTime2 = performance.now();

  console.log(`first Call to substructure search took ${midTime1 - startTime} milliseconds`);
  console.log(`second Call to substructure search took ${midTime2 - midTime1} milliseconds`);
}

async function testSearchType(query: string, searchType: SubstructureSearchType, trueIndices: number[],
  fp?: Fingerprint, similarity?: number) {
  const df = await readDataframe('tests/test_search_types.csv');
  const res = await chemSubstructureSearchLibrary(df.col('canonical_smiles')!, query, '', FILTER_TYPES.substructure, false, true, searchType, similarity, fp);
  const bitset = DG.BitSet.fromBytes(res.buffer.buffer, df.rowCount);
  checkBitSetIndices(bitset, trueIndices);
}

async function testStereoAgnosticSearch(query: string, trueIndices: number[]) {
  const df = await readDataframe('stereo-agnostic-test.csv');
  const res = await chemSubstructureSearchLibrary(df.col('smiles')!, query, '', FILTER_TYPES.substructure, false, true,
    SubstructureSearchType.STEREO_AGNOSTIC);
  const bitset = DG.BitSet.fromBytes(res.buffer.buffer, df.rowCount);
  checkBitSetIndices(bitset, trueIndices);
}
