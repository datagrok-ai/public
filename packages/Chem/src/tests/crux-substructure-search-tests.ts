import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {after, awaitCheck, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {first} from 'rxjs/operators';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {FILTER_TYPES, chemSubstructureSearchLibrary} from '../chem-searches';
import {SubstructureSearchType, getSearchProgressEventName, getSearchQueryAndType,
  getTerminateEventName} from '../constants';
import {Fingerprint} from '../utils/chem-common';
import {SubstructureSearchEngine, cruxReadsQuery, getCruxService, setSubstructureSearchEngine}
  from '../crux/crux-searches';
import {getMolSafe, getQueryMolSafe} from '../utils/mol-creation_rdkit';
import {readDataframe} from './utils';

const QUERIES = ['c1ccccc1', 'c1ccncc1', 'C1CC1', 'C(=O)N', 'c1ccc2[nH]ccc2c1', 'C[C@H](N)C(=O)O', 'C/C=C/C',
  '[CX3](=O)[OX2H1]', '[NX3;H2,H1;!$(NC=O)]', '[F,Cl,Br,I]', '[#7;R]', '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-,=[#8]'];
// Queries Chem sent to RDKit's worker pool while Crux took the query as SMARTS text (getCruxSmarts' refusals): explicit
// [#1] and !#1 (PAINS writes them), R<n>, r<n>, v, h, x, isotopes, radicals, exclusions that may be an H (!C, !a, *!O),
// some of Chem's structural alerts among them; Crux now matches the query molecule Chem's getQueryMolSafe builds.
// On Chem's datasets here each gives the hits of RDKit 2026.03 (RDKit.js, molecules read as getMolSafe reads them) and
// of Chem's RDKit.js, row for row (checked with both before these tests were written).
const FORMERLY_RDKIT_QUERIES = [
  '[!#6;!#1;!#11;!#19]O(S(=O)(=O)(-[C,c]))', '[#7+](:[!#1]:[!#1]:[!#1])-[!#1]=[#8]',
  '[!#6;!#7;!#8;!#16;!#1;!#3;!#9;!#11;!#12;!#15;!#17;!#19;!#20;!#30;!#35]', '[H]N([H])c1ccccc1', '[#1]', '[#6]-[#1]',
  '[#8]-[#1]', '[!#1]:1:[!#1]:[!#1]:[!#1]:[!#1]:[!#1]:1',
  '[#6&R2]', '[O;R1][C;R1][C;R1][O;R1]', '[c;R2]', '[C;R0](=O)N', '[#6;R3]',
  '[C;r6]', '[#6;r5]', '[N;!r3]-[N]', '[#7;r5,r6]', '[c;r10]', '[#6;r16,r17,r18]~[#6]1~[#6]~[#6]~[#6](~[#6])~[#7]1',
  '[C;v4]', '[N;v3]', '[c,C]-[P;v3]', '[#16v3,#16v5]', '[Nv3X3][Nv3X3!H0]', '[S;v6]',
  '[C;!$(C[N])](=O)!@[C;h1,h2;H1,H2][F,Cl,Br,I]', '[C;h3]', '[N;h2]', '[#6;x2]', '[#6;x3]', '[N;x0]',
  '[13CH3]O', '[2H]', '[13#6]', '[2#1,3#1,11C,11c,14C,14c,125I,32P,33P,35S]', '[0#6]',
  'O=C(-[!N])O', '[Si]~[!#6]', '[!a][SX2;H1]', '[Cl,Br,I]C=[!O!R]', 'C[CR0]=[NR0][*!O]', '[!#6][#6](=O)[!#6]',
  '[CH2]', '[O]', '[#6;D4]', '[#6;D1]', '[N;D1]', '[C;X4;H3]',
];
// molecules with what those queries ask about: isotopes, radicals, explicit H, P / S / Si valences, fused, spiro,
// bridged and large rings
const SPECIAL_CASES = ['[13CH3]OC', 'C[2H]', '[2H]C([2H])([2H])O', '[14C]c1ccccc1', 'O[13C](=O)C', '[CH2]C(=O)O',
  '[O]CC', 'C[N+](=O)[O-]', '[N+](=O)[O-]', '[CH3]', 'CS(C)(=O)=O', 'CS(=O)C', 'CP(C)C', 'CP(=O)(C)C', 'C[Si](C)(C)O',
  'c1ccc2ccccc2c1', 'C1CCC2(CC1)CCCC2', 'C1CC2CCC1C2', 'C1CCCCCCCCCCCCCCC1C1CCCCN1', 'OC(=O)C(Cl)C(=O)C', '[H]OC',
  '[H]N([H])c1ccccc1', 'CC(=O)Cl', 'C=CBr', 'NN', 'CSC(=S)SC', '[125I]c1ccccc1', 'C[32P](=O)(O)O', 'CC(C)(C)[2H]',
  'c1cc[nH]c1', '[Na+].[Cl-]', 'CCC(=N)C', 'C#N', ''];
// drawn queries of those kinds, as MOL blocks: a radical (N, M  RAD), an H count (MDL's hhh), an explicit H atom
const RADICAL_MOLBLOCK = ['', '     RDKit          2D', '', '  3  2  0  0  0  0  0  0  0  0999 V2000',
  '    0.0000    0.0000    0.0000 N   0  0  0  0  0  3  0  0  0  0  0  0',
  '    1.2990    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -1.2990    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
  '  1  2  2  0', '  1  3  1  0', 'M  CHG  2   1   1   3  -1', 'M  RAD  1   1   2', 'M  END', ''].join('\n');
const H_COUNT_MOLBLOCK = ['', '  chem-census query', '', '  7  7  0  0  0  0  0  0  0  0999 V2000',
  '    0.0000    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.0000   -1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.0000    3.0000    0.0000 N   0  0  0  3  0  0  0  0  0  0  0  0',
  '  1  2  2  0', '  2  3  1  0', '  3  4  2  0', '  4  5  1  0', '  5  6  2  0', '  6  1  1  0', '  1  7  1  0',
  'M  END'].join('\n');
const EXPLICIT_H_MOLBLOCK = ['', '  chem-census query', '', '  3  2  0  0  0  0  0  0  0  0999 V2000',
  '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    1.3000    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
  '    2.6000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0',
  '  1  2  1  0', '  2  3  1  0', 'M  END'].join('\n');
// RDKit's SMARTS writer writes MDL's single or double with ring topology `-,=&@`, which a SMARTS parser, crux's too,
// reads as `-,(=&@)` (crux-sketch spike query-roundtrip, L6: a custom bond query `-,=;@` passed 969 rows of mol1K for
// 417); as a molecule, the query is the one RDKit matches
const AMBIGUOUS_BOND_MOLBLOCK = ['', '  query-roundtrip', '', '  2  1  0  0  0  0  0  0  0  0999 V2000',
  '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '  1  2  5  0  0  1', 'M  END', ''].join('\n');
const FORMERLY_RDKIT_MOLBLOCKS = [RADICAL_MOLBLOCK, H_COUNT_MOLBLOCK, EXPLICIT_H_MOLBLOCK, AMBIGUOUS_BOND_MOLBLOCK];
// molecules crux read differently from RDKit before crux-core's rdkit-parity-gaps: nitro / N-oxide / azide /
// halogen oxide in their neutral hypervalent form, a radium salt, leading spaces, a cyclic iodine
const SANITIZATION_CASES = ['CC(C)(C)c1ccc(cc1)S(=O)(=O)Nc2ccc(Cl)cc2C(=O)c3ccn(=O)cc3',
  'O=C1NCCN1c2ncc(s2)N(=O)=O', 'CCOc1ccc(cc1)c2cc(C)[n+](CCO)c(c2)c3ccc(OCC)cc3.[O-]Cl(=O)(=O)=O',
  'CNS(=O)(=O)CCCCN=N#N', '[Cl-].[Cl-].[Ra+2]', '  CC[C@@H]1OCC[C@H]1C(=O)N1CCOC(C)(CNC(=O)C2=CC(F)=CN2)C1',
  'I1c2ccccc2c3ccccc13'];
// molecules only RDKit's read without Kekulize accepts, the one Chem retries with
const LENIENT_CASES = ['O=C1N(CCc3c1c2ccccc2n3C)Cc4ncnc4C', 'c2ccc1nncc1c2',
  'Clp1(Cl)np(Cl)(Cl)np2(NNP(=O)(NN2)Oc3ccccc3)n1', 'N1c2ccccc2p3(c4ccccc14)c5ccccc5nc6ccccc36'];

async function search(engine: SubstructureSearchEngine, col: DG.Column, query: string,
  searchType = SubstructureSearchType.CONTAINS, includeMask: BitArray | null = null): Promise<BitArray> {
  setSubstructureSearchEngine(engine);
  return chemSubstructureSearchLibrary(col, query, '', FILTER_TYPES.substructure, false, true, searchType, 0.8,
    Fingerprint.Morgan, includeMask);
}

async function expectSameAsRdkit(col: DG.Column, queries: string[], searchType = SubstructureSearchType.CONTAINS,
  includeMask: BitArray | null = null, viaCrux = true): Promise<void> {
  for (const query of queries) {
    if (viaCrux)
      expect(await cruxReadsQuery(query, ''), true, `crux does not search for ${query}`);
    const rdkit = await search(SubstructureSearchEngine.RDKit, col, query, searchType, includeMask);
    const crux = await search(SubstructureSearchEngine.Crux, col, query, searchType, includeMask);
    expect(crux.equals(rdkit), true, `crux differs from RDKit: ${searchType} ${query} (crux ${crux}, RDKit ${rdkit})`);
  }
}

/** The RDKit worker pool's asks (its init counts every one). */
const rdkitPoolAsks = () => chemCommonRdKit._rdKitService.timesInitialized;

/**
 * The rows RDKit matches, Chem's main-thread RDKit.js reading the molecules (getMolSafe) and the query
 * (getQueryMolSafe) as Chem does: sanitized molecules, their rings the SSSR (RDKit's definition of R<n> and r<n>). The
 * RDKit worker pool is not asked (its search reads the molecules its cache lacks unsanitized, and sanitized for a query
 * that needs it: `substructure search: ring queries`).
 */
function rdkitHits(molecules: (string | null)[], query: string, searchType = SubstructureSearchType.CONTAINS):
  BitArray {
  const rdkit = chemCommonRdKit.getRdKitModule();
  const queryMol = getQueryMolSafe(query, '', rdkit);
  const hits = new BitArray(molecules.length);
  try {
    for (let i = 0; i < molecules.length; i++) {
      const mol = molecules[i] ? getMolSafe(molecules[i]!, {}, rdkit).mol : null;
      if (!mol)
        continue;
      try {
        const matches = !!queryMol && mol.get_substruct_match(queryMol) !== '{}';
        hits.setFast(i, searchType === SubstructureSearchType.CONTAINS ? matches : !matches);
      } finally {
        mol.delete();
      }
    }
  } finally {
    queryMol?.delete();
  }
  return hits;
}

/** Crux searches each query and gives RDKit's hits ({@link rdkitHits}), asking nothing of the RDKit worker pool. */
async function expectRdkitHitsThroughCrux(col: DG.Column, queries: string[],
  searchType = SubstructureSearchType.CONTAINS): Promise<void> {
  const molecules = col.toList();
  const asks = rdkitPoolAsks();
  for (const query of queries) {
    const label = query.includes('M  END') ? `a MOL block (${query.split('\n')[1].trim()})` : query;
    expect(await cruxReadsQuery(query, ''), true, `Crux does not search for ${label}`);
    const crux = await search(SubstructureSearchEngine.Crux, col, query, searchType);
    const rdkit = rdkitHits(molecules, query, searchType);
    expect(crux.equals(rdkit), true, `Crux differs from RDKit: ${col.name} ${searchType} ${label} ` +
      `(Crux ${crux.trueCount()} rows, RDKit ${rdkit.trueCount()}; the first differing row ` +
      `${crux.clone().xor(rdkit).findNext(-1, true)})`);
  }
  expect(rdkitPoolAsks(), asks, `a Crux search on ${col.name} asked for the RDKit workers`);
}

function moleculeColumn(smiles: string[]): DG.Column {
  const col = DG.Column.fromStrings('smiles', smiles);
  col.semType = DG.SEMTYPE.MOLECULE;
  DG.DataFrame.fromColumns([col]);
  return col;
}

function molblock(smiles: string, explicitHydrogens = false): string {
  const mol = chemCommonRdKit.getRdKitModule().get_mol(smiles);
  try {
    return explicitHydrogens ? mol.add_hs() : mol.get_molblock();
  } finally {
    mol.delete();
  }
}

category('crux substructure search', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  after(async () => setSubstructureSearchEngine(null));

  test('matchesRdkit.smilesQueries', async () => {
    const df = await readDataframe('smiles.csv');
    await expectSameAsRdkit(df.col('canonical_smiles')!, QUERIES);
  });

  test('matchesRdkit.molblockQueries', async () => {
    const df = await readDataframe('smiles.csv');
    const queries = ['c1ccccc1', 'C1CNCCN1', 'c1ccc2[nH]ccc2c1', 'NS(=O)(=O)c1ccccc1', 'C[C@H](N)C(=O)O'];
    await expectSameAsRdkit(df.col('canonical_smiles')!, queries.map((q) => molblock(q)));
    await expectSameAsRdkit(df.col('canonical_smiles')!, queries.map((q) => molblock(q, true)));
  });

  test('matchesRdkit.notContainsMalformedAndEmpty', async () => {
    const df = await readDataframe('tests/Test_smiles_with_empty_and_malformed.csv');
    for (const searchType of [SubstructureSearchType.CONTAINS, SubstructureSearchType.NOT_CONTAINS])
      await expectSameAsRdkit(df.col('smiles')!, ['c1ccccc1', 'C1CC1', 'C(=O)N'], searchType);
  });

  test('matchesRdkit.molblockColumn', async () => {
    const df = await readDataframe('mol1K.csv');
    await expectSameAsRdkit(df.col('molecule')!, ['c1ccccc1', 'C(=O)N', 'c1ccc2[nH]ccc2c1', '[#7;R]']);
  });

  test('matchesRdkit.afterEdits', async () => {
    const df = await readDataframe('smiles.csv');
    const col = df.col('canonical_smiles')!;
    await expectSameAsRdkit(col, ['C1CC1']);
    col.set(0, 'C1CC1CCN');
    col.set(500, 'c1ccccc1');
    col.set(998, '');
    await expectSameAsRdkit(col, ['C1CC1', 'c1ccccc1']);
  }, {timeout: 60000});

  test('matchesRdkit.includeMask', async () => {
    const df = await readDataframe('smiles.csv');
    const col = df.col('canonical_smiles')!;
    const mask = await search(SubstructureSearchEngine.RDKit, col, 'c1ccccc1');
    await expectSameAsRdkit(col, ['c1ccncc1', 'C(=O)N'], SubstructureSearchType.CONTAINS, mask);
  });

  test('matchesRdkit.radicalLookingQueries', async () => {
    // typed, [OH] [CH3] [F] [N+] are read as SMARTS by both engines: the atom as written, whatever its neighbours;
    // so is c1cc[n+]cc1, which the radical would make a non-aromatic ring
    expect(await cruxReadsQuery('[OH]', ''), true);
    const col = (await readDataframe('smiles.csv')).col('canonical_smiles')!;
    expect((await search(SubstructureSearchEngine.RDKit, col, '[OH]')).trueCount() > 0, true);
    await expectSameAsRdkit(col, ['[OH]', '[CH3]', '[NH2]', '[F]', '[N+]', '[O-]', 'c1cc[n+]cc1']);
    const pyridinium = await search(SubstructureSearchEngine.RDKit,
      moleculeColumn(['C[n+]1ccccc1', 'c1cc[nH+]cc1', 'c1ccncc1']), 'c1cc[n+]cc1');
    expect([0, 1, 2].map((i) => pyridinium.getBit(i)).join(), 'true,true,false');
  });

  test('matchesRdkit.sanitizationCases', async () => {
    await expectSameAsRdkit(moleculeColumn(SANITIZATION_CASES),
      ['c1ccncc1', '[N+](=O)[O-]', '[O-]', 'N=[N+]=[N-]', 'Cl', '[Ra]', 'C1CCOC1', 'c1ccccc1', 'I']);
  });

  test('matchesRdkit.lenientCases', async () => {
    const col = moleculeColumn(LENIENT_CASES);
    for (const searchType of [SubstructureSearchType.CONTAINS, SubstructureSearchType.NOT_CONTAINS])
      await expectSameAsRdkit(col, ['c1ncnc1', 'c1ccccc1', 'P', 'n', 'C(=O)N'], searchType);
  });

  test('failedCruxFallsBackToRdkit', async () => {
    const df = grok.data.demo.molecules(20000);
    const col = df.col('smiles')!;
    col.semType = DG.SEMTYPE.MOLECULE;
    const rdkit = await search(SubstructureSearchEngine.RDKit, col, 'c1ccccc1');
    setSubstructureSearchEngine(SubstructureSearchEngine.Crux);
    const service = await getCruxService();
    const key = getSearchQueryAndType('c1ccccc1', SubstructureSearchType.CONTAINS, '', 0);
    let ended = false;
    const subs = [
      // the workers go away at the search's first update, as when one runs out of memory
      grok.events.onCustomEvent(getSearchProgressEventName(df.name, col.name)).pipe(first())
        .subscribe(() => service.restart(service.generation)),
      grok.events.onCustomEvent(getTerminateEventName(df.name, col.name)).subscribe((k) => ended ||= k === key),
    ];
    try {
      const result = await chemSubstructureSearchLibrary(col, 'c1ccccc1', 'c1ccccc1', FILTER_TYPES.substructure,
        false, false);
      await awaitCheck(() => ended, 'search has not ended', 60000);
      expect(result.equals(rdkit), true, `crux ${result}, RDKit ${rdkit}`);
    } finally {
      subs.forEach((s) => s.unsubscribe());
    }
  }, {timeout: 90000});

  test('cruxStartsNoRdkitWorker', async () => {
    // a Crux search on SMILES and on molblocks asks nothing of the RDKit worker pool (its init counts every ask), the
    // queries Chem sent to RDKit before included; a query Crux does not read does; the hits are RDKit's
    const smiles = (await readDataframe('smiles.csv')).col('canonical_smiles')!;
    const molblocks = (await readDataframe('mol1K.csv')).col('molecule')!;
    const before = rdkitPoolAsks();
    const crux = [await search(SubstructureSearchEngine.Crux, smiles, 'c1ccccc1'),
      await search(SubstructureSearchEngine.Crux, molblocks, 'c1ccccc1')];
    for (const query of ['[13CH3]O', '[C;v4]', '[#6;x2]', '[#6&R2]', '[C;r6]', '[#1]', '[!#6;!#1][#6]'])
      await search(SubstructureSearchEngine.Crux, smiles, query);
    expect(rdkitPoolAsks(), before, 'the RDKit workers were asked for by a Crux search');
    setSubstructureSearchEngine(SubstructureSearchEngine.Crux);
    // no query in either reading (getQueryMolSafe: none, nor its MOL block failover): RDKit's search reports it
    await chemSubstructureSearchLibrary(smiles, 'C1CC', 'not a molecule').catch(() => null);
    expect(rdkitPoolAsks() > before, true, 'a query Crux does not read asks for the RDKit workers');
    expect(crux[0].equals(await search(SubstructureSearchEngine.RDKit, smiles, 'c1ccccc1')), true, 'SMILES hits');
    expect(crux[1].equals(await search(SubstructureSearchEngine.RDKit, molblocks, 'c1ccccc1')), true, 'molblock hits');
  });

  test('fallsBackToRdkit', async () => {
    // Crux searches for every query its module reads as Chem's getQueryMolSafe builds it; for the others RDKit
    // searches, and its search reports what it reports
    setSubstructureSearchEngine(SubstructureSearchEngine.Crux);
    expect(await cruxReadsQuery('c1ccccc1', ''), true);
    expect(await cruxReadsQuery('', ''), true, 'the pre-search on filter attach');
    for (const [query, failover] of [['C1CC', 'not a molecule'], ['[C;', 'C1CC'], ['C'.repeat(5001), '']])
      expect(await cruxReadsQuery(query, failover), false, `Crux searches for ${query.slice(0, 20)}`);
    const col = (await readDataframe('smiles.csv')).col('canonical_smiles')!;
    const outcome = (s: Promise<BitArray>) => s.then((hits) => hits.toString(), (e) => `error: ${e?.message ?? e}`);
    setSubstructureSearchEngine(SubstructureSearchEngine.Crux);
    const viaCrux = await outcome(chemSubstructureSearchLibrary(col, 'C1CC', 'not a molecule'));
    setSubstructureSearchEngine(SubstructureSearchEngine.RDKit);
    const viaRdkit = await outcome(chemSubstructureSearchLibrary(col, 'C1CC', 'not a molecule'));
    expect(viaCrux, viaRdkit, 'a query Crux does not read');
  });

  test('matchesRdkit.formerlyRdkitQueries', async () => {
    // the queries Chem sent to RDKit while Crux took SMARTS text: through Crux, RDKit's hits, no RDKit worker started
    const smiles = (await readDataframe('smiles.csv')).col('canonical_smiles')!;
    const molblocks = (await readDataframe('mol1K.csv')).col('molecule')!;
    const special = moleculeColumn(SPECIAL_CASES);
    for (const col of [smiles, molblocks, special])
      await expectRdkitHitsThroughCrux(col, [...FORMERLY_RDKIT_QUERIES, ...FORMERLY_RDKIT_MOLBLOCKS]);
    for (const col of [smiles, special]) {
      await expectRdkitHitsThroughCrux(col, ['[#1]', '[#6;r5]', '[C;v4]', '[13#6]', '[CH2]', '[!#6][#6](=O)[!#6]'],
        SubstructureSearchType.NOT_CONTAINS);
    }
  }, {timeout: 180000});

  test('matchesRdkit.explicitHydrogenMolblocks', async () => {
    // a MOL block column whose molecules keep their H atoms (mol1K's 67 that have them, and 200 of smiles.csv with
    // all of them): Crux reads each MOL block as Chem's getMolSafe does, its H atoms kept, so H atoms, degrees (D),
    // implicit H (h) and connections (X) are RDKit's
    const mol1K = (await readDataframe('mol1K.csv')).col('molecule')!.toList() as string[];
    const smiles = (await readDataframe('smiles.csv')).col('canonical_smiles')!.toList() as string[];
    const withH = [...mol1K.filter((m) => / H {2}/.test(m)), ...smiles.slice(0, 200).map((s) => molblock(s, true))];
    expect(withH.length, 267, 'MOL blocks with H atoms');
    const col = moleculeColumn(withH);
    await expectRdkitHitsThroughCrux(col, ['[#6;D1]', '[#6;D4]', '[N;D1]', '[C;X4;H3]', 'c1ccccc1', 'C(=O)N', '[#1]',
      '[#6]-[#1]', '[#8]-[#1]', '[C;h3]', '[N;h2]', '[CH2]', '[#6;x2]', '[H]N([H])c1ccccc1', EXPLICIT_H_MOLBLOCK,
      molblock('CC(C)O', true)]);
    await expectRdkitHitsThroughCrux(col, ['[#1]', '[#6;D1]'], SubstructureSearchType.NOT_CONTAINS);
  }, {timeout: 120000});

  test('matchesRdkit.cxsmilesQueries', async () => {
    // CXSMILES extensions the Crux module reads as RDKit: coordinates, stereo groups, wedges, radicals, R-group,
    // attachment-point and alias labels; Crux searches them, with RDKit's hits
    const read = ['C[C@H](N)C(=O)O |&1:1|', 'C[C@H](N)C(=O)O |o1:1|', 'C[C@H](N)C(=O)O |&1:1,r|', 'C[C@H](N)O |w:1.0|',
      'CCO |(0,0,;1.5,0,;3,0,)|', 'CC(=O)N |(0,0,;1,0,;2,0,;3,0,),$;;;_AP1$|', 'CCO |$_R1;;$|', 'CC*.C |$;;OMe;$|',
      '[CH3]C |^1:0|'];
    const smiles = (await readDataframe('smiles.csv')).col('canonical_smiles')!;
    for (const col of [smiles, moleculeColumn(SPECIAL_CASES)])
      await expectRdkitHitsThroughCrux(col, read);
    // query features in the extension its reader does not read yet (query atom labels, positional variation, ring bond
    // count, unsaturation, substitution, bond lists): RDKit searches them, Crux mode gives RDKit's search
    const notRead = ['Cc1ccccc1 |$;;;;Q_e;;$|', 'Cc1ccccc1 |$;;;;AH_p;;$|', 'ClC |$X_p;$|',
      'c1ccccc1* |$;;;;;;star_e$|', '*c1ccccc1 |m:0:1.2.3|', 'C=CC |$;;$,rb:1:*|', 'CCC |u:1|', 'CCC |s:1:1|',
      'c1ccccc1 |c:0,2,4|'];
    for (const query of notRead) {
      expect(await cruxReadsQuery(query, ''), false, `Crux searches for ${query}`);
      const crux = await search(SubstructureSearchEngine.Crux, smiles, query);
      const rdkit = await search(SubstructureSearchEngine.RDKit, smiles, query);
      expect(crux.equals(rdkit), true, `Crux mode differs from RDKit: ${query} (${crux.trueCount()}, ` +
        `${rdkit.trueCount()})`);
    }
  }, {timeout: 120000});

  test('matchesRdkit.emptyQueryMolecule', async () => {
    // a query of no atoms (a MOL block of none; a query that does not read, with no failover: an empty molecule)
    // matches no molecule in RDKit; the pre-search (an empty query string) is not one
    const empty = ['', '  empty', '', '  0  0  0  0  0  0  0  0  0  0999 V2000', 'M  END', ''].join('\n');
    const col = moleculeColumn(SPECIAL_CASES);
    for (const searchType of [SubstructureSearchType.CONTAINS, SubstructureSearchType.NOT_CONTAINS])
      await expectRdkitHitsThroughCrux(col, [empty, 'C1CC'], searchType);
  });

  test('supersededSearchStops', async () => {
    setSubstructureSearchEngine(SubstructureSearchEngine.Crux);
    const df = grok.data.demo.molecules(20000);
    const col = df.col('smiles')!;
    col.semType = DG.SEMTYPE.MOLECULE;
    const keys = ['C1CC1', 'c1ccccc1'].map((q) => getSearchQueryAndType(q, SubstructureSearchType.CONTAINS, '', 0));
    const terminated = new Set<string>();
    let completed = 0;
    const subs = [
      grok.events.onCustomEvent(getSearchProgressEventName(df.name, col.name))
        .subscribe((p) => completed += Number(p) === 100 ? 1 : 0),
      grok.events.onCustomEvent(getTerminateEventName(df.name, col.name)).subscribe((key) => terminated.add(key)),
    ];
    try {
      // the second search replaces the first one before it gets to the workers
      await Promise.all(['C1CC1', 'c1ccccc1'].map((q) =>
        chemSubstructureSearchLibrary(col, q, q, FILTER_TYPES.substructure, false, false)));
      await awaitCheck(() => keys.every((k) => terminated.has(k)), 'searches have not ended', 30000);
    } finally {
      subs.forEach((s) => s.unsubscribe());
    }
    expect(completed, 1);
  });

  test('atMostSevenFilterUpdates', async () => {
    setSubstructureSearchEngine(SubstructureSearchEngine.Crux);
    const df = grok.data.demo.molecules(50000);
    const col = df.col('smiles')!;
    col.semType = DG.SEMTYPE.MOLECULE;
    const progress: number[] = [];
    let done = false;
    const progressEvents = grok.events.onCustomEvent(getSearchProgressEventName(df.name, col.name));
    const subs = [
      progressEvents.subscribe((p) => progress.push(Number(p))),
      grok.events.onCustomEvent(getTerminateEventName(df.name, col.name)).subscribe((key) => done ||= key !== ''),
    ];
    try {
      // not awaited to the end: the bit array fills in as the (cold, index building) search goes
      await chemSubstructureSearchLibrary(col, 'c1ccccc1', 'c1ccccc1', FILTER_TYPES.substructure, false, false);
      await awaitCheck(() => done, 'search has not finished', 60000);
    } finally {
      subs.forEach((s) => s.unsubscribe());
    }
    expect(progress.length >= 1 && progress.length <= 7, true, `${progress.length} filter updates`);
    expect(progress[progress.length - 1], 100);
  }, {timeout: 90000});

  test('substructureFilterWithEdits', async () => {
    const df = await readDataframe('sar-small.csv');
    await grok.data.detectSemanticTypes(df);
    const query = molblock('Clc1ccccc1');
    expect(await cruxReadsQuery(query, ''), true);
    const expected = (await search(SubstructureSearchEngine.RDKit, df.col('smiles')!, query)).trueCount();
    setSubstructureSearchEngine(SubstructureSearchEngine.Crux);
    const tv = grok.shell.addTableView(df);
    try {
      tv.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
        type: DG.FILTER_TYPE.SUBSTRUCTURE, column: 'smiles', columnName: 'smiles', molBlock: query}, false);
      await awaitCheck(() => df.filter.trueCount === expected, `df hasn't been filtered to ${expected} rows`, 15000);
      const row = df.filter.findNext(-1, false);
      df.col('smiles')!.set(row, 'NCCc1ccccc1Cl');
      await awaitCheck(() => df.filter.get(row) && df.filter.trueCount === expected + 1,
        'the edited molecule has not been re-checked', 15000);
    } finally {
      tv.close();
    }
  }, {timeout: 45000});

  // RDKit's SMARTS writer writes MDL's single or double with ring topology `-,=&@`, which a SMARTS parser, crux's too,
  // reads as `-,(=&@)` (crux-sketch spike query-roundtrip, L6: a custom bond query `-,=;@` passed 969 rows of mol1K for
  // 417); Crux takes the query as the molecule RDKit matches, no SMARTS text between them
  test('anAndOverAnOrBondQueryAsRdkitMatchesIt', async () => {
    setSubstructureSearchEngine(SubstructureSearchEngine.Crux);
    const molblock = AMBIGUOUS_BOND_MOLBLOCK;
    const df = await readDataframe('smiles.csv');
    await expectRdkitHitsThroughCrux(df.col('canonical_smiles')!, [molblock]);
    await expectSameAsRdkit(df.col('canonical_smiles')!, [molblock]);
    const rdkit = chemCommonRdKit.getRdKitModule();
    const intended = rdkit.get_qmol('[#6]-,=;@[#6]');
    const found = (q: any) => df.col('canonical_smiles')!.toList().filter((s: string) => {
      const m = rdkit.get_mol(s);
      try {
        return m.get_substruct_match(q) !== '{}';
      } finally {
        m.delete();
      }
    }).length;
    const query = rdkit.get_qmol(molblock);
    try {
      expect(found(query), found(intended), 'the molblock\'s query against [#6]-,=;@[#6]');
    } finally {
      query.delete();
      intended.delete();
    }
  });

  test('workerPoolIsBounded', async () => {
    // each worker holds a crux WebAssembly instance, and a page has room for only so many WebAssembly memories:
    // sized by the threads alone, crux's pool and RDKit's filled it on a 32-thread machine
    const service = await getCruxService();
    expect(service.workerCount, Math.max(1, Math.min(navigator.hardwareConcurrency - 2, 16)), 'crux workers');
  });
});
