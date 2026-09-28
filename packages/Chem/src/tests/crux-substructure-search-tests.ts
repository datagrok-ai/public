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
import {SubstructureSearchEngine, getCruxQuery, getCruxService, setSubstructureSearchEngine}
  from '../crux/crux-searches';
import {readDataframe} from './utils';

const QUERIES = ['c1ccccc1', 'c1ccncc1', 'C1CC1', 'C(=O)N', 'c1ccc2[nH]ccc2c1', 'C[C@H](N)C(=O)O', 'C/C=C/C',
  '[CX3](=O)[OX2H1]', '[NX3;H2,H1;!$(NC=O)]', '[F,Cl,Br,I]', '[#7;R]', '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-,=[#8]'];
// isotopes and SMARTS primitives crux does not evaluate like RDKit
const RDKIT_ONLY_QUERIES = ['[13CH3]O', '[C;v4]', '[#6;x2]', '[#6&R2]', '[C;r6]'];
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
      expect(await getCruxQuery(query, '') !== null, true, `crux does not search for ${query}`);
    const rdkit = await search(SubstructureSearchEngine.RDKit, col, query, searchType, includeMask);
    const crux = await search(SubstructureSearchEngine.Crux, col, query, searchType, includeMask);
    expect(crux.equals(rdkit), true, `crux differs from RDKit: ${searchType} ${query} (crux ${crux}, RDKit ${rdkit})`);
  }
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
    expect(await getCruxQuery('[OH]', ''), '[O&H1]');
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

  test('fallsBackToRdkit', async () => {
    setSubstructureSearchEngine(SubstructureSearchEngine.Crux);
    for (const query of RDKIT_ONLY_QUERIES)
      expect(await getCruxQuery(query, ''), null, `crux should not handle ${query}`);
    expect(await getCruxQuery('c1ccccc1', ''), '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1');
    const df = await readDataframe('smiles.csv');
    await expectSameAsRdkit(df.col('canonical_smiles')!, RDKIT_ONLY_QUERIES, SubstructureSearchType.CONTAINS, null,
      false);
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
    expect(await getCruxQuery(query, '') !== null, true);
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
});
