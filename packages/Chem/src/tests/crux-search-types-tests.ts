import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {after, waitUntil, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {FILTER_TYPES, chemSubstructureSearchLibrary} from '../chem-searches';
import {SubstructureSearchType, getSearchProgressEventName, getSearchQueryAndType,
  getTerminateEventName} from '../constants';
import {Fingerprint} from '../utils/chem-common';
import {SubstructureSearchEngine, cruxReadsQuery, setSubstructureSearchEngine} from '../crux/crux-searches';
import {ChemicalEngine} from '../engine/chem-engine';
import {readDataframe} from './utils';
import {inBothEngines, rdkitPoolAsks} from './engine-utils';

// The search types beside Contains / Not contains on Crux in Crux mode (crux/crux-search-types.ts,
// crux/jobs/search-types.ts): the hits of Chem's RDKit search (on getMolSafe's molecules, as its pool caches them),
// no RDKit worker asked for.
const T = SubstructureSearchType;
const TYPES = [T.EXACT_MATCH, T.STEREO_AGNOSTIC, T.INCLUDED_IN, T.NOT_INCLUDED_IN];

/** The search engine of a chemical engine's mode. */
const searchEngine = (engine: ChemicalEngine): SubstructureSearchEngine =>
  engine === ChemicalEngine.Crux ? SubstructureSearchEngine.Crux : SubstructureSearchEngine.RDKit;

async function search(engine: SubstructureSearchEngine, col: DG.Column, query: string,
  searchType: SubstructureSearchType, cutoff = 0.8, fp = Fingerprint.Morgan): Promise<BitArray> {
  setSubstructureSearchEngine(engine);
  return chemSubstructureSearchLibrary(col, query, '', FILTER_TYPES.substructure, false, true, searchType, cutoff, fp);
}

const bits = (b: BitArray): string => Array.from({length: b.length}, (_, i) => b.getBit(i) ? '1' : '0').join('');

function molblock(smiles: string): string {
  const mol = chemCommonRdKit.getRdKitModule().get_mol(smiles);
  try {
    return mol.get_molblock();
  } finally {
    mol.delete();
  }
}

/** Each query of each type through Crux gives RDKit's hits, and the Crux searches ask nothing of the RDKit pool. */
async function expectRdkitHits(col: DG.Column, queries: string[], types = TYPES): Promise<void> {
  for (const searchType of types) {
    for (const query of queries) {
      const label = query.includes('M  END') ? `a MOL block (${query.split('\n').length} lines)` : query;
      expect(await cruxReadsQuery(query, ''), true, `Crux does not search for ${label}`);
      const rdkit = await search(SubstructureSearchEngine.RDKit, col, query, searchType);
      const asks = rdkitPoolAsks();
      const crux = await search(SubstructureSearchEngine.Crux, col, query, searchType);
      expect(rdkitPoolAsks(), asks, `${searchType} ${label}: the Crux search asked for the RDKit workers`);
      expect(crux.equals(rdkit), true, `${searchType} ${label} on ${col.name}: Crux ${crux.trueCount()} rows, RDKit ` +
        `${rdkit.trueCount()} (the first differing row ${crux.clone().xor(rdkit).findNext(-1, true)})`);
    }
  }
}

category('crux substructure search: search types', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  after(async () => setSubstructureSearchEngine(null));

  test('matchesRdkit.smilesColumn', async () => {
    // rows of the column as SMILES and as MOL blocks (stereo from their wedges), smaller queries, a SMARTS (Exact: each
    // a substructure of the other)
    const col = (await readDataframe('smiles.csv')).col('canonical_smiles')!;
    const rows = (col.toList() as string[]).filter((s, i) => i % 97 === 3).slice(0, 6);
    await expectRdkitHits(col, [...rows, ...rows.slice(0, 3).map((s) => molblock(s)), 'c1ccccc1', 'C[C@H](N)C(=O)O',
      'CC(=O)Nc1ccc(O)cc1', '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1']);
  }, {timeout: 120000});

  test('matchesRdkit.molblockColumn', async () => {
    const col = (await readDataframe('mol1K.csv')).col('molecule')!;
    const rows = (col.toList() as string[]).filter((_, i) => i % 211 === 5).slice(0, 4);
    await expectRdkitHits(col, rows, [T.EXACT_MATCH, T.STEREO_AGNOSTIC, T.INCLUDED_IN]);
  }, {timeout: 120000});

  test('matchesRdkit.stereo', async () => {
    // stereoisomers and their flat forms: Exact tells them apart, Stereo agnostic does not
    const col = (await readDataframe('stereo-agnostic-test.csv')).col('smiles')!;
    await expectRdkitHits(col, ['C[C@@H](O)F', 'F/C=C/Cl', 'O[C@@H]1CC[C@H](O)CC1', molblock('C[C@H](O)F'), 'CC(O)F'],
      [T.EXACT_MATCH, T.STEREO_AGNOSTIC]);
  });

  test('rowsCruxDoesNotReadByRdkit', async () => {
    // a CXSMILES whose extension Crux does not read (an atom value): RDKit answers for its row; empty cells match none
    const col = DG.Column.fromStrings('smiles', ['CCO |$_AV:1;;$|', 'CCO', '', 'c1ccccc1O', 'CCN']);
    col.semType = DG.SEMTYPE.MOLECULE;
    DG.DataFrame.fromColumns([col]);
    await expectRdkitHits(col, ['CCO', 'c1ccccc1OCCO']);
    const exact = await search(SubstructureSearchEngine.Crux, col, 'CCO', T.EXACT_MATCH);
    expect(bits(exact), '11000');
  });

  test('keysFollowExtensionEdits', async () => {
    // an edit of a CXSMILES extension alone keeps the index's segment (it leaves extensions out): the keys follow it
    const col = DG.Column.fromStrings('smiles', ['C', 'CC']);
    col.semType = DG.SEMTYPE.MOLECULE;
    DG.DataFrame.fromColumns([col]);
    expect(bits(await search(SubstructureSearchEngine.Crux, col, 'C', T.EXACT_MATCH)), '10');
    col.set(0, 'C |^1:0|'); // a radical: [CH3]
    await expectRdkitHits(col, ['C'], [T.EXACT_MATCH]);
    expect(bits(await search(SubstructureSearchEngine.Crux, col, 'C', T.EXACT_MATCH)), '00');
  });

  test('similarOnCrux', async () => {
    // Crux mode's fingerprints (the engine's cache) and Chem's Tanimoto: RDKit mode's hits for Morgan and MACCS
    const df = await readDataframe('smiles.csv');
    const col = df.col('canonical_smiles')!;
    const query = col.get(250);
    const settings: [Fingerprint, number][] = [[Fingerprint.Morgan, 0.3], [Fingerprint.MACCS, 0.6]];
    const {rdkit, crux} = await inBothEngines(async (engine) => {
      const hits: BitArray[] = [];
      for (const [fp, cutoff] of settings)
        hits.push(await search(searchEngine(engine), col, query, T.IS_SIMILAR, cutoff, fp));
      return hits;
    }, {noRdkitPool: true});
    expect(!!df.col('~canonical_smiles.Morgan.Crux'), true, 'Crux\'s fingerprints cached');
    settings.forEach(([fp], k) => {
      expect(crux[k].trueCount() > 0, true, `${fp}: no hits`);
      expect(bits(crux[k]), bits(rdkit[k]), `${fp}: Crux ${crux[k].trueCount()} rows, RDKit ${rdkit[k].trueCount()}`);
    });
  }, {timeout: 120000});

  test('filterSearchEnds', async () => {
    // not awaited, as the filter searches: the bit array fills in, the search reports 100% and ends with its key
    setSubstructureSearchEngine(SubstructureSearchEngine.Crux);
    const df = grok.data.demo.molecules(20000);
    const col = df.col('smiles')!;
    col.semType = DG.SEMTYPE.MOLECULE;
    const query = col.get(77);
    const expected = await search(SubstructureSearchEngine.Crux, col, query, T.INCLUDED_IN);
    const key = getSearchQueryAndType(query, T.INCLUDED_IN, Fingerprint.Morgan, 0.8);
    let ended = false;
    let completed = false;
    const subs = [
      grok.events.onCustomEvent(getSearchProgressEventName(df.name, col.name))
        .subscribe((p) => completed ||= p === 100),
      grok.events.onCustomEvent(getTerminateEventName(df.name, col.name)).subscribe((k) => ended ||= k === key),
    ];
    try {
      const result = await chemSubstructureSearchLibrary(col, query, query, FILTER_TYPES.substructure, false, false,
        T.INCLUDED_IN);
      await waitUntil(() => ended && completed, 'the search has not ended', 30000);
      expect(result.equals(expected), true, `${result.trueCount()} rows, awaited ${expected.trueCount()}`);
    } finally {
      subs.forEach((s) => s.unsubscribe());
    }
  }, {timeout: 60000});

  test('timing on 100,000 rows', async () => {
    // each type's first search (Crux: the segments' keys; RDKit: its fingerprints) and a second one, both engines, on a
    // fresh 100,000-row column each; a measurement: with --benchmark only
    if (!DG.Test.isInBenchmark)
      return;
    const lines: string[] = [];
    await inBothEngines(async (engine) => {
      for (const searchType of [...TYPES, T.IS_SIMILAR]) {
        const col = grok.data.demo.molecules(100000).col('smiles')!;
        col.semType = DG.SEMTYPE.MOLECULE;
        const queries = [col.get(1234), col.get(56789)];
        const times = [];
        for (const query of queries) {
          const t0 = performance.now();
          const hits = await search(searchEngine(engine), col, query, searchType, 0.6);
          times.push(`${Math.round(performance.now() - t0)} ms (${hits.trueCount()} hits)`);
        }
        lines.push(`${engine} ${searchType}: first ${times[0]}, second ${times[1]}`);
      }
    });
    console.log(`search types timing | ${lines.join('; ')}`);
    return lines.join('; ');
  }, {timeout: 600000, benchmark: true, benchmarkTimeout: 1800000});
});
