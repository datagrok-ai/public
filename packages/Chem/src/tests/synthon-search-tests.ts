import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';
import {readDataframe} from './utils';
import {ChemicalEngine} from '../engine/chem-engine';
import {expectNoFallback, inEngine} from './engine-utils';
import {PackageFunctions} from '../package';

const MOLECULE = 'Cc1ccc2c(C(NC)=O)n[nH]c2c1';
const SPACE = 'Syntons_5567.csv';

function compareDfColumns(result: DG.DataFrame, expected: DG.DataFrame): void {
  for (const col of expected.columns)
    expect(result.columns.contains(col.name), true, `missing column: ${col.name}`);
}

function compareProductColumn(result: DG.DataFrame, expected: DG.DataFrame): void {
  const resProducts = new Set<string>();
  const resCol = result.getCol('product');
  for (let i = 0; i < result.rowCount; i++)
    resProducts.add(resCol.get(i));

  const expCol = expected.getCol('product');
  for (let i = 0; i < expected.rowCount; i++)
    expect(resProducts.has(expCol.get(i)), true, `missing product: ${expCol.get(i)}`);
}

/**
 * Chem:SynthonSearch's hits against synthon_search.py's (files/tests/Synthon_*_Search_Results.csv, the script on its
 * RDKit 2025.09.3): the products in order, their reactions and synthons, the similarities (Datagrok's float).
 */
function compareWithScript(result: DG.DataFrame, expected: DG.DataFrame, label: string): void {
  expect(result.rowCount, expected.rowCount, `${label}: rows`);
  const pairs: [string, string][] = [['smiles', 'product'], ['reaction_id', 'reaction_id'],
    ['synthon_1_id', 'precursor_1_id'], ['synthon_2_id', 'precursor_2_id']];
  for (const [got, want] of pairs) {
    for (let i = 0; i < expected.rowCount; i++)
      expect(String(result.get(got, i)), String(expected.get(want, i)), `${label}: ${got} of row ${i}`);
  }
  if (expected.columns.contains('similarity')) {
    for (let i = 0; i < expected.rowCount; i++) {
      const [got, want] = [result.get('similarity', i), expected.get('similarity', i)];
      expect(Math.abs(got - want) < 1e-6, true, `${label}: similarity of row ${i}: ${got} (the script's ${want})`);
    }
  }
}

/** Chem:SynthonSearch (and Search Synthons through it) in the engine given, against the script's hits; each call's time
 * logged. */
async function checkSynthonSearch(engine: ChemicalEngine): Promise<void> {
  await inEngine(engine, async () => {
    const text = await _package.files.readAsText(`synthon-data/${SPACE}`);
    const library = DG.FileInfo.fromString(SPACE, text);
    const times: string[] = [];
    for (const [searchType, file] of [['substructure', 'Synthon_Substructure_Search_Results.csv'],
      ['similarity', 'Synthon_Similarity_Search_Results.csv']]) {
      const started = performance.now();
      // the bundle's own function: the engine set here is this bundle's (grok.functions.call runs package.js's)
      const result = await expectNoFallback(['synthonSearch'], () =>
        // eslint-disable-next-line new-cap -- the function keeps the script's name
        PackageFunctions.SynthonSearch(MOLECULE, library, SPACE, 100, searchType, 0.5, true));
      times.push(`${searchType} ${Math.round(performance.now() - started)} ms`);
      compareWithScript(result, await readDataframe(`tests/${file}`), searchType);
      console.log(`synthon search | ${engine}, ${searchType} columns: ${result.columns.toList()
        .map((c) => `${c.name} ${c.type} ${c.semType}`).join(', ')}`);
      // a synthon's structure is the library's
      const id = String(result.get('synthon_1_id', 0));
      const line = text.split('\n').find((l) => l.split(',')[1] === id);
      expect(result.get('synthon_1', 0), line?.split(',')[0], `${searchType}: synthon_1 of row 0`);
    }
    // the panels' call, Search Synthons, through it: package.js's, which follows the stand's engine (not this bundle's)
    if (await grok.functions.call('Chem:chemicalEngine') === engine) {
      const started = performance.now();
      const panel: DG.DataFrame = await grok.functions.call('Chem:SynthonSearchFunc', {spaceName: SPACE,
        molecule: MOLECULE, maxHits: 100, searchType: 'similarity', similarityCutoff: 0.5, includeSynthons: true});
      times.push(`Search Synthons (similarity) ${Math.round(performance.now() - started)} ms`);
      compareWithScript(panel, await readDataframe('tests/Synthon_Similarity_Search_Results.csv'), 'Search Synthons');
    }
    console.log(`synthon search | ${engine}, Chem's function end to end: ${times.join(', ')}`);
  });
}

category('synthon search', () => {
  test('Chem:SynthonSearch is Chem\'s function, the script SynthonSearchPython', async () => {
    const func = DG.Func.find({package: 'Chem', name: 'SynthonSearch'});
    expect(func.length, 1, 'one Chem:SynthonSearch');
    expect(func[0] instanceof DG.Script, false, 'Chem:SynthonSearch is the package function');
    const script = DG.Func.find({package: 'Chem', name: 'SynthonSearchPython'});
    expect(script.length === 1 && script[0] instanceof DG.Script, true, 'the script is SynthonSearchPython');
  });

  // red / green for the Chemical Engine: in Crux mode the Crux workers answer (no fallback to the script)
  test('Crux: Chem:SynthonSearch gives the script\'s hits', async () => {
    await checkSynthonSearch(ChemicalEngine.Crux);
  }, {timeout: 600000});

  // the script itself (RDKit mode, and Crux's fallback): its kernel times out on the stand, as the tests below
  test('RDKit: Chem:SynthonSearch is the script', async () => {
    await checkSynthonSearch(ChemicalEngine.RDKit);
  }, {timeout: 600000, skipReason: 'GROK-19863'});

  test('substructure search', async () => {
    const result: DG.DataFrame = await grok.functions.call('Chem:synthonSearchFunc', {
      spaceName: SPACE,
      molecule: MOLECULE,
      maxHits: 100,
      searchType: 'substructure',
    });
    expect(result instanceof DG.DataFrame, true, 'result should be a DataFrame');
    expect(result.rowCount > 0, true, 'should return at least one hit');

    const expected = await readDataframe('tests/Synthon_Substructure_Search_Results.csv');
    compareDfColumns(result, expected);
    expect(result.rowCount, expected.rowCount, 'row count mismatch');
    compareProductColumn(result, expected);
  }, {timeout: 120000, skipReason: 'GROK-19863'});

  test('similarity search', async () => {
    const result: DG.DataFrame = await grok.functions.call('Chem:synthonSearchFunc', {
      spaceName: SPACE,
      molecule: MOLECULE,
      maxHits: 100,
      searchType: 'similarity',
      similarityCutoff: 0.5,
    });
    expect(result instanceof DG.DataFrame, true, 'result should be a DataFrame');
    expect(result.rowCount > 0, true, 'should return at least one hit');

    const expected = await readDataframe('tests/Synthon_Similarity_Search_Results.csv');
    compareDfColumns(result, expected);
    expect(result.rowCount, expected.rowCount, 'row count mismatch');
    compareProductColumn(result, expected);

    const simCol = result.getCol('similarity');
    for (let i = 0; i < result.rowCount; i++)
      expect(simCol.get(i) >= 0.5, true, `similarity below cutoff at row ${i}`);
  }, {timeout: 120000, skipReason: 'GROK-19863'});
});
