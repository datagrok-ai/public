import {category, test, expect, before, testEvent, expectArray, waitUntil} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {SubstructureFilter} from '../widgets/chem-substructure-filter';
import {readDataframe} from './utils';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import { chemSimilaritySearch } from '../analysis/chem-similarity-viewer';
import { BitArrayMetrics } from '@datagrok-libraries/ml/src/typed-metrics';
import { Fingerprint } from '../utils/chem-common';
import { SubstructureSearchType } from '../constants';
import { sketchersWarmUp } from './sketcher-tests';
import {SubstructureSearchEngine, setSubstructureSearchEngine} from '../crux/crux-searches';

const expectedResults: {[key: string]: any} = {
  'oneColumn': [737141248, 593097, 3256025153, 4],
  'malformed': [3029331455, 16127],
  'empty': [524214],
};

const molblock1 = `
Actelion Java MolfileCreator 1.0

  7  7  0  0  0  0  0  0  0  0999 V2000
   12.9197  -11.7188   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.9197  -13.2187   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   14.2188  -13.9688   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   15.5178  -13.2187   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   15.5178  -11.7188   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   14.2188  -10.9688   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   14.2188   -9.4687   -0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  6  7  1  0  0  0  0
M  END
`;

const molblock2 = `
Actelion Java MolfileCreator 1.0

 23 25  0  0  0  0  0  0  0  0999 V2000
    8.2116  -11.0868    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.0566  -12.3261    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.5524  -12.2140    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.3974  -13.4534    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.8932  -13.3413    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   13.6432  -14.6403    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.7982  -15.8797    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   15.1264  -14.8639    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   16.2259  -13.8436    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   16.1138  -12.3478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   14.8745  -11.5028    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   15.2083  -10.0404    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   14.1088   -9.0202    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   14.4425   -7.5578    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   12.6754   -9.4623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.3416  -10.9247    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   13.4412  -11.9450    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   17.4129  -11.5978    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   18.7119  -12.3478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   20.0109  -11.5978    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   20.0109  -10.0978    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   18.7119   -9.3478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   17.4129  -10.0978    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  6  7  2  0  0  0  0
  6  8  1  0  0  0  0
  8  9  1  0  0  0  0
  9 10  2  0  0  0  0
 10 11  1  0  0  0  0
 11 12  2  0  0  0  0
 12 13  1  0  0  0  0
 13 14  1  0  0  0  0
 13 15  2  0  0  0  0
 15 16  1  0  0  0  0
 16 17  2  0  0  0  0
 10 18  1  0  0  0  0
 18 19  1  0  0  0  0
 19 20  1  0  0  0  0
 20 21  1  0  0  0  0
 21 22  1  0  0  0  0
 22 23  1  0  0  0  0
 17  5  1  0  0  0  0
 17 11  1  0  0  0  0
 23 18  1  0  0  0  0
M  END
`;

const molFileForCloneTest2 = `
MJ201900                      

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.1786    0.8920    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8930    0.4795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8930   -0.3455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1786   -0.7580    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5358   -0.3455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5358    0.4795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
M  END
`;

const molFileForCloneTest1 = `
MJ201900                      

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.6919    0.4455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4064    0.0330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4064   -0.7920    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6919   -1.2044    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0225   -0.7920    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0225    0.0330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  6  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
M  END
`;

/** Toluene drawn with its ring Kekulé and its methyl bond aromatic (MDL type 4), as a sketcher writes it once the user
 * marks that bond aromatic (the product owner's report of 2026-10-07, crux-sketch spike query-roundtrip). */
const TOLUENE_AROMATIC_METHYL_BOND = `
  query-roundtrip

  7  7  0  0  0  0  0  0  0  0999 V2000
    2.5714    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0714    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3214   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1786   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9286    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1786    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3214    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  2  0
  7  2  1  0
M  END
`;

category('substructure filters', async () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    const funcs = DG.Func.find({meta: {role: DG.FUNC_TYPES.MOLECULE_SKETCHER}});
    await sketchersWarmUp(funcs);
  });

  // The report's rows: with Filter as you draw on, toluene filtered first, then its methyl bond marked aromatic. The
  // filter named a query by RDKit's sanitized SMARTS, which reads that bond single, so the second query had the first's
  // name and its search was skipped: toluene's 822 rows of mol1K stayed, where the aromatic bond passes 293.
  test('a query that only an unsanitized reading tells from the last one is searched again: toluene, ' +
    'then its methyl bond aromatic', async () => {
    const df = await readDataframe('mol1K.csv');
    await grok.data.detectSemanticTypes(df);
    const dialogs: DG.Dialog[] = [];
    const filter = await createFilter('molecule', df, dialogs);
    try {
      const rdkit = chemCommonRdKit.getRdKitModule();
      const matches = (smarts: string): number => {
        const q = rdkit.get_qmol(smarts);
        let n = 0;
        try {
          for (let i = 0; i < df.rowCount; i++) {
            const m = rdkit.get_mol(df.get('molecule', i));
            try {
              if (m.get_substruct_match(q) !== '{}')
                n++;
            } finally {
              m.delete();
            }
          }
        } finally {
          q.delete();
        }
        return n;
      };
      const toluene = matches('[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1');
      const marked = matches('[#6]:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1');
      expect(toluene !== marked && marked > 0, true, `mol1K tells the two apart (${toluene}, ${marked})`);
      // the same written V3000 (RDKit's writer would kekulize the bond away)
      const v3000 = ['', '  query-roundtrip', '', '  0  0  0     0  0            999 V3000', 'M  V30 BEGIN CTAB',
        'M  V30 COUNTS 7 7 0 0 0', 'M  V30 BEGIN ATOM', 'M  V30 1 C 2.5714 0 0 0', 'M  V30 2 C 1.0714 0 0 0',
        'M  V30 3 C 0.3214 -1.299 0 0', 'M  V30 4 C -1.1786 -1.299 0 0', 'M  V30 5 C -1.9286 0 0 0',
        'M  V30 6 C -1.1786 1.299 0 0', 'M  V30 7 C 0.3214 1.299 0 0', 'M  V30 END ATOM', 'M  V30 BEGIN BOND',
        'M  V30 1 4 1 2', 'M  V30 2 2 2 3', 'M  V30 3 1 3 4', 'M  V30 4 2 4 5', 'M  V30 5 1 5 6', 'M  V30 6 2 6 7',
        'M  V30 7 1 7 2', 'M  V30 END BOND', 'M  V30 END CTAB', 'M  END', ''].join('\n');
      // with either engine of the search, and the same query written V3000 (the search reads it unsanitized too)
      for (const engine of [SubstructureSearchEngine.RDKit, SubstructureSearchEngine.Crux]) {
        setSubstructureSearchEngine(engine);
        for (const [what, query] of [['V2000', TOLUENE_AROMATIC_METHYL_BOND], ['V3000', v3000]]) {
          filter.sketcher.setMolFile(TOLUENE_AROMATIC_METHYL_BOND.replace('  1  2  4  0', '  1  2  1  0'));
          await waitUntil(() => df.filter.trueCount === toluene,
            `${engine}: toluene: ${toluene} rows contain it`, 15000);
          filter.sketcher.setMolFile(query);
          await waitUntil(() => df.filter.trueCount === marked, 'see below', 15000).catch(() => {
            const kept = filter.moleculeToSmarts(filter.currentMolecule);
            const given = filter.moleculeToSmarts(filter.getSketcherMolecule());
            throw new Error(`with ${engine} and the molblock ${what}, ${df.filter.trueCount} rows pass, ` +
              `${marked} contain it; the filter keeps ${kept}, its sketcher ` +
              `(${filter.sketcher.sketcher?.constructor.name}) gives ${given}; ` +
              `searches ${filter.currentSearches.size}, sync ${filter.syncEvent}, calculating ${filter.calculating}`);
          });
          expect(filter.currentMolecule, query, `${engine}, ${what}: the molecule the filter keeps`);
        }
      }
    } finally {
      setSubstructureSearchEngine(null);
      dialogs.forEach((d) => d.close());
      filter.detach();
    }
  });

  test('filterBy2Columns', async () => {
    const df = await readDataframe('tests/smiles_2_columns.csv');
    await grok.data.detectSemanticTypes(df);
    const sketcherDialogs: DG.Dialog[] = [];

    const filter1 = await createFilter('smiles1', df, sketcherDialogs);
    const filter2 = await createFilter('smiles2', df, sketcherDialogs);
    const molfile1 = `
      MJ201900                      

  5  5  0  0  0  0  0  0  0  0999 V2000
   -0.1785   -0.1946    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8930    0.2178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8930    1.0428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5358    1.0428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5358    0.2178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  5  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
M  END`;
    const molfile2 = `
      MJ201900                      

  5  5  0  0  0  0  0  0  0  0999 V2000
   -0.1785   -0.1946    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8930    0.2178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8930    1.0428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5358    1.0428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5358    0.2178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  1  2  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
`;

    filter1.sketcher.setMolFile(molfile1);
    filter2.sketcher.setMolFile(molfile2);
    await waitUntil(() => df.filter.trueCount === 2 && df.filter.get(4) && df.filter.get(5) && !df.filter.get(0),
      'df hasn\'t been filtered', 7000);
    await searchEnded(filter1);
    await searchEnded(filter2);
    sketcherDialogs.forEach((it) => it.close());
    filter1.detach();
    filter2.detach();
  });

  test('filterByOneColumn', async () => {
    await testOneColumn('tests/spgi-100.csv', 'Structure', 'c1ccccc1',
      'terminate_substructure_search-tests/spgi-100-Structure', 'oneColumn', 32);
  });

  test('malformed_filterByOneColumn', async () => {
    await testOneColumn('tests/Test_smiles_malformed.csv', 'canonical_smiles', 'c1ccccc1',
      'terminate_substructure_search-tests/Test_smiles_malformed-canonical_smiles', 'malformed', 36);
  });

  test('empty_filterByOneColumn', async () => {
    await testOneColumn('tests/sar-small_empty_vals.csv', 'smiles', 'C1CCCCC1',
      'terminate_substructure_search-tests/sar-small_empty_vals-smiles', 'empty', 16);
  });

  test('terminateOneSearchByAnother', async () => {
    const df = await readDataframe('tests/smi10K.csv');
    await grok.data.detectSemanticTypes(df);
    const sketcherDialogs: DG.Dialog[] = [];

    const filter = await createFilter('smiles', df, sketcherDialogs);
    const substr1 = 'C1CCCCC1';
    const substr2 = 'CC1CCCCC1';

    //starting filtering by 1st structure
    filter.sketcher.setSmiles(substr1);
    //terminating filtering by 1st structure
    await waitUntil(() => df.filter.trueCount !== df.rowCount,
      'filtering by 1st structure hasn\'t been started', 60000);
    filter.sketcher.setSmiles(substr2);
    //finishing filtering by 2nd structure
    await waitUntil(() => df.filter.trueCount === 286, 'filtering by 1st structure hasn\'t been started', 60000);

    await searchEnded(filter);
    sketcherDialogs.forEach((it) => it.close());
    filter.detach();
  }, {timeout: 60000});

  test('filteringMultipleDfs', async () => {
    const df1 = await readDataframe('tests/smi10K.csv');
    const df2 = await readDataframe('tests/smi10K.csv');
    df2.name = 'tests/smi10K (2)';
    await grok.data.detectSemanticTypes(df1);
    await grok.data.detectSemanticTypes(df2);
    const sketcherDialogs: DG.Dialog[] = [];

    const filter1 = await createFilter('smiles', df1, sketcherDialogs);
    const filter2 = await createFilter('smiles', df2, sketcherDialogs);

    const substr1 = 'C1CCCCC1';
    const substr2 = 'CC1CCCCC1';

    //starting filtering by 1st structure
    filter1.sketcher.setSmiles(substr1);
    await waitUntil(() => filter1.calculating || df1.filter.trueCount !== df1.rowCount,
      'filtering of df1 didn\'t start', 30000);
    //opening 2nd filter, sketching a molecule and waiting for filters to complete search
    filter2.sketcher.setSmiles(substr2);
    await waitUntil(() => df1.filter.trueCount === 462, 'df1 hasn\'t been filtered', 30000);
    await waitUntil(() => df2.filter.trueCount === 286, 'df2 hasn\'t been filtered', 30000);
    await searchEnded(filter1);
    await searchEnded(filter2);
    sketcherDialogs.forEach((it) => it.close());
    filter1.detach();
    filter2.detach();
  }, {timeout: 60000});

  test('multipleDfsWithTerminatedSearch', async () => {
    DG.chem.currentSketcherType = 'OpenChemLib';
    const df1 = await readDataframe('tests/smi10K.csv');
    const df2 = await readDataframe('tests/smi10K.csv');
    df2.name = 'tests/smi10K (2)';
    await grok.data.detectSemanticTypes(df1);
    await grok.data.detectSemanticTypes(df2);
    const sketcherDialogs: DG.Dialog[] = [];

    const filter1 = await createFilter('smiles', df1, sketcherDialogs);
    const filter2 = await createFilter('smiles', df2, sketcherDialogs);
    const substr1 = 'C1CCCCC1';
    const substr2 = 'c1ccccc1';
    const substr3 = 'CC1CCCCC1';

    //starting filtering by 1st structure
    filter1.sketcher.setSmiles(substr1);
    await waitUntil(() => df1.filter.trueCount !== 10001, 'filtering of df1 didn\'t start', 30000);
    //setting 1st structure to the 2nd filter
    filter2.sketcher.setSmiles(substr2);
    //setting 2nd structure to the 2nd filter
    filter2.sketcher.setSmiles(substr3);
    //waiting for 2nd filter to complete
    await waitUntil(() => df1.filter.trueCount === 462, 'df1 hasn\'t been filtered correctly', 30000);
    await waitUntil(() => df2.filter.trueCount === 286, 'df2 hasn\'t been filtered correctly', 30000);

    sketcherDialogs.forEach((it) => it.close());
    filter1.detach();
    filter2.detach();
  }, {timeout: 60000});

  test('similaritySearchAfterTerminatedSearch', async () => { //#2533 (https://github.com/datagrok-ai/public/issues/2533)
    const df = await readDataframe('tests/smi10K.csv');
    await grok.data.detectSemanticTypes(df);
    const sketcherDialogs: DG.Dialog[] = [];

    const filter = await createFilter('smiles', df, sketcherDialogs);
    const substr1 = 'C1CCCCC1';
    const substr2 = DG.WHITE_MOLBLOCK;

    //start filtering by 1st structure
    filter.sketcher.setSmiles(substr1);
    //waiting for first filtering results
    await waitUntil(() => df.filter.trueCount !== df.rowCount, 'filtering hasn\'t been started', 10000);
    //terminate filtering
    filter.sketcher.setMolFile(substr2);
    //finish filtering
    await waitUntil(() => df.filter.trueCount === df.rowCount, 'filter hasn\'t been reset', 10000);
    const simResults = await chemSimilaritySearch(df.col('smiles')!, df.get('smiles', 0),
      'Tanimoto' as BitArrayMetrics, 12, 0.01, 'Morgan' as Fingerprint);
    expect(simResults?.get('indexes', 0), 0);
    expect(simResults?.get('indexes', 5), 4002);
    expect(simResults?.get('indexes', 11), 731);
    sketcherDialogs.forEach((it) => it.close());
    filter.detach();
  });

  test('filterOptionsSynchronization', async () => { //#2512 (https://github.com/datagrok-ai/public/issues/2512)
    const df = await readDataframe('tests/spgi-100.csv');
    await grok.data.detectSemanticTypes(df);
    const sketcherDialogs: DG.Dialog[] = [];

    const filter1 = await createFilter('Structure', df, sketcherDialogs);
    const filter2 = await createFilter('Structure', df, sketcherDialogs);
    const substr = 'C1CCCCC1';
    const terminateFlag = 'terminate_substructure_search-tests/spgi-100-Structure';

    //filter by structure and wait for results
    filter1.sketcher.setSmiles(substr);
    await waitUntil(() => df.filter.trueCount === 5, `df hasn't been filtered during 5000 ms`, 5000);
    await searchEnded(filter1);
    //check that filter structure is synchronized in second filter
    await waitUntil(() => filter2.sketcher.getSmiles() === 'C1CCCCC1',
      () => `structure is not synchronized between filters: ${filter2.sketcher.getSmiles()}`);
    //change filtering options in second filter and wait for results
    filter2.searchTypeInput.value = SubstructureSearchType.NOT_CONTAINS;
    await waitUntil(() => df.filter.trueCount === 95 && !filter2.calculating && filter2.currentSearches.size === 0,
      () => `df hasn't been filtered by Not contains: ${df.filter.trueCount} rows pass`);
    //check that filter options are synchronized
    expect(filter1.searchType, SubstructureSearchType.NOT_CONTAINS, 'filter options not synchronized between filters');
    sketcherDialogs.forEach((it) => it.close());
    filter1.detach();
    filter2.detach();
  });

  test('properSearchFinish', async () => { //#2400 (https://github.com/datagrok-ai/public/issues/2400)
    const df = await readDataframe('tests/spgi-100.csv');
    await grok.data.detectSemanticTypes(df);
    const sketcherDialogs: DG.Dialog[] = [];

    DG.chem.currentSketcherType = 'Ketcher';
    const filter1 = await createFilter('Structure', df, sketcherDialogs);

    //filter by structure and wait for results
    filter1.sketcher.setSmiles('C1CCCCC1');
    await waitUntil(() => df.filter.trueCount === 5, 'df hasn\'t been filtered', 3000);

    //check that search is finished and loader is disabled (the rows of the last batch are filtered before it finishes)
    await waitUntil(() => !filter1.calculating, 'search hasn\'t been finished properly, loader is active', 3000);
    sketcherDialogs.forEach((it) => it.close());
    filter1.detach();
    DG.chem.currentSketcherType = 'OpenChemLib';
  });

  test('detachingFilters', async () => { //check that in case multiple filters are applied to the same column (like when cloning the view), we can detach one, but the df will stay filtered
    const df = await readDataframe('tests/spgi-100.csv');
    await grok.data.detectSemanticTypes(df);
    grok.shell.addTableView(df);
    const sketcherDialogs: DG.Dialog[] = [];

    DG.chem.currentSketcherType = 'OpenChemLib';
    const filter1 = await createFilter('Structure', df, sketcherDialogs);
    const filter2 = await createFilter('Structure', df, sketcherDialogs);

    //filter by structure and wait for results
    filter1.sketcher.setSmiles('C1CCCCC1');
    await waitUntil(() => df.filter.trueCount === 5, 'df hasn\'t been filtered 1', 5000);
    await waitUntil(() => filter2.bitset?.trueCount === 5, 'filters haven\'t been synchronized 1', 5000);
    await searchEnded(filter1);
    //filter1 is active filter, detaching filter2 and check that df is still filtered (detach refilters synchronously)
    filter2.detach();
    expect(df.filter.trueCount, 5, 'filter has been reset 1');

    const filter3 = await createFilter('Structure', df, sketcherDialogs);
    //filter by structure and wait for results
    filter1.sketcher.setSmiles('c1ccccc1');
    await waitUntil(() => df.filter.trueCount === 32, 'df hasn\'t been filtered 2', 5000);
    await waitUntil(() => filter3.bitset?.trueCount === 32, 'filters haven\'t been synchronized 2', 5000);
    await searchEnded(filter1);
    //filter1 is active filter, detaching active filter and check that df is still filtered
    filter1.detach();
    expect(df.filter.trueCount, 32, 'filter has been reset 2');

    const filter4 = await createFilter('Structure', df, sketcherDialogs);

    //detaching active filter3 while bitset hasn't yet been synchronized with filter4
    filter3.detach();
    expect(df.filter.trueCount, 32, 'filter has been reset 3');

    filter4.detach();
    sketcherDialogs.forEach((it) => it.close());
  });

  test('removedDisabledFilterReleasesPeers', async () => { //Remove All deactivates a card before it removes it
    const df = await readDataframe('tests/spgi-100.csv');
    await grok.data.detectSemanticTypes(df);
    grok.shell.addTableView(df);
    const sketcherDialogs: DG.Dialog[] = [];

    DG.chem.currentSketcherType = 'OpenChemLib';
    const filter1 = await createFilter('Structure', df, sketcherDialogs);
    const filter2 = await createFilter('Structure', df, sketcherDialogs);

    filter1.onActiveChanged(false);
    expect(filter2._peerFilterDisabled, true, 'the peer of a deactivated filter should stop applying');
    filter1.detach();
    expect(filter2._peerFilterDisabled, false, 'a removed filter should release its peers');

    filter2.sketcher.setSmiles('C1CCCCC1');
    await waitUntil(() => df.filter.trueCount === 5, 'the peer of a removed filter doesn\'t filter', 5000);
    await searchEnded(filter2);
    filter2.detach();
    sketcherDialogs.forEach((it) => it.close());
  });

  // a card restored from a layout gets its molecule from applyState in a view of its own: both cards keep the SMILES
  // the query was typed as (after a project was reopened, they named it by SMARTS)
  test('a filter card restored with a SMILES keeps it, and so does the card it syncs with', async () => {
    const df = await readDataframe('tests/sar-small_test.csv');
    await grok.data.detectSemanticTypes(df);
    const tv = grok.shell.addTableView(df);
    let tv2: DG.TableView | undefined;
    const card = async (view: DG.TableView): Promise<SubstructureFilter> => {
      const fg = view.getFiltersGroup({createDefaultFilters: false});
      let found: SubstructureFilter | undefined;
      await waitUntil(() => (found = fg.filters.find((f: any) => f?.columnName === 'smiles') as SubstructureFilter) !== undefined,
        'the substructure card was not added', 10000);
      return found!;
    };
    try {
      tv.getFiltersGroup({createDefaultFilters: false}).updateOrAdd(
        {type: DG.FILTER_TYPE.SUBSTRUCTURE, column: 'smiles', columnName: 'smiles', molBlock: 'c1ccncc1'}, false);
      const first = await card(tv);
      await waitUntil(() => !first.calculating && df.filter.trueCount < df.rowCount, 'the query did not filter', 15000);
      const filtered = df.filter.trueCount;
      const layout = tv.saveLayout();
      tv2 = grok.shell.addTableView(df);
      tv2.loadLayout(layout);
      const second = await card(tv2);
      await waitUntil(() => !second.calculating && second.currentMolecule !== '',
        'the restored card has no query', 15000);
      // the card takes its query from a timer applyState starts, kept in its subs (the only one that closes, once fired)
      await waitUntil(() => second.subs.some((s) => s.closed) && !second.calculating &&
        second.currentSearches.size === 0, 'the restored card did not take its query', 15000);
      expect(df.filter.trueCount, filtered);
      for (const [name, f] of [['restored', second], ['first', first]] as [string, SubstructureFilter][]) {
        expect(f.currentMolecule, 'c1ccncc1', `the ${name} card's molecule`);
        expect(f.filterSummary, 'c1ccncc1_Contains', `the ${name} card's summary`);
      }
    } finally {
      tv2?.close();
      tv.close();
    }
  });

  // the typed molecule kept for its own molblock must not swallow a query feature drawn on it, which leaves its SMILES
  test('a ring bond count drawn on a typed molecule is an edit', async () => {
    const df = await readDataframe('tests/sar-small_test.csv');
    await grok.data.detectSemanticTypes(df);
    const dialogs: DG.Dialog[] = [];
    const filter = await createFilter('smiles', df, dialogs);
    try {
      filter.sketcher.setSmiles('c1ccncc1');
      await waitUntil(() => filter.currentMolecule === 'c1ccncc1' && !filter.calculating,
        'the typed query did not filter', 10000);
      const mol = chemCommonRdKit.getRdKitModule().get_mol('c1ccncc1');
      const ringBondCount = mol.get_molblock().replace('M  END', 'M  RBC  1   1   2\nM  END');
      mol.delete();
      filter.sketcher.setMolFile(ringBondCount);
      await waitUntil(() => filter.currentMolecule.includes('M  RBC'),
        `the drawn ring bond count was taken for no edit: the filter keeps ${filter.currentMolecule}`, 10000);
    } finally {
      dialogs.forEach((it) => it.close());
      filter.detach();
    }
  });

  test('runMultipleUseASFilter', async () => { //#2628 (https://github.com/datagrok-ai/public/issues/2628)
    const df = await readDataframe('tests/sar-small_test.csv');
    await grok.data.detectSemanticTypes(df);
    const tv = grok.shell.addTableView(df);

    //filtering by 1st structure and wait for results
    tv.getFiltersGroup({ createDefaultFilters: false }).updateOrAdd({
      type: DG.FILTER_TYPE.SUBSTRUCTURE,
      column: 'smiles',
      columnName: 'smiles',
      molBlock: molblock1,
    }, false);
    await waitUntil(() => df.filter.trueCount === 4, 'df hasn\'t been filtered by 1st structure', 3000);

    //filtering by 2nd structure and wait for results
    tv.getFiltersGroup({ createDefaultFilters: false }).updateOrAdd({
      type: DG.FILTER_TYPE.SUBSTRUCTURE,
      column: 'smiles',
      columnName: 'smiles',
      molBlock: molblock2,
    }, false);
    await waitUntil(() => df.filter.trueCount === 3, 'df hasn\'t been filtered by 2nd structure', 15000);

    DG.chem.currentSketcherType = 'OpenChemLib';
  });

  // #2628: with a second view of the table whose Structure card was removed, Use as filter in the first view
  // (an updateOrAdd of its card) highlighted the value and filtered nothing
  test('Use as filter filters after a second view removed its Structure card', async () => {
    const df = await readDataframe('tests/spgi-100.csv');
    await grok.data.detectSemanticTypes(df);
    const tv = grok.shell.addTableView(df);
    let tv2: DG.TableView | undefined;
    const card = async (view: DG.TableView): Promise<SubstructureFilter> => {
      const fg = view.getFiltersGroup({createDefaultFilters: false});
      let found: SubstructureFilter | undefined;
      await waitUntil(() => (found = fg.filters.find((f: any) => f?.columnName === 'Structure') as SubstructureFilter) !== undefined,
        'the substructure card was not added', 10000);
      return found!;
    };
    const useAsFilter = (molBlock: string) => tv.getFiltersGroup({createDefaultFilters: false}).updateOrAdd(
      {type: DG.FILTER_TYPE.SUBSTRUCTURE, column: 'Structure', columnName: 'Structure', molBlock: molBlock}, false);
    try {
      useAsFilter('c1ccncc1');
      const first = await card(tv);
      await waitUntil(() => !first.calculating && df.filter.trueCount === 17,
        () => `pyridine: ${df.filter.trueCount} rows pass`, 15000);
      tv2 = grok.shell.addTableView(df);
      tv2.loadLayout(tv.saveLayout());
      const second = await card(tv2);
      tv2.getFiltersGroup({createDefaultFilters: false}).remove(second);
      grok.shell.v = tv;
      useAsFilter(df.col('Structure')!.get(0));
      await waitUntil(() => df.filter.trueCount === 1,
        () => `the molecule of row 1: ${df.filter.trueCount} rows pass`, 15000);
    }
    finally {
      tv2?.close();
      tv.close();
    }
  });

});


async function createFilter(colName: string, df: DG.DataFrame, sketcherDialogs: DG.Dialog[]):
  Promise<SubstructureFilter> {
  const filter = new SubstructureFilter();
  filter.attach(df);
  filter.applyState({columnName: colName});
  sketcherDialogs.push(ui.dialog().add(filter.root).show());
  await ui.tools.waitForElementInDom(filter.sketcher.root);
  filter.column = df.col(colName);
  filter.columnName = colName;
  filter.tableName = df.name;
  await filter.sketcher.sketcherReady();
  return filter;
}

async function testOneColumn(dfName: string, colName: string, substructure: string, terminateFlag: string,
  expectedKey: string, expectedTrueCount: number) {
  const df = await readDataframe(dfName);
  await grok.data.detectSemanticTypes(df);
  const sketcherDialogs: DG.Dialog[] = [];

  const filter = await createFilter(colName, df, sketcherDialogs);

  filter.sketcher.setSmiles(substructure);
  await waitUntil(() => df.filter.trueCount === expectedTrueCount, 'df hasn\'t been filtered', 10000);
  expectArray(df.filter.getBuffer(), expectedResults[expectedKey]);
  sketcherDialogs.forEach((it) => it.close());
  filter.detach();
}

/** Waits until the filter has no search under way: its loader off and its progress closed. */
async function searchEnded(filter: SubstructureFilter): Promise<void> {
  await waitUntil(() => !filter.calculating && filter.currentSearches.size === 0,
    () => `the search did not end (searches ${filter.currentSearches.size}, calculating ${filter.calculating})`, 15000);
}

