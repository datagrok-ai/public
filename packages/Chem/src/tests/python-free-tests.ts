/* eslint-disable max-len, new-cap -- PackageFunctions keep the scripts' names (Curate, Mutate, ...) */
/**
 * Chem's Python scripts and container conversions in Crux mode (python-free.ts; crux-core decision
 * chem-python-free): Curate, Mutate, Generate Conformers, the Gasteiger Partial Charges panel, Butina clustering,
 * the InChI converter, Calculate logP, Murcko Scaffolds, Filter by Catalogs and Map Identifiers' chemistry by Crux in
 * Crux mode, against the scripts' and the chem-chem container's own results on the pinned RDKit (crux-core
 * tools/chem-wasm/chem_scripts_oracle.py); by the script or the container in RDKit mode. The functions are this bundle's
 * (`PackageFunctions`), so that the engine set here is the one they run on.
 */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {waitUntil, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {ensureContainerRunning} from '@datagrok-libraries/test/src/test-container-utils';
import {_package} from '../package-test';
import {PackageFunctions} from '../package';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {ChemicalEngine} from '../engine/chem-engine';
import {loadCrux} from '../engine/crux-module';
import {getCruxService} from '../crux/crux-searches';
import {CurateOptions, PYTHON_SCRIPTS, butinaClusters, gasteigerPicture, inchiToInchiKey, inchiToSmiles, pythonScript,
  smilesToCanonical, smilesToInchi, smilesToInchiKey} from '../python-free';
import {getMapIdentifiers} from '../widgets/identifiers';
import {CONTAINER_TIMEOUT, readDataframe} from './utils';
import {expectNoFallback, inBothEngines, inEngine} from './engine-utils';

/** Chem's files/chem_standards.csv (the Curate sample) and what curate.py gives it with every step on (RDKit
 * 2026.03.6; Chem's `curate.smiles` test's values). */
const STANDARDS = ['CCC(=O)O[Na]', 'CCC(=O)[O-].[Na+]', '[Na]OC(=O)c1ccccc1', 'O=C([O-])c1ccccc1', 'C[N+](C)=CC=C[O-]',
  'CN(C)C=CC=O', 'C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O', 'O=S(O)c1ccc(S(=O)(=O)[O-])cc1', 'O=C([O-])c1ccccc1',
  'O=C(O)c1ccccc1', 'C1(=CCCCC1)O', 'O=C1CCCCC1', 'CCC1=C(C)C=CC(O)=N1.OC(=O)CCC(=O)O', 'CCC1=C(C)C=CC(O)=N1'];
const STANDARDS_CURATED = ['CCC(=O)[O-]', 'CCC(=O)O', 'O=C([O-])C1=CC=CC=C1', 'O=C(O)C1=CC=CC=C1', 'CN(C)C=CC=O',
  'CN(C)C=CC=O', 'O=S(O)C1=CC=C(S(=O)(=O)O)C=C1', 'O=S(O)C1=CC=C(S(=O)(=O)O)C=C1', 'O=C(O)C1=CC=CC=C1',
  'O=C(O)C1=CC=CC=C1', 'O=C1CCCCC1', 'O=C1CCCCC1', 'CCC1=C(C)C=CC(=O)N1', 'CCC1=C(C)C=CC(=O)N1'];
const EVERY_STEP = {kekulization: true, normalization: true, reionization: true, neutralization: true,
  tautomerization: true, mainFragment: true};
/** mutate.py's default molecule. */
const MUTATE_DEFAULT = 'CN1C(CC(O)C1=O)C1=CN=CC=C1';
/** butina_cluster.py's `cluster (Butina)` of files/tests/smiles-50.csv at the cutoffs 0.4 (the default) and 0.7. */
const BUTINA_50 = [4, 7, 4, 17, 8, 5, 16, 5, 5, 8, 7, 4, 15, 14, 13, 12, 3, 3, 3, 3, 11, 0, 0, 2, 2, 6, 1, 0, 0, 0, 4,
  0, 2, 0, 0, 6, 0, 0, 0, 0, 10, 0, 0, 0, 0, 9, 0, 0, 0, 0];
const BUTINA_50_07 = [1, 1, 1, 7, 2, 1, 6, 1, 1, 2, 3, 1, 12, 5, 5, 11, 2, 2, 2, 2, 10, 0, 0, 0, 0, 4, 0, 0, 0, 0, 1,
  0, 0, 0, 0, 4, 0, 0, 0, 0, 9, 0, 0, 0, 0, 8, 0, 0, 0, 0];

/** Molecules for Calculate logP, Murcko Scaffolds and Filter by Catalogs, and what each script gives them on the pinned RDKit 2026.03.6 (the
 * oracle): a SMILES, an acyclic one, a blank cell, a broken SMILES, a stereo ring, a quinone, a catechol, an azo dye, a
 * pyrrole without its H (MolFromSmiles reads none), a MOL block (MolFromSmiles reads none; calc_logP.py reads it). */
const MOLECULES = ['Cc1cc(Oc2nccc(CCC)c2)ccc1', 'CCO', '', 'C1CC(', 'C[C@H]1CC[C@@H](C)N1Cc1ccccc1', 'O=C1C=CC(=O)C=C1',
  'Oc1ccccc1O', 'CN(C)c1ccc(N=Nc2ccccc2)cc1', 'c1cccn1',
  '\n     RDKit          2D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2990    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\nM  END\n'];
/** calc_logP.py's `clogP` (null: None). */
const LOGP = [4.134820000000003, -0.0014000000000000123, null, null, 3.059400000000002, 0.25059999999999993,
  1.0977999999999999, 4.168000000000003, null, -0.3915];
/** murcko_scaffold.py's `scaffolds` ('': the script's None as well as its ''). */
const SCAFFOLDS = ['c1ccc(Oc2ccccn2)cc1', '', '', '', 'c1ccc(CN2CCCC2)cc1', 'O=C1C=CC(=O)C=C1', 'c1ccccc1',
  'c1ccc(N=Nc2ccccc2)cc1', '', ''];
/** filter_catalogs.py's flags. */
const PAINS_A = [false, false, false, false, false, true, false, true, false, false];
const ZINC = [false, false, false, false, false, true, false, false, false, false];
const ALANINE = 'C[C@H](N)C(=O)O';
const ALANINE_INCHI = 'InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1';

/**
 * A function's output column `name` as its `join(<table>)` output leaves it: the script (RDKit mode) joins its output
 * into the input table itself, called without its outputs processed too, and returns its output table emptied; Crux
 * returns the table with the column, which the platform joins.
 */
const outputColumn = (result: DG.DataFrame, table: DG.DataFrame, name: string): DG.Column | null =>
  table.col(name) ?? result.col(name);

/** Curate's `curated_molecule` of the table's column (the options not given: the function's defaults). */
async function curated(df: DG.DataFrame, column: string, options: Partial<CurateOptions>): Promise<string[]> {
  const o = options;
  const result = await PackageFunctions.Curate(df, df.col(column)!, o.kekulization, o.normalization, o.reionization,
    o.neutralization, o.tautomerization, o.mainFragment);
  const col = outputColumn(result, df, 'curated_molecule');
  expect(col !== null, true, 'Curate gives curated_molecule');
  return col!.toList();
}

/** `df` with the molecules in a column `smiles`, its semantic type detected. */
async function moleculesTable(values: string[]): Promise<DG.DataFrame> {
  const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('smiles', values)]);
  await grok.data.detectSemanticTypes(df);
  return df;
}

/** Equal within a float column's precision (Datagrok keeps a float column's values in 32 bits). */
const close = (a: number | null, b: number | null) => a === null || b === null ? a === b :
  Math.abs(a - b) <= 1e-6 * Math.max(1, Math.abs(b));

/** A column's values, null where empty (`toList()` leaves holes there). */
const cells = (col: DG.Column): any[] => Array.from({length: col.length}, (_, i) => col.isNone(i) ? null : col.get(i));

/** The width of a panel picture's image (its background), 0 where it does not load. */
function pictureWidth(picture: HTMLElement): Promise<number> {
  const url = /url\(["']?(.*?)["']?\)/.exec(picture.style.backgroundImage)?.[1] ?? '';
  return new Promise((resolve) => {
    const img = new Image();
    img.onload = () => resolve(img.naturalWidth);
    img.onerror = () => resolve(0);
    img.src = url;
  });
}

category('chemical engine: python scripts', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('Callers find them: the package functions under the scripts\' names, the scripts under theirs', async () => {
    // Flow lists and loads functions by name, case-sensitively (a saved flow's node: `DG Functions/<role>/<name>`), and
    // its allowlist names them so (`Chem:CalculateLogP`, ...)
    for (const name of ['Curate', 'Mutate', 'GenerateConformers', 'ButinaMoleculesClustering', 'InchiToMol',
      'CalculateLogP', 'MurckoScaffolds', 'FilterByCatalogs']) {
      const f = DG.Func.find({package: 'Chem', name}).find((g) => g.nqName === `Chem:${name}`);
      expect(!!f && !(f instanceof DG.Script), true, `the package function Chem:${name}`);
    }
    for (const name of Object.values(PYTHON_SCRIPTS)) {
      const script = pythonScript(name);
      expect(script instanceof DG.Script && script.language === 'python', true, `the script ${name}`);
    }
    // the Gasteiger script stays a panel (Flow's node type) that no molecule shows
    const gasteiger = pythonScript(PYTHON_SCRIPTS.gasteiger);
    expect(gasteiger.options['role'], 'panel');
    expect(gasteiger.options['condition'], 'false');
    // the Chemical Properties dialog lists its calculators by this tag and shows `logP` (the friendly name less
    // "Calculate "): once
    const logPs = DG.Func.find({meta: {function_family: 'biochem-calculator'}}).filter((f) => f.friendlyName === 'Calculate logP');
    expect(logPs.map((f) => f.nqName).join(), 'Chem:CalculateLogP');
  });

  test('Gasteiger Partial Charges: one pane, the package\'s panel', async () => {
    // the script is a panel for Flow's sake that no molecule shows (`#condition: false`)
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('smiles', ['CC(=O)O'])]);
    await grok.data.detectSemanticTypes(df);
    const tv = grok.shell.addTableView(df);
    grok.shell.windows.showContextPanel = true;
    try {
      df.currentCell = df.cell(0, 'smiles');
      const headers = () => Array.from(document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header')) as HTMLElement[];
      const chemistry = () => headers().find((h) => h.textContent === 'Chemistry');
      await waitUntil(() => !!chemistry(), 'no Chemistry pane', 10000);
      if (!chemistry()!.classList.contains('expanded'))
        chemistry()!.click();
      const panes = () => headers().filter((h) => (h.textContent ?? '').startsWith('Gasteiger Partial Charges')).length;
      await waitUntil(() => panes() > 0, 'no Gasteiger Partial Charges pane', 10000).catch(() => {
        throw new Error(`no Gasteiger Partial Charges pane among ${headers().map((h) => h.textContent).join(', ')}`);
      });
      expect(panes(), 1);
    } finally {
      tv.close();
    }
  }, {timeout: 30000});

  test('Curate: Crux gives the script\'s SMILES, every step on', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('smiles', STANDARDS)]);
    const values = await inEngine(ChemicalEngine.Crux, () =>
      expectNoFallback(['Curate'], () => curated(df, 'smiles', EVERY_STEP)));
    expect(JSON.stringify(values), JSON.stringify(STANDARDS_CURATED));
  }, {timeout: 60000});

  test('Curate: empty and malformed molecules, both engines', async () => {
    // RDKit mode: the script (the stand's Python RDKit); Crux mode: Crux; '' where the molecule does not read
    const values = await inBothEngines(async (engine) => {
      const df = await readDataframe('tests/Test_smiles_malformed.csv');
      await grok.data.detectSemanticTypes(df);
      const curatedValues = await curated(df, 'canonical_smiles', {});
      expect(curatedValues.length, df.rowCount);
      expect(curatedValues.filter((v) => v).length, 43, `${engine}: the readable molecules`);
      return curatedValues;
    }, {noFallback: ['Curate']});
    // the script on the stand's RDKit and Crux (RDKit 2026.03) agree on these
    expect(JSON.stringify(values.crux), JSON.stringify(values.rdkit));
  }, {timeout: 180000});

  test('Mutate: Crux draws the script\'s mutations under a seed, and makes every one', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      // mutate.py after random.seed(42) (CPython's Mersenne Twister)
      const seeded = await (await getCruxService()).idleWorkerJob<string[]>('chemMutate',
        [MUTATE_DEFAULT, JSON.stringify({seed: 42, maxRandomResults: 8})], 30000);
      expect(JSON.stringify(seeded), JSON.stringify(['CNC(=O)C(O)CCC1=CN=CC=C1', 'CN1C(=O)CCC1C1=CN=CC=C1',
        'CN1C(=O)C(O)CC1(N)C1=CN=CC=C1', 'C=CC(=CN=C)C1CC(O)C(=O)N1C', 'CN1CC(O)CC1C1=CN=CC=C1',
        'CN1C(=O)C(C)(O)CC1C1=CN=CC=C1', 'O=C1C(O)CC(C2=CN=CC=C2)N1CCl', 'CC(C1=CN=CC=C1)N(C)C(=O)CO']));
      await expectNoFallback(['Mutate'], async () => {
        const all = await PackageFunctions.Mutate(MUTATE_DEFAULT, 1, false, 100);
        expect(all.rowCount, 95, 'the script\'s 95 single mutations');
        expect(all.col('mutations')!.semType, DG.SEMTYPE.MOLECULE);
        expect((await PackageFunctions.Mutate(MUTATE_DEFAULT, 2, true, 20)).rowCount, 20);
      });
    });
  }, {timeout: 60000});

  test('Generate Conformers: Crux gives the script\'s table', async () => {
    const df = await inEngine(ChemicalEngine.Crux, () => expectNoFallback(['Generate Conformers'], () =>
      PackageFunctions.GenerateConformers('CCCC', 50, true, 0.1, 5000, 42)));
    expect(df.columns.names().join(','), 'smiles,molblock,conformer,energy,rmsd');
    // generate_conformers.py on the pinned RDKit: 4 conformers, MMFF94 energies -4.2938, -5.076, -4.2938, -4.2938
    expect(df.rowCount, 4);
    const energies = df.col('energy')!.toList();
    [-4.2938, -5.076, -4.2938, -4.2938].forEach((e, i) => expect(Math.abs(energies[i] - e) < 1e-3, true,
      `energy ${i + 1}: ${energies[i]}`));
    expect(df.col('rmsd')!.get(0), 0);
    expect(df.col('smiles')!.get(0), 'CCCC');
    expect((df.col('molblock')!.get(0) as string).includes(' H '), true, 'the conformers keep their hydrogens');
  }, {timeout: 60000});

  test('Gasteiger Partial Charges: Crux draws the charges, the script in RDKit mode', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      const crux = await loadCrux('tests');
      // RDKit 2026.03.6's `_GasteigerCharge` of acetic acid
      const map = JSON.parse(crux.chemGasteigerMap('CC(=O)O', 10, 400, 400));
      [0.03376759340677524, 0.2996845531582142, -0.2528204334965765, -0.481433475172105].forEach((c, i) =>
        expect(Math.abs(map.charges[i] - c) < 1e-12, true, `charge ${i}: ${map.charges[i]}`));
      expect(map.svg.includes('width=\'400px\''), true, 'the script\'s 400 x 400');
      const picture = await expectNoFallback(['Gasteiger Partial Charges'], () => gasteigerPicture('CC(=O)O', 10));
      expect(picture.style.backgroundImage.includes('image/png'), true, 'Crux\'s picture, as a PNG');
      expect(await pictureWidth(picture), 400, 'Crux\'s picture loads');
      expect(JSON.parse(crux.chemGasteigerMap('not a molecule', 10)), null, 'no picture, as the script draws none');
    });
    // RDKit mode (and Crux's fallback): the script `Chemistry | Gasteiger Partial Charges` draws it, a PNG
    // (reported: the platform hands the script's graphics over as an element, which was read as the PNG's base64)
    const png = await inEngine(ChemicalEngine.RDKit, () => gasteigerPicture('CC(=O)O', 10));
    expect(await pictureWidth(png), 400, `the script's picture loads (${png.style.backgroundImage.slice(0, 60)})`);
  }, {timeout: 120000});

  test('inchiToMol: the converter gives the script\'s molecule, both engines', async () => {
    // RDKit mode: the script; Crux mode: the Crux sketcher's InChI module (RDKit's MolFromInchi), no script
    await inBothEngines(async (engine) => {
      expect(await PackageFunctions.InchiToMol('InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H'), 'c1ccccc1', engine);
      expect(await PackageFunctions.InchiToMol(ALANINE_INCHI), ALANINE, engine);
    }, {noFallback: ['inchiToMol']});
  }, {timeout: 120000});

  test('Butina: Crux gives the script\'s clusters', async () => {
    await inEngine(ChemicalEngine.Crux, () => expectNoFallback(['Butina Molecules Clustering'], async () => {
      for (const [cutoff, want] of [[0.4, BUTINA_50], [0.7, BUTINA_50_07]] as const) {
        const df = await readDataframe('tests/smiles-50.csv');
        const clusters = await PackageFunctions.ButinaMoleculesClustering(df, df.col('canonical_smiles')!, cutoff);
        expect(JSON.stringify(clusters.col('cluster (Butina)')!.toList()), JSON.stringify(want), `cutoff ${cutoff}`);
      }
    }));
    // RDKit's ClusterData on two points within the cutoff and a lone one: the pair first (the later point its
    // centroid), then the singleton
    expect(Array.from(butinaClusters(Int32Array.from([2, 2, 1]), (i) => [[0, 1], [0, 1], [2]][i])).join(), '0,0,1');
  }, {timeout: 60000});

  test('Calculate logP: Crux gives the script\'s values, the script in RDKit mode, one column each', async () => {
    const values = await inBothEngines(async (engine) => {
      const df = await moleculesTable(MOLECULES);
      const result = await PackageFunctions.CalculateLogP(df, df.col('smiles')!);
      expect([...df.columns.names(), ...result.columns.names()].filter((n) => n.startsWith('clogP')).join(), 'clogP',
        `${engine}: one clogP column`);
      return cells(outputColumn(result, df, 'clogP')!);
    }, {noFallback: ['Calculate logP']});
    LOGP.forEach((v, i) => expect(close(values.crux[i], v), true, `Crux, row ${i}: ${values.crux[i]}, the script ${v}`));
    // the stand's script (its kernel's RDKit) and Crux (RDKit 2026.03) agree on these
    values.rdkit.forEach((v, i) => expect(close(v, values.crux[i]), true,
      `RDKit mode, row ${i}: ${v}, Crux ${values.crux[i]}`));
  }, {timeout: 120000});

  test('Murcko Scaffolds and Filter by Catalogs: Crux gives the scripts\' columns, both engines', async () => {
    await inBothEngines(async (engine) => {
      const df = await moleculesTable(MOLECULES);
      const scaffolds = await PackageFunctions.MurckoScaffolds(df, df.col('smiles')!);
      expect(JSON.stringify(cells(outputColumn(scaffolds, df, 'scaffolds')!).map((s) => s ?? '')),
        JSON.stringify(SCAFFOLDS), `${engine}: scaffolds`);
      const zinc = await PackageFunctions.FilterByCatalogs(df, df.col('smiles')!, 'ZINC');
      expect(JSON.stringify(cells(outputColumn(zinc, df, 'ZINC')!)), JSON.stringify(ZINC), `${engine}: ZINC`);
    }, {noFallback: ['Murcko Scaffolds', 'Filter by Catalogs']});
    const df = await moleculesTable(MOLECULES);
    const pains = await inEngine(ChemicalEngine.Crux, () => PackageFunctions.FilterByCatalogs(df, df.col('smiles')!,
      'PAINS_A'));
    expect(JSON.stringify(cells(pains.col('PAINS_A')!)), JSON.stringify(PAINS_A), 'Crux: PAINS_A');
  }, {timeout: 120000});

  test('Map Identifiers: Crux gives the container\'s SMILES, InChIs and keys, the container in RDKit mode', async () => {
    await ensureContainerRunning('name = "chem-chem"', CONTAINER_TIMEOUT);
    const smiles = ['OCC', 'C1CC(', '', ALANINE, 'F/C=C/F', 'C[C@H]1CC[C@@H](C)N1Cc1ccccc1'];
    // the InChI library hashes a malformed InChI's layers into a key all the same; an InChIKey is no InChI
    const inchis = [ALANINE_INCHI, '', 'InChI=1S/garbage', 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+',
      'QNAYBMKLOCPYGJ-REOHCLBHSA-N'];
    const got = await inBothEngines(async () => JSON.stringify([await smilesToCanonical(smiles),
      await smilesToInchi(smiles), await smilesToInchiKey(smiles), await inchiToSmiles(inchis),
      await inchiToInchiKey(inchis)]), {noFallback: ['Map Identifiers']});
    // the container's chem.py on the pinned RDKit (inchi_to_inchi_key: None for '' and an InChI the library refuses)
    expect(got.crux, JSON.stringify([
      ['CCO', '', '', ALANINE, 'F/C=C/F', 'C[C@@H]1CC[C@H](C)N1Cc1ccccc1'],
      ['InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', '', '', ALANINE_INCHI, 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+',
        'InChI=1S/C13H19N/c1-11-8-9-12(2)14(11)10-13-6-4-3-5-7-13/h3-7,11-12H,8-10H2,1-2H3/t11-,12+'],
      ['LFQSCWFLJHTTHZ-UHFFFAOYSA-N', '', '', 'QNAYBMKLOCPYGJ-REOHCLBHSA-N', 'WFLOTYSKFUPZQB-OWOJBTEDSA-N',
        'KKRNGUFZPFNIFS-TXEJJXNPSA-N'],
      [ALANINE, '', '', 'F/C=C/F', ''],
      ['QNAYBMKLOCPYGJ-REOHCLBHSA-N', null, 'LKNHRTZJKNITMJ-UHFFFAOYSA-N', 'WFLOTYSKFUPZQB-OWOJBTEDSA-N', null],
    ]));
    expect(got.rdkit, got.crux, 'the container (its RDKit) and Crux agree on these');
    // the menu's path (widgets/identifiers.ts) on Crux
    await inEngine(ChemicalEngine.Crux, () => expectNoFallback(['Map Identifiers'], async () => {
      const df = await moleculesTable(smiles);
      await getMapIdentifiers(df, df.col('smiles')!, 'smiles', 'inchi_key');
      expect(df.col('inchi_key')?.get(3), 'QNAYBMKLOCPYGJ-REOHCLBHSA-N');
    }));
  }, {timeout: 60000 + CONTAINER_TIMEOUT});
});
