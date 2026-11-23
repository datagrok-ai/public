/* eslint-disable max-len */
/**
 * What moved to Crux in Crux mode after the first features (src/engine/README.md, "What follows the setting"): the column
 * jobs of engine/molecule-jobs.ts on the Crux workers (CoordGen coordinates, flattening, beautifying, elemental
 * analysis, the R-group dialog's cores), the reactions' canonical SMILES, `getCLogP`, the atom picker's 2D→3D mapping;
 * drawing, properties and QED of a CXSMILES with query features; the fallback report; `Chem:removeDuplicates`. Each test
 * runs both engines where the feature has both, and shows the Crux path ran by a molecule marked unreadable for Crux
 * (`withCruxUnreadable`), which RDKit then does, counted as a fallback.
 */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {PackageFunctions} from '../package';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {ChemicalEngine} from '../engine/chem-engine';
import {cruxModuleSync} from '../engine/crux-module';
import {cruxFallbackCount, cruxFallbackReason, quietFallbacksForTests, setFallbackReporterForTests, warnCruxFallback}
  from '../engine/fallback';
import {atomCountsByEngine, atomsColumnByEngine, beautifyMolsByEngine, coordGenCoordsByEngine, coreMolblocksByCrux,
  flattenMoleculesByEngine, validationTextsByEngine} from '../engine/molecule-jobs';
import {validationText} from '../utils/mol-creation_rdkit';
import {MAX_SMILES_LENGTH} from '../utils/chem-constants';
import {canonicalSmiles, standardizeFragment} from '../analysis/sar-matrix/build/sar-matrix-columns';
import {cachedAtomCount} from '../analysis/sar-matrix/build/sar-matrix-decompose';
import {getAtomsColumn} from '../utils/elemental-analysis-utils';
import {removeWaterAndSalts} from '../utils/reactions/reactions';
import {mapAtomIndices2Dto3DByEngine} from '../utils/atom-index-mapper';
import {RDKitCellRenderer} from '../rendering/rdkit-cell-renderer';
import {getPropertiesAsColumns} from '../widgets/properties';
import {cruxQed} from '../widgets/drug-likeness';
import {readDataframe} from './utils';
import {expectNoFallback, inBothEngines, inEngine, withCruxUnreadable} from './engine-utils';
import {_package} from '../package-test';

/** Element counts as comparable text: each element with an atom, its counts; the invalid rows. */
function countsText([elements, invalid]: [Map<string, Int32Array>, number[]]): string {
  const used = [...elements.entries()].filter(([, counts]) => counts.some((n) => n !== 0))
    .sort(([a], [b]) => a < b ? -1 : a > b ? 1 : 0);
  return used.map(([symbol, counts]) => `${symbol}: ${counts.join(',')}`).join('\n') +
    `\ninvalid: ${invalid.filter((i) => i !== undefined).join(',')}`;
}

const SALT = 'c1ccccc1[NH3+].[Cl-]';
const UNREADABLE_BY_CRUX = 'CC(=O)Nc1ccc(O)cc1'; // marked unreadable for the Crux workers in the tests below

category('chemical engine: column jobs', () => {
  let rdkit: any;

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    rdkit = chemCommonRdKit.getRdKitModule();
  });

  test('elemental analysis: Crux mode counts as RDKit mode on SMILES, V2000 and V3000', async () => {
    const spgi = await readDataframe('tests/spgi-100.csv');
    const drugs = await readDataframe('tests/approved-drugs-100.csv');
    const smiles = DG.Column.fromStrings('smiles', [...grok.data.demo.molecules(100).col('smiles')!.toList(), '',
      'C1CC(', UNREADABLE_BY_CRUX, 'CC(=O)[O-].[Na+]']);
    smiles.meta.units = DG.UNITS.Molecule.SMILES;
    for (const col of [smiles, spgi.col('Structure')!, drugs.col('molecule')!]) {
      const {rdkit, crux} = await inBothEngines(() => atomsColumnByEngine(col, getAtomsColumn),
        {noFallback: ['elemental analysis']});
      expect(countsText(crux), countsText(rdkit), `${col.name}: Crux's counts`);
    }
    // a molecule Crux does not read: RDKit counts it, the fallback counted
    await withCruxUnreadable([UNREADABLE_BY_CRUX], async () => {
      const rdkitMode = await inEngine(ChemicalEngine.RDKit, () => atomsColumnByEngine(smiles, getAtomsColumn));
      const fallbacks = cruxFallbackCount('elemental analysis');
      const cruxMode = await inEngine(ChemicalEngine.Crux, () => atomsColumnByEngine(smiles, getAtomsColumn));
      expect(countsText(cruxMode), countsText(rdkitMode), 'with the fallback');
      expect(cruxFallbackCount('elemental analysis'), fallbacks + 1, 'RDKit counted the molecule Crux did not read');
    });
  });

  test('CoordGen coordinates, flattening and beautifying on the Crux workers', async () => {
    const molecules = [SALT, 'C[C@H](N)C(=O)O', 'C/C=C/C', '', 'C1CC(', UNREADABLE_BY_CRUX];
    const jobs: [string, (m: string[]) => Promise<string[]>, (mol: any) => string][] = [
      ['CoordGen coordinates', coordGenCoordsByEngine, (mol) => mol.get_new_coords(true)],
      ['flattenMolecules', (m) => flattenMoleculesByEngine(m, true),
        (mol) => mol.get_smiles(JSON.stringify({doIsomericSmiles: false}))],
      ['beautifyMols', beautifyMolsByEngine, (mol) => {
        mol.set_new_coords();
        mol.normalize_depiction(1);
        mol.straighten_depiction(true);
        return mol.get_v3Kmolblock();
      }],
    ];
    const by = (module: any, molecule: string, write: (mol: any) => string): string => {
      const mol = module.get_mol(molecule);
      try {
        return write(mol);
      } finally {
        mol.delete();
      }
    };
    await withCruxUnreadable([UNREADABLE_BY_CRUX], async () => {
      for (const [call, job, write] of jobs) {
        const rdkitMode = await inEngine(ChemicalEngine.RDKit, () => job(molecules));
        const fallbacks = cruxFallbackCount(call);
        const cruxMode = await inEngine(ChemicalEngine.Crux, () => job(molecules));
        // Crux's answers are the main thread's Crux module's (RDKit 2026.03.6's): CoordGen's MOL block of the salt
        // has no valence fields where Chem's RDKit 2024.09 writes N+'s and Cl-'s, and another Kekulé form
        for (let i = 0; i < 3; i++) {
          expect(cruxMode[i], by(cruxModuleSync(), molecules[i], write), `${call}: Crux's for ${molecules[i]}`);
          expect(rdkitMode[i], by(rdkit, molecules[i], write), `${call}: RDKit's for ${molecules[i]}`);
        }
        expect(cruxMode[3], rdkitMode[3], `${call}: an empty cell`);
        expect(cruxMode[4], '', `${call}: a molecule neither reads`);
        expect(cruxMode[5], rdkitMode[5], `${call}: RDKit did the molecule Crux did not read`);
        expect(cruxFallbackCount(call), fallbacks + 1, `${call}: one fallback`);
      }
    });
  });

  test('flattening keeps E/Z unless asked, in both engines', async () => {
    const molecules = ['C/C=C/[C@H](N)C', '[2H][C@](F)(Cl)Br', 'C1CC/C=C/CCC1', 'C[C@H](N)C(=O)O |&1:1|',
      'C[S@@](=O)c1ccccc1', 'C/C=C(/F)[C@@H]1CC[C@H](C)CC1', '', 'C1CC('];
    const expected: [boolean, string[]][] = [
      [false, ['C/C=C/C(C)N', '[2H]C(F)(Cl)Br', 'C1=C/CCCCCC/1', 'CC(N)C(=O)O', 'CS(=O)c1ccccc1',
        'C/C=C(/F)C1CCC(C)CC1', '', '']],
      [true, ['CC=CC(C)N', '[2H]C(F)(Cl)Br', 'C1=CCCCCCC1', 'CC(N)C(=O)O', 'CS(=O)c1ccccc1', 'CC=C(F)C1CCC(C)CC1',
        '', '']],
    ];
    for (const [flattenDoubleBonds, smiles] of expected) {
      const {rdkit, crux} = await inBothEngines(() => flattenMoleculesByEngine(molecules, flattenDoubleBonds),
        {noFallback: ['flattenMolecules']});
      expectArray(rdkit, smiles);
      expectArray(crux, smiles);
    }
  });

  test('getCLogP: RDKit 2026.03\'s Crippen logP in Crux mode', async () => {
    const {rdkit, crux} = await inBothEngines(async () => PackageFunctions.getCLogP('C[C@H](N)C(=O)O'));
    // the same Crippen sum; Chem's RDKit.js 2024.09 writes get_descriptors' numbers cut to 5 decimals, 2026.03 in full
    expect(rdkit, -0.58179);
    expect(Math.abs(crux + 0.5818) < 1e-12, true, `Crux mode: ${crux}`);
  });

  test('reactions\' canonical SMILES and Remove Water and Salts on the Crux workers', async () => {
    const molecules = ['OC(=O)c1ccccc1.O', 'Cl.CCO', UNREADABLE_BY_CRUX + '.Cl', 'C1CC('];
    await withCruxUnreadable([UNREADABLE_BY_CRUX + '.Cl'], async () => {
      const rdkitMode = await inEngine(ChemicalEngine.RDKit, () => removeWaterAndSalts(molecules));
      const fallbacks = cruxFallbackCount('convertMolNotation');
      const cruxMode = await inEngine(ChemicalEngine.Crux, () => removeWaterAndSalts(molecules));
      expect(cruxMode.join(' '), rdkitMode.join(' '));
      expect(cruxMode.slice(0, 3).join(' '), 'O=C(O)c1ccccc1 CCO CC(=O)Nc1ccc(O)cc1');
      expect(cruxFallbackCount('convertMolNotation'), fallbacks + 1, 'RDKit did the molecule Crux did not read');
    });
  });

  test('the R-group dialog\'s cores from the Crux workers', async () => {
    const cores = ['c1ccc([*:1])cc1', '', 'C'.repeat(6000), UNREADABLE_BY_CRUX];
    expect(await inEngine(ChemicalEngine.RDKit, () => coreMolblocksByCrux(cores)), null, 'RDKit mode: the dialog\'s code');
    await withCruxUnreadable([UNREADABLE_BY_CRUX], async () => {
      const blocks = await inEngine(ChemicalEngine.Crux, () => coreMolblocksByCrux(cores));
      const mol = cruxModuleSync()!.get_mol(cores[0])!;
      const expected = mol.get_molblock().replace('ISO', 'RGP');
      mol.delete();
      expect(blocks![0], expected, 'Crux\'s MOL block, its ISO line written as RGP');
      expect(blocks!.slice(1).map((b) => b === undefined ? 'undefined' : String(b)).join(','), 'null,null,undefined',
        'an empty core, an over-long SMILES (none), a core the dialog\'s RDKit code does');
    });
  });

  test('the SAR matrix\'s keys and atom counts by the engine', async () => {
    // a molecule Chem's RDKit 2024.09 writes with other stereo marks than RDKit 2026.03 (Crux): the SAR matrix compares
    // its keys with the workers' canonical SMILES (Crux's in Crux mode), so its own are Crux's there too
    const styrene = 'C(=C\\C1CCCN(CCc2ccccc2)C1)\\c1ccccc1';
    const {rdkit: rdkitMode, crux: cruxMode} = await inBothEngines(async () => [canonicalSmiles(styrene),
      standardizeFragment(styrene + '.[*:1]C').value]);
    expect(rdkitMode[0], 'C(=C\\C1CCCN(CCc2ccccc2)C1)\\c1ccccc1', 'RDKit 2024.09\'s');
    expect(cruxMode[0], 'C(=C/C1CCCN(CCc2ccccc2)C1)/c1ccccc1', 'RDKit 2026.03\'s, Crux\'s');
    expect(cruxMode[1].includes('C(=C/C1CCCN(CCc2ccccc2)C1)/c1ccccc1'), true, `Crux's fragment: ${cruxMode[1]}`);
    // the workers' counts, the main thread's for a structure they do not read
    const structures = ['c1ccccc1', 'CC(=O)O', '', 'C1CC(', UNREADABLE_BY_CRUX];
    await withCruxUnreadable([UNREADABLE_BY_CRUX], async () => {
      const fallbacks = cruxFallbackCount('atom counts');
      const counts = await inEngine(ChemicalEngine.Crux, () => atomCountsByEngine(structures,
        (s) => cachedAtomCount(s)));
      expect(counts.join(','), '6,4,0,0,11');
      expect(cruxFallbackCount('atom counts'), fallbacks + 1, 'the main thread counted the one the workers did not');
    });
  });

  test('the atom picker\'s 2D→3D mapping by the engine', async () => {
    const pdb = [
      'ATOM      1  O   LIG     1       2.460   0.000   0.000  1.00  0.00           O  ',
      'ATOM      2  C2  LIG     1       1.230   0.000   0.000  1.00  0.00           C  ',
      'ATOM      3  C1  LIG     1       0.000   0.000   0.000  1.00  0.00           C  ',
    ].join('\n');
    const mapped = await inBothEngines(async () => JSON.stringify(mapAtomIndices2Dto3DByEngine(rdkit, 'CCO', pdb)),
      {noFallback: ['atom index mapping']});
    expect(mapped.crux, mapped.rdkit);
  });

  test('a CXSMILES query row: Crux draws it and gives its properties and QED, as RDKit', async () => {
    // Crux reads every field of a CXSMILES extension as RDKit 2026.03 (crux-core e881b8e98): its query features
    // (rb:, m:, ...) no longer send the row to RDKit's drawing or leave Crux's properties and QED empty
    const queries = ['C1CCCCC1 |rb:0:2|', '*c1ccccc1 |m:0:1.2.3|'];
    await inEngine(ChemicalEngine.Crux, () => expectNoFallback(['rendering', 'properties', 'drugLikeness'], async () => {
      const renderer = new RDKitCellRenderer(rdkit) as any;
      for (const query of queries) {
        const drawn = renderer._fetchMol(query, [], false, false, {}, false);
        expect(drawn.engine, ChemicalEngine.Crux, `drawn: ${query}`);
        expect(drawn.molCtx.mol !== null, true, `read: ${query}`);
      }
      expect(renderer._fetchMol('C1CCCCC1', [], false, false, {}, false).engine, ChemicalEngine.Crux, 'a plain row');
      expect(renderer._fetchMol('c1ccccc1CC', [{molecule: queries[1]}], false, false, {}, false).engine,
        ChemicalEngine.Crux, 'a CXSMILES query scaffold');
      // the query features leave the molecule as it is: MW and QED those of the molecule without them
      const col = DG.Column.fromStrings('smiles', ['CCO', ...queries, 'C1CCCCC1', '*c1ccccc1']);
      const [mw] = await getPropertiesAsColumns(col, ['MW']);
      expect(mw.isNone(0), false, 'CCO\'s MW');
      expect([1, 2].map((i) => mw.isNone(i)).join(','), 'false,false', 'Crux\'s MW for a query row');
      expect(`${mw.get(1)} ${mw.get(2)}`, `${mw.get(3)} ${mw.get(4)}`, 'the molecules\' MW');
      expect(JSON.stringify(cruxQed(queries[0])), JSON.stringify(cruxQed('C1CCCCC1')), 'Crux\'s QED of a query row');
      expect(cruxQed(queries[0]) !== null, true, 'a QED');
      expect(cruxQed('CCO') !== null, true, 'Crux\'s QED of CCO');
    }));
  });

  test('Chem:removeDuplicates, as js-api\'s Sketcher calls it; the cell editor\'s "edited?" check', async () => {
    const molecules = ['OCC', 'c1ccccc1', 'C(O)C', 'C1CC('];
    // the published package's function, in the stand's engine (setChemicalEngine reaches this bundle's code only)
    const kept = await grok.functions.call('Chem:removeDuplicates', {molecules, molecule: 'CCO'}) as string[];
    expect(kept.join(' '), 'c1ccccc1 C1CC(', 'Chem:removeDuplicates');
    await inBothEngines(async (engine) => {
      expect(PackageFunctions.removeDuplicates(molecules, 'CCO').join(' '), 'c1ccccc1 C1CC(', `${engine} mode`);
      expect(chemCommonRdKit.sameMolecule('OCC', 'CCO'), true, `${engine} mode: one molecule`);
      // the OK of a cell whose molecule does not read did nothing: the check threw
      expect(chemCommonRdKit.sameMolecule('C1CC(', 'CCO'), false, `${engine} mode: a cell that does not read is edited`);
    }, {noFallback: ['removeDuplicates']});
  });

  test('the first fallback of a call is reported to the platform\'s log, once', async () => {
    const reported: string[] = [];
    const previous = setFallbackReporterForTests((call, message) => reported.push(`${call}: ${message}`));
    try {
      const call = `test report ${Date.now()}`;
      warnCruxFallback(call, 'first');
      warnCruxFallback(call, new Error('second'));
      expect(reported.join(' | '), `${call}: first`);
      expect(cruxFallbackCount(call), 2, 'both counted');
    } finally {
      setFallbackReporterForTests(previous);
    }
  });

  test('a fallback a test makes on purpose is counted, not reported, and the next real one is still reported',
    async () => {
      // the tests' own fallbacks (a molecule marked unreadable) reached the platform's log as real use's, and took the
      // call's one report of the session: a real fallback of that call after them was never logged
      const reported: string[] = [];
      const previous = setFallbackReporterForTests((call, message) => reported.push(`${call}: ${message}`));
      try {
        const call = `test quiet ${Date.now()}`;
        await quietFallbacksForTests(async () => warnCruxFallback(call, 'on purpose'));
        expect(reported.join(' | '), '', 'the test\'s own fallback: not reported');
        expect(cruxFallbackCount(call), 1, 'counted');
        expect(cruxFallbackReason(call), 'on purpose', 'its reason kept');
        warnCruxFallback(call, new Error('real'));
        expect(reported.join(' | '), `${call}: real`, 'the first real one reported');
        expect(cruxFallbackCount(call), 2, 'both counted');
      } finally {
        setFallbackReporterForTests(previous);
      }
    });

  test('validateMolecules: each row\'s validateMolecule text, in both engines, across chunks', async () => {
    // valid SMILES, V2000 and V3000 blocks, every kind of message, empty and blank cells, an over-long SMILES, a CXSMILES
    // query (RDKit's in both modes), SMARTS; 2,500 rows, so the Crux workers' chunks and the RDKit pool's segments cut
    // between them
    const valid = (await readDataframe('smiles.csv')).col('canonical_smiles')!.toList().filter((m) => m) as string[];
    const blocks = [...(await readDataframe('mol1K.csv')).col('molecule')!.toList().slice(0, 40),
      ...(await readDataframe('v3000_sample.csv')).col('molecule')!.toList().slice(0, 5)] as string[];
    const brokenBlock = blocks[0].replace(/^ *\d+( +\d+ .*V2000)$/m, 'xxx$1'); // its counts line's atom count broken
    const kinds = ['C1CC(', 'C1CC', 'CC)', 'c1cccc1', 'C(C)(C)(C)(C)C', '[NH4]', 'N(C)(C)(C)C', 'not a molecule', 'C[Xx]',
      'C%1', 'C1=CC=CC=C1 |junk', '', '  ', 'C'.repeat(MAX_SMILES_LENGTH + 1), '*c1ccccc1 |$Q_e;;;;;;$|', '[#6]~[#7]',
      '[NX3][CX3]=[OX1]', 'Cl[Fe](Cl)(Cl)(Cl)(Cl)(Cl)Cl', 'c1ccccc1O |$_AV:;;$|', brokenBlock];
    const molecules = Array.from({length: 2500}, (_, i) => i % 10 === 3 ? kinds[Math.floor(i / 10) % kinds.length] :
      i % 10 === 7 ? blocks[Math.floor(i / 10) % blocks.length] : valid[i % valid.length]);
    const UNREADABLE_INVALID = 'CC(=O'; // a second molecule marked unreadable for Crux below, one with a message
    for (const [row, m] of [[0, UNREADABLE_BY_CRUX], [999, UNREADABLE_INVALID], [1000, UNREADABLE_BY_CRUX],
      [1001, ''], [2499, UNREADABLE_INVALID]] as [number, string][])
      molecules[row] = m;
    const column = () => validationTextsByEngine(molecules);
    const loop = async () => molecules.map((m) => validationText(PackageFunctions.validateMolecule(m)));
    for (const engine of [ChemicalEngine.RDKit, ChemicalEngine.Crux]) {
      const [texts, perRow] = await inEngine(engine, async () => [await column(), await loop()]);
      expect(texts.length, molecules.length, `${engine}: one text per row`);
      const differ = molecules.flatMap((m, i) => texts[i] === perRow[i] ? [] :
        [`row ${i} ${m.slice(0, 30)}: ${texts[i]} vs ${perRow[i]}`]);
      expect(differ.slice(0, 3).join('\n'), '', `${engine}: validateMolecule's text on every row (${differ.length} differ)`);
      // what the texts are, so that both sides cannot be wrong alike
      expect(texts[1001], '', `${engine}: an empty cell is validateMolecule('')'s: valid`);
      expect(texts[molecules.indexOf('C'.repeat(MAX_SMILES_LENGTH + 1))], 'unknown error', `${engine}: an over-long SMILES`);
      expect(texts[molecules.indexOf('C1CC(')].startsWith('SMILES Parse Error: '), true, `${engine}: a syntax error`);
      expect(texts[molecules.indexOf('c1cccc1')].includes('kekulize'), true, `${engine}: a ring that does not kekulize`);
      expect(texts[molecules.indexOf(blocks[0])], '', `${engine}: a MOL block`);
      expect(texts.some((t) => /^\[\d\d:\d\d:\d\d\] |\n$/m.test(t)), false, `${engine}: no times, no last line break`);
      expect(brokenBlock !== blocks[0], true, 'the broken block differs');
      expect(texts.filter((t) => t === '').length > 1900, true, `${engine}: the valid molecules are valid`);
    }
    // a row the Crux workers leave (marked unreadable): RDKit's text, counted as a fallback; the other rows Crux mode's
    const rdkitMode = await inEngine(ChemicalEngine.RDKit, loop);
    const cruxMode = await inEngine(ChemicalEngine.Crux, loop);
    await withCruxUnreadable([UNREADABLE_BY_CRUX, UNREADABLE_INVALID], async () => {
      const fallbacks = cruxFallbackCount('validateMolecules');
      const texts = await inEngine(ChemicalEngine.Crux, column);
      for (const row of [0, 999, 1000, 2499])
        expect(texts[row], rdkitMode[row], `row ${row}: RDKit's text for a molecule Crux does not read`);
      expect(texts[999].includes('SMILES Parse Error'), true, 'the unreadable molecule\'s message');
      const unmarked = molecules.flatMap((m, i) => m === UNREADABLE_BY_CRUX || m === UNREADABLE_INVALID ||
        texts[i] === cruxMode[i] ? [] : [i]);
      expect(unmarked.join(','), '', 'every other row is Crux mode\'s');
      expect(cruxFallbackCount('validateMolecules'), fallbacks + 1, 'one fallback for the call');
    });
  }, {timeout: 300000});

  test('InChIs and InChI keys: a row for every molecule, the empty and malformed ones empty', async () => {
    // reported: in Crux mode Chem | Calculate | To InChI on a column with empty cells failed (row counts don't match),
    // and RDKit's InChIs for the malformed rows landed on other rows
    const molecules = ['CCO', '', 'C1CC(', 'CC(=O)O', ''];
    const inchis = async () => {
      const col = DG.Column.fromStrings('smiles', molecules);
      return [await PackageFunctions.getInchis(col), await PackageFunctions.getInchiKeys(col)]
        .map((c) => JSON.stringify(c.toList()));
    };
    const {rdkit, crux} = await inBothEngines(inchis);
    expect(crux.join('\n'), rdkit.join('\n'), 'Crux mode gives RDKit mode\'s column');
    expect(JSON.parse(crux[0]).length, molecules.length, 'one InChI per row');
    // an empty cell is no molecule for RDKit to try: no fallback
    await inEngine(ChemicalEngine.Crux, () => expectNoFallback(['InChI'], async () => {
      const col = DG.Column.fromStrings('smiles', ['CCO', '', 'CC(=O)O']);
      expect(JSON.stringify(Array.from((await PackageFunctions.getInchis(col)).toList(), (v) => !!v)), '[true,false,true]');
    }));
  });
});
