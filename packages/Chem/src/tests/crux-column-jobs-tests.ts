/* eslint-disable max-len */
/**
 * What moved to Crux in Crux mode after the first features (src/engine/README.md, "What follows the setting"): the column
 * jobs of engine/molecule-jobs.ts on the Crux workers (CoordGen coordinates, flattening, beautifying, elemental
 * analysis, the R-group dialog's cores), the reactions' canonical SMILES, `getCLogP`, the atom picker's 2D→3D mapping;
 * the CXSMILES query guard of drawing, properties and QED; the fallback report; `Chem:removeDuplicates`. Each test
 * runs both engines where the feature has both, and shows the Crux path ran by a molecule marked unreadable for Crux
 * (`withCruxUnreadable`), which RDKit then does, counted as a fallback.
 */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {PackageFunctions} from '../package';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {ChemicalEngine} from '../engine/chem-engine';
import {cruxModuleSync} from '../engine/crux-module';
import {cruxFallbackCount, setFallbackReporterForTests, warnCruxFallback} from '../engine/fallback';
import {atomCountsByEngine, atomsColumnByEngine, beautifyMolsByEngine, coordGenCoordsByEngine, coreMolblocksByCrux,
  flattenMoleculesByEngine} from '../engine/molecule-jobs';
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
      ['flattenMolecules', flattenMoleculesByEngine, (mol) => mol.get_smiles(JSON.stringify({doIsomericSmiles: false}))],
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

  test('a CXSMILES query row: RDKit draws it, Crux\'s properties and QED leave it', async () => {
    const queries = ['C1CCCCC1 |rb:0:2|', '*c1ccccc1 |m:0:1.2.3|'];
    await inEngine(ChemicalEngine.Crux, async () => {
      const renderer = new RDKitCellRenderer(rdkit) as any;
      for (const query of queries)
        expect(renderer._fetchMol(query, [], false, false, {}, false).engine, ChemicalEngine.RDKit, `drawn: ${query}`);
      expect(renderer._fetchMol('C1CCCCC1', [], false, false, {}, false).engine, ChemicalEngine.Crux, 'a plain row');
      expect(renderer._fetchMol('c1ccccc1CC', [{molecule: queries[1]}], false, false, {}, false).engine,
        ChemicalEngine.RDKit, 'a CXSMILES query scaffold');
      const col = DG.Column.fromStrings('smiles', ['CCO', ...queries]);
      const [mw] = await getPropertiesAsColumns(col, ['MW']);
      expect(mw.isNone(0), false, 'CCO\'s MW');
      expect([1, 2].map((i) => mw.isNone(i)).join(','), 'true,true', 'no Crux MW for a query row');
      expect(cruxQed(queries[0]), null, 'no Crux QED for a query row');
      expect(cruxQed('CCO') !== null, true, 'Crux\'s QED of CCO');
    });
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
