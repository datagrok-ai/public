/**
 * A CXSMILES with query features in its extension (query atom labels, ring bond count rb:, positional variation m:,
 * unsaturation u:, substitution s:, ...) in Crux mode. Crux reads every field of a CXSMILES extension as RDKit 2026.03
 * does (crux-core e881b8e98, its decision chem-wasm "CXSMILES extensions read as RDKit 2026.03 reads them"); before,
 * Chem sent such a molecule to RDKit (crux-service.ts `cruxReadsCxExtension`, removed). Each test is one kind of call
 * that moved: Crux does it now, touching no RDKit (no fallback counted, RDKit 2026.03's answer where 2024.09's differ),
 * with RDKit mode's result where both versions agree. The search, the cell renderer, properties and QED: the
 * `matchesRdkit.cxsmilesQueries` (crux-substructure-search-tests.ts) and `a CXSMILES query row` (crux-column-jobs-
 * tests.ts) tests.
 */
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';
import {PackageFunctions} from '../package';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {ChemicalEngine} from '../engine/chem-engine';
import {Fingerprint} from '../utils/chem-common';
import {fingerprintsByEngine, inchisByEngine, structuralAlertsByEngine} from '../engine/worker-jobs';
import {cruxInchi} from '../crux/crux-inchi';
import {convertNotationForColumn} from '../utils/convert-notation-utils';
import {calculateDescriptors} from '../descriptors/descriptors-engine';
import {compareValidation} from './validation-parity';
import {expectNoFallback, inBothEngines, inEngine} from './engine-utils';

/** One CXSMILES per kind of query feature, and each one's molecule without its extension. */
const QUERIES = ['*c1ccccc1 |$Q_e;;;;;;$|', 'C1CCCCC1 |rb:0:2|', '*c1ccccc1 |m:0:1.2.3|', 'CCC |u:1|', 'CCC |s:1:1|'];
const MOLECULES = ['*c1ccccc1', 'C1CCCCC1', '*c1ccccc1', 'CCC', 'CCC'];
const ASPIRIN = 'CC(=O)Oc1ccccc1C(=O)O';

/** The number of pixels a drawing inked. */
function inked(image: ImageData): number {
  let n = 0;
  for (let i = 3; i < image.data.length; i += 4)
    n += image.data[i] > 0 ? 1 : 0;
  return n;
}

category('chemical engine: CXSMILES query features', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('validateMolecule: Crux writes the log of a CXSMILES with query features, RDKit 2026.03\'s', async () => {
    // a SMILES error before the extension: RDKit 2026.03's lines (RDKit-js 2026.03.6) say where it is, Chem's RDKit
    // 2024.09's do not; the queries read without a message in both
    const broken = 'C1CC( |rb:0:2|';
    const {rdkit, crux} = await inBothEngines(async () =>
      [broken, ...QUERIES].map((s) => PackageFunctions.validateMolecule(s)), {noFallback: ['validateMolecule']});
    expect((crux[0] ?? '').includes('SMILES Parse Error: check for mistakes around position 5:'), true,
      `Crux's log, RDKit 2026.03's lines: ${crux[0]}`);
    expect(compareValidation(rdkit[0], crux[0]) !== null, true, 'RDKit mode\'s but for the versions\' differences');
    expect(crux.slice(1).map((m) => m ?? '').join('|'), QUERIES.map(() => '').join('|'), 'the queries read');
  });

  test('conversions: Crux converts a CXSMILES with query features, as RDKit 2026.03 writes it', async () => {
    const column = async () => convertNotationForColumn(DG.Column.fromStrings('m', QUERIES), DG.chem.Notation.Smiles);
    const {rdkit, crux} = await inBothEngines(async () => ({
      smiles: QUERIES.map((s) =>
        PackageFunctions.convertMolNotation(s, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles)),
      column: await column(),
      block: PackageFunctions.convertMolNotation(QUERIES[2], DG.chem.Notation.Unknown, DG.chem.Notation.MolBlock),
    }), {noFallback: ['convertMolNotation']});
    expect(crux.smiles.join(' '), MOLECULES.join(' '), 'Crux\'s SMILES: the molecules');
    expect(crux.smiles.join(' '), rdkit.smiles.join(' '), 'RDKit mode\'s SMILES');
    expect(crux.column.join(' '), rdkit.column.join(' '), 'the column: RDKit mode\'s SMILES, by the Crux workers');
    // the positional variation's dummy as RDKit 2026.03 writes it in a MOL block (with an atom list, `M  ALS`), where
    // Chem's RDKit 2024.09 writes the symbol `A` and no list
    expect(crux.block.includes('M  ALS'), true, `Crux's MOL block, RDKit 2026.03's: ${crux.block}`);
    expect(rdkit.block.includes('M  ALS'), false, 'RDKit mode\'s MOL block, RDKit 2024.09\'s');
  });

  test('InChI, fingerprints and alerts: the Crux workers do a CXSMILES with query features, as RDKit', async () => {
    const molecules = [...QUERIES, ASPIRIN];
    const {rdkit, crux} = await inBothEngines(async () => {
      const fps = await fingerprintsByEngine(Fingerprint.Pattern, molecules, false);
      const alerts = await structuralAlertsByEngine({'PAINS': ['c1ccccc1', '[CH2]']}, molecules);
      return {
        inchis: await inchisByEngine(molecules, false), keys: await inchisByEngine(molecules, true),
        engine: fps.engine, fps: fps.fps.map((fp) => fp?.toString() ?? ''),
        alerts: alerts.map(([rule, flags]) => `${rule} ${flags.join(',')}`).join('; '),
      };
    }, {noFallback: ['InChI', 'fingerprints', 'structural alerts']});
    expect(crux.inchis.join('|'), rdkit.inchis.join('|'), 'RDKit\'s InChIs');
    expect(crux.keys.join('|'), rdkit.keys.join('|'), 'RDKit\'s InChIKeys');
    expect(crux.inchis[1], 'InChI=1S/C6H12/c1-2-4-6-5-3-1/h1-6H2', 'the ring bond count query\'s molecule');
    expect(crux.engine, ChemicalEngine.Crux, 'Crux computed the fingerprints');
    expect(crux.fps.join('|'), rdkit.fps.join('|'), 'RDKit\'s pattern fingerprints');
    expect(crux.alerts, rdkit.alerts, 'RDKit\'s matches');
    await inEngine(ChemicalEngine.Crux, () => expectNoFallback(['InChI'], async () =>
      expect((await cruxInchi(QUERIES[1]))?.inchi, rdkit.inchis[1], 'the Identifiers panel\'s InChI: Crux\'s')));
  });

  test('descriptors and Murcko scaffolds: Crux\'s for a CXSMILES with query features', async () => {
    await inEngine(ChemicalEngine.Crux, () => expectNoFallback(['calcDescriptors', 'Murcko Scaffolds'], async () => {
      // the query feature leaves the molecule as it is: the values of the molecule without it
      const col = DG.Column.fromStrings('smiles', [QUERIES[1], MOLECULES[1], QUERIES[3], MOLECULES[3]]);
      col.semType = DG.SEMTYPE.MOLECULE;
      const values = await calculateDescriptors(col, ['MolWt', 'HeavyAtomCount', 'RingCount']);
      for (const c of values) {
        expect(`${c.get(0)} ${c.get(2)}`, `${c.get(1)} ${c.get(3)}`, `${c.name} of the queries`);
        expect(c.isNone(0) || c.isNone(2), false, `${c.name}: Crux's values`);
      }
      const df = DG.DataFrame.fromColumns([col]);
      const result = await PackageFunctions.MurckoScaffolds(df, col); // eslint-disable-line new-cap
      const scaffolds = df.col('scaffolds') ?? result.col('scaffolds');
      expect(`${scaffolds?.get(0)}`, `${scaffolds?.get(1)}`, 'the ring bond count query\'s scaffold');
      expect(`${scaffolds?.get(0)}`, 'C1CCCCC1', 'cyclohexane');
    }));
  });

  test('drawing: canvasMol draws a CXSMILES with query features with Crux', async () => {
    await inEngine(ChemicalEngine.Crux, () => expectNoFallback(['canvasMol', 'rendering', 'drawing'], async () => {
      for (const query of QUERIES) {
        const canvas = document.createElement('canvas');
        chemCommonRdKit.drawMoleculeToCanvas(0, 0, 200, 100, canvas, query, null);
        expect(inked(canvas.getContext('2d')!.getImageData(0, 0, canvas.width, canvas.height)) > 0, true,
          `canvasMol inked ${query}`);
      }
      const canvas = document.createElement('canvas');
      chemCommonRdKit.drawMoleculeToCanvas(0, 0, 200, 100, canvas, 'c1ccccc1CC', QUERIES[2]);
      expect(inked(canvas.getContext('2d')!.getImageData(0, 0, canvas.width, canvas.height)) > 0, true,
        'canvasMol inked a molecule with a query scaffold');
    }));
  });
});
