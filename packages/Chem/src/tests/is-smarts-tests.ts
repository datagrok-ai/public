import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';
import {PackageFunctions} from '../package';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {before, expect, category, test, awaitCheck} from '@datagrok-libraries/test/src/test';
import {isFragment, _isSmarts} from '../utils/mol-creation_rdkit';
// the module itself for the functions of 2026-10-14 (the prefilter), absent in the code before them
import * as molCreation from '../utils/mol-creation_rdkit';
import {CustomSketcherInput} from '../widgets/col-highlights';
import {renderMolecule} from '../rendering/render-molecule';
import {prepareReactionDepiction} from '../rendering/reaction-depiction';
import {SMARTS_MOLECULE_MESSAGE} from '../constants';
import {readDataframe} from './utils';

const smarts = [
  `
    RDKit          2D

 0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 "NOT [C,N]" -1.942000 0.958900 0.000000 0
M  V30 2 C -2.656400 0.546400 0.000000 0
M  V30 3 C -2.656400 -0.278600 0.000000 0
M  V30 4 C -1.942000 -0.691100 0.000000 0
M  V30 5 C -1.227500 -0.278600 0.000000 0
M  V30 6 C -1.227500 0.546400 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 2 3 4
M  V30 3 1 4 5
M  V30 4 2 5 6
M  V30 5 2 1 2
M  V30 6 1 6 1
M  V30 END BOND
M  V30 END CTAB
M  END
  `,
  `
MJ201900

  6  6  1  0  0  0  0  0  0  0999 V2000
   -1.9420    0.9589    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6564    0.5464    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6564   -0.2786    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9420   -0.6911    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2275   -0.2786    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2275    0.5464    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  1  2  2  0  0  0  0
  6  1  1  0  0  0  0
  1 T    2   6   7
M  ALS   1  2 T C   N
M  END
`,
  `[!#6&!#7]1:[#6]:[#6]:[#6]:[#6]:[#6]:1`,
];

const fragments = [
  `
  RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 R# 1.500000 0.000000 0.000000 0
M  V30 2 C 0.750000 -1.299038 0.000000 0
M  V30 3 C -0.750000 -1.299038 0.000000 0
M  V30 4 C -1.500000 0.000000 0.000000 0
M  V30 5 C -0.750000 1.299038 0.000000 0
M  V30 6 C 0.750000 1.299038 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 END BOND
M  V30 END CTAB
M  END
`,
  `
  RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 * 1.500000 0.000000 0.000000 0
M  V30 2 C 0.750000 -1.299038 0.000000 0
M  V30 3 C -0.750000 -1.299038 0.000000 0
M  V30 4 C -1.500000 0.000000 0.000000 0
M  V30 5 C -0.750000 1.299038 0.000000 0
M  V30 6 C 0.750000 1.299038 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 END BOND
M  V30 END CTAB
M  END
`,
  `
MJ201900                      

  6  6  0  0  0  0  0  0  0  0999 V2000
    7.5467   -7.6478    0.0000 R   0  0  0  0  0  0  0  0  0  0  0  0
    8.9731   -7.6474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2613   -7.2355    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.9731   -8.4726    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5467   -8.4763    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2631   -8.8844    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  3  1  2  0  0  0  0
  4  2  2  0  0  0  0
  1  5  1  0  0  0  0
  2  3  1  0  0  0  0
  5  6  2  0  0  0  0
  6  4  1  0  0  0  0
M  END`,
  `
MJ201900                      

  6  6  0  0  0  0  0  0  0  0999 V2000
    7.5467   -7.6478    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
    8.9731   -7.6474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2613   -7.2355    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.9731   -8.4726    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5467   -8.4763    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2631   -8.8844    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  3  1  2  0  0  0  0
  4  2  2  0  0  0  0
  1  5  1  0  0  0  0
  2  3  1  0  0  0  0
  5  6  2  0  0  0  0
  6  4  1  0  0  0  0
M  END`,
  `CC1(N[*:1])CCN([*:2])C1`,
];

const notSmarts = [
  'O=C1CN=C(c2ccccc2N1)C3CCCCC3',
  `
-ISIS-  03062316482D

 27 29  0  0  1  0  0  0  0  0999 V2000
    4.5000    0.3958    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    4.9333    1.2958    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.9083    1.5167    0.0000 N   0  0  3  0  0  0  0  0  0  0  0  0
    5.9083   -0.7292    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    4.9333   -0.5083    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    6.6875    0.8958    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.6875   -0.1042    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5000    0.3958    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3083   -1.2875    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    1.2583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3083    2.0750    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5000    2.1250    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6750   -2.2208    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3208   -1.1417    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1292    2.4917    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0625   -2.8542    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6958   -1.9208    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0500   -3.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6167    4.1417    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    7.5833    1.4292    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4375   -3.6333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.5833   -0.6417    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    1.2583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3958    3.1708    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8042   -4.5667    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.4875    0.9167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.4875   -0.1250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
  3  2  1  0  0  0  0
  4  5  1  0  0  0  0
  5  1  1  0  0  0  0
  6  3  1  0  0  0  0
  7  4  1  0  0  0  0
  1  8  1  1  0  0  0
  5  9  1  1  0  0  0
 10  8  1  0  0  0  0
 11  2  2  0  0  0  0
 12 10  2  0  0  0  0
 13  9  1  0  0  0  0
 14  9  2  0  0  0  0
 15  3  1  0  0  0  0
 16 18  1  0  0  0  0
 17 14  1  0  0  0  0
 18 13  2  0  0  0  0
 19 24  1  0  0  0  0
 20  6  1  0  0  0  0
 21 16  1  0  0  0  0
 22  7  1  0  0  0  0
 23 10  1  0  0  0  0
 24 15  1  0  0  0  0
 25 21  1  0  0  0  0
 26 20  2  0  0  0  0
 27 22  2  0  0  0  0
  7  6  2  0  0  0  0
 16 17  2  0  0  0  0
 27 26  1  0  0  0  0
M  END
`,
  `
    RDKit          2D

 0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 22 24 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -5.494284 1.062496 0.000000 0
M  V30 2 C -5.052151 -0.370863 0.000000 0
M  V30 3 C -3.589759 -0.704645 0.000000 0
M  V30 4 C -3.147627 -2.138004 0.000000 0
M  V30 5 N -1.685235 -2.471785 0.000000 0
M  V30 6 C -1.351453 -3.934177 0.000000 0
M  V30 7 O -2.524201 -4.869412 0.000000 0
M  V30 8 C -0.000000 -4.585003 0.000000 0
M  V30 9 N 1.351453 -3.934177 0.000000 0
M  V30 10 C 1.685235 -2.471785 0.000000 0
M  V30 11 C 0.750000 -1.299038 0.000000 0
M  V30 12 C 1.500000 0.000000 0.000000 0
M  V30 13 C 0.750000 1.299038 0.000000 0
M  V30 14 C -0.750000 1.299038 0.000000 0
M  V30 15 C -1.500000 0.000000 0.000000 0
M  V30 16 C -0.750000 -1.299038 0.000000 0
M  V30 17 C 3.147627 -2.138004 0.000000 0
M  V30 18 C 4.167886 -3.237582 0.000000 0
M  V30 19 C 5.630278 -2.903800 0.000000 0
M  V30 20 C 6.072410 -1.470441 0.000000 0
M  V30 21 C 5.052151 -0.370863 0.000000 0
M  V30 22 C 3.589759 -0.704645 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 2 6 7
M  V30 7 1 6 8
M  V30 8 1 8 9
M  V30 9 2 9 10
M  V30 10 1 10 11
M  V30 11 2 11 12
M  V30 12 1 12 13
M  V30 13 2 13 14
M  V30 14 1 14 15
M  V30 15 2 15 16
M  V30 16 1 10 17
M  V30 17 1 17 18
M  V30 18 1 18 19
M  V30 19 1 19 20
M  V30 20 1 20 21
M  V30 21 1 21 22
M  V30 22 1 16 5
M  V30 23 1 16 11
M  V30 24 1 22 17
M  V30 END BOND
M  V30 END CTAB
M  END
`,
  `
     RDKit          2D

 22 24  0  0  0  0  0  0  0  0999 V2000
   -5.4943    1.0625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0522   -0.3709    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5898   -0.7046    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1476   -2.1380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6852   -2.4718    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3515   -3.9342    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5242   -4.8694    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -4.5850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3515   -3.9342    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.6852   -2.4718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1476   -2.1380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1679   -3.2376    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6303   -2.9038    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0724   -1.4704    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0522   -0.3709    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5898   -0.7046    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  5  6  1  0
  6  7  2  0
  6  8  1  0
  8  9  1  0
  9 10  2  0
 10 11  1  0
 11 12  2  0
 12 13  1  0
 13 14  2  0
 14 15  1  0
 15 16  2  0
 10 17  1  0
 17 18  1  0
 18 19  1  0
 19 20  1  0
 20 21  1  0
 21 22  1  0
 16  5  1  0
 16 11  1  0
 22 17  1  0
M  END
`,
];

category('isSmarts', async () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('isSmarts', async () => {
    smarts.forEach((mol) => {
      expect(_isSmarts(mol), true);
    });
    notSmarts.forEach((mol) => expect(_isSmarts(mol), false));
  });

  test('isFragment', async () => {
    fragments.forEach((mol) => {
      expect(isFragment(mol), true);
    });
    notSmarts.forEach((mol) => expect(isFragment(mol), false));
  });

  test('benchmark.isSmarts.100000', async () => {
    const molfiles = await readDataframe('tests/molfiles_100.csv');
    const molecules = molfiles.col('molfile')!.categories;
    console.profile('isSmarts');
    for (let i = 0; i < 100000; i++)
      _isSmarts(molecules[i%100]);
    console.profileEnd('isSmarts');
  });

  // Chem:isSmarts (PackageFunctions.isSmarts) and its trails, 2026-10-14: the product owner's call on fix A
  // ("let's follow the trails and make sure all is good"; no heavy loops on an RDKit parse). Each test names its trail.

  test('isSmarts: a SMARTS no syntax tells from SMILES is SMARTS', async () => {
    for (const s of SMARTS_OF_PRIMITIVES)
      expect(PackageFunctions.isSmarts(s), true, `${s} is SMARTS`);
  });

  test('isSmarts: a SMILES, a CXSMILES (stereo groups, atom labels, coordinates) is no SMARTS', async () => {
    for (const s of [...SMILES, ...CXSMILES])
      expect(PackageFunctions.isSmarts(s), false, `${s} is no SMARTS`);
  });

  // js-api: DG.chem.isSmarts finds Chem:isSmarts and calls it, for DG.chem.Sketcher.setMolecule and Retrosynthesis
  test('isSmarts: DG.chem.isSmarts answers as Chem:isSmarts', async () => {
    for (const s of [...SMARTS_OF_PRIMITIVES, ...SMILES, ...CXSMILES])
      expect(DG.chem.isSmarts(s), SMARTS_OF_PRIMITIVES.includes(s), s);
  });

  // RDKit is asked only about a text outside SMILES syntax: it must let through every string RDKit reads as a query
  // only (a false negative would set a SMARTS as SMILES) and no SMILES of the demo files (each would cost a parse)
  test('isSmarts: the SMILES-syntax prefilter skips no SMARTS and parses no SMILES of the demo files', async () => {
    const isSmilesSyntax = (molCreation as any).isSmilesSyntax as ((s: string) => boolean) | undefined;
    const readsAsQueryOnly = (molCreation as any).readsAsQueryOnly as ((s: string, m: any) => boolean) | undefined;
    if (!isSmilesSyntax || !readsAsQueryOnly)
      throw new Error('no isSmilesSyntax or readsAsQueryOnly in mol-creation_rdkit');
    const rdkit = chemCommonRdKit.getRdKitModule();
    const smarts = [...SMARTS_OF_PRIMITIVES,
      ...(await readDataframe('alert-collection.csv')).col('smarts')!.toList(),
      ...(await readDataframe('pharmacophore-features.csv')).col('smarts')!.toList()];
    const skipped = smarts.filter((s) => !_isSmarts(s) && readsAsQueryOnly(s, rdkit) && isSmilesSyntax(s));
    expect(skipped.length, 0, `read as SMILES, being SMARTS: ${skipped.slice(0, 5).join(' ')}`);
    const demo = await demoMoleculeStrings();
    const parsed = demo.filter((s) => !DG.chem.isMolBlock(s) && !_isSmarts(s) &&
      !isSmilesSyntax(s.replace(/\s+\|[^|]*\|\s*$/, '')));
    expect(parsed.length, 0, `SMILES parsed by RDKit: ${parsed.slice(0, 5).join(' ')}`);
  });

  test('benchmark.PackageFunctions.isSmarts.100000', async () => {
    const strings = await demoMoleculeStrings();
    const t0 = performance.now();
    for (let i = 0; i < 100000; i++)
      PackageFunctions.isSmarts(strings[i % strings.length]);
    const ms = performance.now() - t0;
    // the same strings without the prefilter, every string outside SMARTS syntax parsed by RDKit: 5,000, scaled
    const readsAsQueryOnly = (molCreation as any).readsAsQueryOnly as ((s: string, m: any) => boolean) | undefined;
    let unfiltered = NaN;
    if (readsAsQueryOnly) {
      const rdkit = chemCommonRdKit.getRdKitModule();
      const t1 = performance.now();
      for (let i = 0; i < 5000; i++) {
        const s = strings[i % strings.length];
        void (!/M  END|V2000|V3000/.test(s) && (_isSmarts(s) || readsAsQueryOnly(s, rdkit)));
      }
      unfiltered = (performance.now() - t1) * 20;
    }
    expect(ms < 3000, true, `100,000 calls took ${ms.toFixed(0)} ms`);
    return `100,000 strings of ${strings.length}: ${ms.toFixed(0)} ms; ` +
      `without the prefilter ~${unfiltered.toFixed(0)} ms (5,000 timed, x 20)`;
  });

  // js-api: DG.chem.Sketcher.setMolecule, behind the field's Enter, Recent and Favourites, the substructure filter, the
  // similarity viewer, the molecule cell editor and the pinned column's editor
  test('isSmarts: the sketcher is set with a SMARTS of primitives as SMARTS, a CXSMILES as a molecule', async () => {
    const sketcher = new DG.chem.Sketcher();
    try {
      for (const s of SMARTS_OF_PRIMITIVES) {
        sketcher.setMolecule(s);
        expect((sketcher as any)._smarts, s, `${s} set as SMARTS`);
      }
      for (const s of CXSMILES) {
        sketcher.setMolecule(s);
        expect((sketcher as any)._smiles, s, `${s} set as SMILES`);
      }
    } finally {
      sketcher.detach();
    }
  });

  // the column highlights' sketcher input, which also asks at each change
  test('isSmarts: a column highlight\'s SMARTS of primitives is set as SMARTS', async () => {
    const input = new CustomSketcherInput('[NX3][CX3]=[OX1]');
    try {
      expect((input.sketcher as any)._smarts, '[NX3][CX3]=[OX1]');
    } finally {
      input.sketcher.detach();
    }
  });

  // _isSmarts trails: the panels that do not apply to a query (SMARTS_MOLECULE_MESSAGE); two that compute in the
  // browser (Descriptors asks a container, which a dummy atom fails)
  test('isSmarts: a panel computes a CXSMILES with a stereo group or labels, not a SMARTS of any bonds', async () => {
    const said = (w: DG.Widget) => w.root.textContent?.includes(SMARTS_MOLECULE_MESSAGE) ?? false;
    const molecule = (s: string) => DG.SemanticValue.fromValueType(s, DG.SEMTYPE.MOLECULE);
    for (const s of ['C[C@H](N)C(=O)O |&1:1|', 'Cc1ccccc1 |$;;;;Q_e;;$|']) {
      expect(said(PackageFunctions.properties(molecule(s))), false, `Properties of ${s}`);
      expect(said(PackageFunctions.toxicity(molecule(s))), false, `Toxicity of ${s}`);
    }
    expect(said(PackageFunctions.properties(molecule('c1ccccc1~O'))), true, 'Properties of c1ccccc1~O');
  });

  // _isSmarts trails: Copy as SMILES (the cell's action) and Copy SMILES (the molecule menu of a rendered molecule)
  // give a SMILES: a CXSMILES loses its extension
  test('isSmarts: Copy as SMILES and the molecule menu\'s Copy SMILES drop a CXSMILES extension', async () => {
    const rdkit = chemCommonRdKit.getRdKitModule();
    const smilesOf = (s: string) => {
      const mol = rdkit.get_mol(s);
      try {
        return mol.get_smiles();
      } finally {
        mol.delete();
      }
    };
    await withClipboard(async (copied) => {
      for (const s of COPIED_CXSMILES) {
        PackageFunctions.copyAsSmiles(DG.SemanticValue.fromValueType(s, DG.SEMTYPE.MOLECULE));
        expect(copied.pop(), smilesOf(s), `Copy as SMILES of ${s}`);
      }
      for (const s of COPIED_CXSMILES) {
        const box = renderMolecule(s, {width: 200, height: 100});
        document.body.append(box);
        try {
          box.querySelector('canvas')!.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true}));
          const item = () => Array.from(document.querySelectorAll<HTMLElement>('.d4-menu-item-label'))
            .find((el) => el.textContent?.trim() === 'Copy SMILES');
          await awaitCheck(() => item() !== undefined, 'the molecule menu has no Copy SMILES', 5000);
          item()!.click();
          await awaitCheck(() => copied.length > 0, 'Copy SMILES copied nothing', 5000);
          expect(copied.pop(), smilesOf(s), `Copy SMILES of ${s}`);
        } finally {
          box.remove();
          for (const el of Array.from(document.querySelectorAll('.d4-menu-popup')))
            el.remove();
        }
      }
    });
  });

  // _isSmarts trail: Copy as SMARTS gives RDKit's SMARTS of a CXSMILES (it gave the CXSMILES as it was, its `&` or `$`
  // read as SMARTS)
  test('isSmarts: Copy as SMARTS gives a CXSMILES with a stereo group as RDKit\'s SMARTS', async () => {
    await withClipboard(async (copied) => {
      const s = 'C[C@H](N)C(=O)O |&1:1|';
      PackageFunctions.copyAsSmarts(DG.SemanticValue.fromValueType(s, DG.SEMTYPE.MOLECULE));
      const smarts = copied.pop() ?? '';
      expect(smarts.includes('|'), false, `Copy as SMARTS of ${s}: ${smarts}`);
    });
  });

  // _isSmarts trail: getQueryMolSafe's treatAsSmarts, every substructure search of a string: a CXSMILES with a stereo
  // group was searched as SMARTS (its [C@H] an H count, its C aliphatic only), now as its SMILES is
  test('isSmarts: a substructure search by a CXSMILES with a stereo group finds what its SMILES finds', async () => {
    const df = await readDataframe('smiles.csv');
    const col = df.col('canonical_smiles')!;
    const plain = await grok.chem.searchSubstructure(col, 'C[C@H](O)[C@@H](C)CC');
    const cx = await grok.chem.searchSubstructure(col, 'C[C@H](O)[C@@H](C)CC |&1:1,3|');
    expect(plain.trueCount > 0, true, 'the SMILES finds rows');
    expect(cx.trueCount, plain.trueCount, 'the CXSMILES finds what its SMILES finds');
  });

  // _isSmarts trail: a reaction is depicted from its molecules only when no side is SMARTS (once per reaction drawn): a
  // side with an any bond was read as a molecule (RDKit's SMILES reader takes ~), now as the SMARTS it is
  test('isSmarts: a reaction with an any bond is not depicted as molecules', async () => {
    expect(prepareReactionDepiction(chemCommonRdKit.getRdKitModule(), 'C~C>>CC', new Map()), null);
  });
});

/** SMARTS that hold no `[#n]`, `$()`, `&`, `;`, `,`, `!` or `~`: RDKit reads each as a query only. */
const SMARTS_OF_PRIMITIVES = ['[NX3][CX3]=[OX1]', '[R]', 'A-A', 'C@C', 'a1aaaaa1', '[CH2X4]', '[N+H3]', '[C-H]', '[H2]',
  'C=\\C'];
const SMILES = ['c1ccccc1', '[NH4+]', '[2H]C', 'C->[Fe]', '[*:1]C', 'CC[C@@H](C)O', ''];
const CXSMILES = ['C[C@H](N)C(=O)O |&1:1|', 'C[C@H](O)[C@@H](C)CC |&1:1,3|', 'Cc1ccccc1 |$;;;;Q_e;;$|',
  'CCO |(0,0,;1,0,;2,0,)|', 'C[C@H](N)C(=O)O |o1:1|', 'C1CC1 |r|'];
/** CXSMILES whose extension Copy SMILES dropped (`;`, `,`, `&` or `$` in it read as SMARTS) or kept (`|o1:1|`,
 * `|r|`). */
const COPIED_CXSMILES = ['C[C@H](N)C(=O)O |&1:1|', 'Cc1ccccc1 |$;;;;Q_e;;$|', 'CCO |(0,0,;1,0,;2,0,)|',
  'C[C@H](N)C(=O)O |o1:1|', 'C1CC1 |r|'];

/** The demo files' molecule strings: smiles.csv's SMILES, mol1K's molblocks and RDKit's SMILES of them, spgi-100's
 * molblocks. */
async function demoMoleculeStrings(): Promise<string[]> {
  const rdkit = chemCommonRdKit.getRdKitModule();
  const mol1K = (await readDataframe('mol1K.csv')).col('molecule')!.toList();
  const mol1KSmiles = mol1K.map((mb) => {
    const mol = rdkit.get_mol(mb);
    try {
      return mol?.get_smiles() ?? '';
    } finally {
      mol?.delete();
    }
  }).filter((s) => s);
  return [...(await readDataframe('smiles.csv')).col('canonical_smiles')!.toList(), ...mol1K, ...mol1KSmiles,
    ...(await readDataframe('tests/spgi-100.csv')).col('Structure')!.toList()].filter((s) => s);
}

/** Runs `f` with navigator.clipboard.writeText collecting what is written, in place of the clipboard. */
async function withClipboard(f: (copied: string[]) => Promise<void>): Promise<void> {
  const copied: string[] = [];
  const clipboard = navigator.clipboard as any;
  Object.defineProperty(clipboard, 'writeText', {configurable: true, writable: true,
    value: async (text: string) => {copied.push(text);}});
  try {
    await f(copied);
  } finally {
    delete clipboard.writeText;
  }
}
