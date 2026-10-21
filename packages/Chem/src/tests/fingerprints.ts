import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {before, category, expect, expectFloat, test} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {chemGetFingerprints} from '../chem-searches';
import {Fingerprint, rdKitFingerprintToBitArray} from '../utils/chem-common';
import {chemSimilaritySearch} from '../analysis/chem-similarity-viewer';
import {tanimotoSimilarity} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import {BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import { readDataframe } from './utils';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {ChemicalEngine, cruxReady, isCruxRequested, setChemicalEngine} from '../engine/chem-engine';
import {fingerprintsKey} from '../engine/worker-jobs';
import {FingerprintOptions} from '../utils/fingerprint-options';

const ASPIRIN = 'CC(=O)Oc1ccccc1C(=O)O';
const CAFFEINE = 'Cn1cnc2c1c(=O)n(C)c(=O)n2C';
const L_ALANINE = 'C[C@H](N)C(=O)O';
const D_ALANINE = 'C[C@@H](N)C(=O)O';
/** The bits RDKit 2026.03.6 sets (Python, rdFingerprintGenerator; the script: chem-crux-more's py/fcfp_expected.py). */
const RDKIT_BITS = {
  aspirinEcfp6: '389,456,541,589,650,695,740,807,909,1017,1027,1035,1047,1057,1088,1106,1199,1380,1410,1447,1467,' +
    '1468,1616,1729,1750,1775,1865,1873,1917,1970,1991',
  aspirinEcfp4At1024: '11,23,33,64,175,356,386,389,423,444,456,592,650,695,705,726,751,807,849,893,909,946,967,1017',
  aspirinFcfp4: '0,1,2,4,32,33,321,535,598,626,809,1001,1048,1083,1085,1405,1502,1950,1979,1990,2017,2040',
  caffeineFcfp4: '0,2,4,6,45,317,411,465,476,497,570,886,966,984,1155,1177,1220,1223,1344,1397,1863,1932',
  lAlanineChiral: '1,283,389,650,786,807,1057,1112,1171,1187,1844,1917',
  dAlanineChiral: '1,46,283,389,650,786,807,1057,1113,1171,1844,1917',
};

/** The indexes of the set bits. */
function onBits(fp: BitArray | null): string {
  const bits: number[] = [];
  for (let i = -1; fp && (i = fp.findNext(i)) !== -1;)
    bits.push(i);
  return bits.join(',');
}

/** A one-column table of molecules. */
function molTable(smiles: string[]): DG.Column {
  const col = DG.Column.fromStrings('smiles', smiles);
  DG.DataFrame.fromColumns([col]).name = 'fingerprint options';
  col.semType = DG.SEMTYPE.MOLECULE;
  return col;
}

/** The fingerprint caches of the column's table, sorted. */
function cachedColumns(col: DG.Column): string {
  return col.dataFrame.columns.names().filter((n) => n.startsWith('~')).sort().join(', ');
}

/** Column names as {@link cachedColumns} lists them. */
function sorted(...names: string[]): string {
  return names.sort().join(', ');
}

/** The fingerprints of molecules in a fresh table (no caches). */
async function freshBits(smiles: string[], type: Fingerprint, options?: FingerprintOptions): Promise<string[]> {
  return (await chemGetFingerprints(molTable(smiles), type, false, options)).map(onBits);
}

/** The cache key suffix of the engine in effect. */
function engineKey(): string {
  return fingerprintsKey(isCruxRequested() ? ChemicalEngine.Crux : ChemicalEngine.RDKit);
}

/** `check` with the engine pinned to each engine in turn (the stand's setting after). */
async function inBothEngines(check: (engine: ChemicalEngine) => Promise<void>): Promise<void> {
  try {
    for (const engine of [ChemicalEngine.RDKit, ChemicalEngine.Crux]) {
      setChemicalEngine(engine);
      if (engine === ChemicalEngine.Crux)
        expect(await cruxReady(), true, 'the Crux module loads');
      await check(engine);
    }
  } finally {
    setChemicalEngine(null);
  }
}

category('fingerprints', async () => {
  let molecules: DG.DataFrame;
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    molecules = DG.Test.isInBenchmark ? await readDataframe('tests/smi10K.csv') :
      grok.data.demo.molecules(100);
  });

  test('Pattern', async () => {
    await chemGetFingerprints(molecules.col('smiles')!, Fingerprint.Pattern, false)!;
  }, {benchmark: true});

  test('Morgan', async () => {
    await chemGetFingerprints(molecules.col('smiles')!, Fingerprint.Morgan, false)!;
  }, {benchmark: true});

  test('similarity search on the other fingerprints', async () => {
    const col = molecules.col('smiles')!;
    const rdkit = chemCommonRdKit.getRdKitModule();
    const fingerprintOf: {[fp: string]: (mol: RDMol) => Uint8Array} = {
      [Fingerprint.RDKit]: (mol) => mol.get_rdkit_fp_as_uint8array(),
      [Fingerprint.MACCS]: (mol) => mol.get_maccs_fp_as_uint8array(),
      [Fingerprint.AtomPair]: (mol) => mol.get_atom_pair_fp_as_uint8array(),
      [Fingerprint.TopologicalTorsion]: (mol) => mol.get_topological_torsion_fp_as_uint8array(),
    };
    const bits = (smiles: string, fp: string) => {
      const mol = rdkit.get_mol(smiles);
      try {
        return rdKitFingerprintToBitArray(fingerprintOf[fp](mol))!;
      } finally {
        mol.delete();
      }
    };
    const scores = async (fp: Fingerprint) =>
      (await chemSimilaritySearch(col, col.get(0), BitArrayMetricsNames.Tanimoto, 12, 0.01, fp))!;
    const morgan = (await scores(Fingerprint.Morgan)).getCol('score').toList();
    for (const fp of Object.keys(fingerprintOf) as Fingerprint[]) {
      const fps = await chemGetFingerprints(col, fp, false);
      expect(fps.filter((f) => f == null).length, 0, `${fp}: molecules without a fingerprint`);
      const res = await scores(fp);
      expect(res.rowCount, 12, fp);
      expect(res.get('score', 0), 1, fp);
      const query = bits(col.get(0), fp);
      for (let i = 0; i < res.rowCount; i++)
        expectFloat(res.get('score', i), tanimotoSimilarity(query, bits(col.get(res.get('indexes', i)), fp)), 1e-4,
          `${fp}: the score of row ${res.get('indexes', i)}`);
      expect(res.getCol('score').toList().some((s, i) => Math.abs(s - morgan[i]) > 1e-4), true,
        `${fp} scores the same as Morgan`);
    }
  });

  test('an edit reaches every cached fingerprint column', async () => {
    // each cache has its own edited rows: the first cache read used up the column's one change set, so the others kept
    // the old molecule's fingerprint (here Morgan, read after AtomPair)
    const col = molTable([ASPIRIN, 'c1ccccc1', CAFFEINE]);
    await chemGetFingerprints(col, Fingerprint.Morgan, false);
    await chemGetFingerprints(col, Fingerprint.AtomPair, false);
    col.set(1, 'c1ccncc1');
    const atomPairs = (await chemGetFingerprints(col, Fingerprint.AtomPair, false)).map(onBits);
    const morgan = (await chemGetFingerprints(col, Fingerprint.Morgan, false)).map(onBits);
    const molecules = col.toList();
    expect(atomPairs[1], (await freshBits(molecules, Fingerprint.AtomPair))[1], 'AtomPair: the edited row');
    expect(morgan[1], (await freshBits(molecules, Fingerprint.Morgan))[1], 'Morgan: the edited row AtomPair read first');
    expect(morgan[0], (await freshBits(molecules, Fingerprint.Morgan))[0], 'Morgan: a row not edited');
  });

  test('an edit reaches a cached column saved before it was read (a project\'s)', async () => {
    // a table holding a hidden ~smiles.AtomPair column, as a project brings it, that no read has tracked yet: Morgan
    // is read first (the column's edits are tracked from then on), a row is edited, then AtomPair is read
    const molecules = [ASPIRIN, 'c1ccccc1', CAFFEINE];
    const source = molTable(molecules);
    await chemGetFingerprints(source, Fingerprint.AtomPair, false);
    const saved = source.dataFrame.columns.toList().find((c) => c.name.startsWith('~smiles.AtomPair'))!;
    const col = DG.Column.fromStrings('smiles', molecules);
    const copy = DG.Column.fromType(DG.COLUMN_TYPE.BYTE_ARRAY, saved.name, saved.length).init((i) => saved.get(i));
    DG.DataFrame.fromColumns([col, copy]).name = 'a project';
    col.semType = DG.SEMTYPE.MOLECULE;
    await chemGetFingerprints(col, Fingerprint.Morgan, false);
    col.set(1, 'c1ccncc1');
    const atomPairs = (await chemGetFingerprints(col, Fingerprint.AtomPair, false)).map(onBits);
    const fresh = await freshBits(col.toList(), Fingerprint.AtomPair);
    expect(atomPairs[1], fresh[1], 'the saved AtomPair cache: the edited row');
    expect(atomPairs[0], fresh[0], 'the saved AtomPair cache: a row not edited');
  });

  test('options: the defaults are today\'s fingerprints and cache', async () => {
    const col = molTable([ASPIRIN, L_ALANINE, CAFFEINE]);
    const today = (await chemGetFingerprints(col, Fingerprint.Morgan, false)).map(onBits);
    for (const defaults of [{}, {radius: 2, size: 2048, useChirality: false, useFeatures: false}, {radius: null}]) {
      expect((await chemGetFingerprints(col, Fingerprint.Morgan, false, defaults)).map(onBits).join(';'),
        today.join(';'), JSON.stringify(defaults));
    }
    expect(cachedColumns(col), `~smiles.Morgan${engineKey()}`, 'one cache, named as before (saved projects read it)');
    // options a type does not read change nothing for it
    const pattern = (await chemGetFingerprints(col, Fingerprint.Pattern, false)).map(onBits).join(';');
    expect((await chemGetFingerprints(col, Fingerprint.Pattern, false, {radius: 3, useFeatures: true}))
      .map(onBits).join(';'), pattern, 'Pattern reads no radius or features');
  });

  test('options: Morgan radius and size are RDKit\'s in both engines', async () => {
    await inBothEngines(async (engine) => {
      const col = molTable([ASPIRIN, CAFFEINE]);
      const ecfp6 = await chemGetFingerprints(col, Fingerprint.Morgan, false, {radius: 3});
      expect(onBits(ecfp6[0]), RDKIT_BITS.aspirinEcfp6, `${engine}: aspirin, radius 3`);
      const small = await chemGetFingerprints(col, Fingerprint.Morgan, false, {size: 1024});
      expect(small[0]?.length, 1024, `${engine}: the size`);
      expect(onBits(small[0]), RDKIT_BITS.aspirinEcfp4At1024, `${engine}: aspirin, 1024 bits`);
      expect(cachedColumns(col), `~smiles.Morgan(size=1024)${fingerprintsKey(engine)}`,
        `${engine}: the last options' cache (it replaced radius 3's)`);
    });
  });

  test('options: size for every type that reads it, as RDKit', async () => {
    // RDKit mode: a details string Chem's RDKit.js refused would leave every fingerprint null (caught per molecule);
    // aspirin's 1024-bit fingerprints as Chem's RDKit.js 2024.09 gives them (RDKit.js 2026.03.6 and Crux the same):
    // the number of bits set and the first ten
    const expected: {[type: string]: [number, string]} = {
      [Fingerprint.RDKit]: [319, '2,7,9,19,23,24,25,30,35,36'],
      [Fingerprint.Pattern]: [169, '9,13,28,36,37,39,41,42,46,52'],
      [Fingerprint.AtomPair]: [68, '16,32,56,57,68,192,220,221,222,224'],
      [Fingerprint.TopologicalTorsion]: [17, '12,28,48,49,148,376,377,408,576,612'],
    };
    setChemicalEngine(ChemicalEngine.RDKit);
    try {
      for (const [type, [count, first]] of Object.entries(expected)) {
        const [fp] = await chemGetFingerprints(molTable([ASPIRIN]), type as Fingerprint, false, {size: 1024});
        expect(fp?.length, 1024, `${type}: the size`);
        const bits = onBits(fp).split(',');
        expect(`${bits.length}: ${bits.slice(0, 10).join(',')}`, `${count}: ${first}`, `${type}: aspirin's bits`);
      }
    } finally {
      setChemicalEngine(null);
    }
  });

  test('options: each option set has its own cached column', async () => {
    const col = molTable([ASPIRIN, CAFFEINE]);
    const key = engineKey();
    const radius3 = (await chemGetFingerprints(col, Fingerprint.Morgan, false, {radius: 3})).map(onBits);
    const radius2 = (await chemGetFingerprints(col, Fingerprint.Morgan, false)).map(onBits);
    expect(radius3[0], RDKIT_BITS.aspirinEcfp6, 'radius 3');
    expect(radius2[0] !== radius3[0], true, 'radius 2 is not read from radius 3\'s cache');
    expect((await chemGetFingerprints(col, Fingerprint.Morgan, false, {radius: 3})).map(onBits).join(';'),
      radius3.join(';'), 'radius 3 again, from its cache');
    expect(cachedColumns(col), sorted(`~smiles.Morgan${key}`, `~smiles.Morgan(radius=3)${key}`));
    // a column keeps the defaults' cache and the last options asked for
    await chemGetFingerprints(col, Fingerprint.Morgan, false, {radius: 4});
    expect(cachedColumns(col), sorted(`~smiles.Morgan${key}`, `~smiles.Morgan(radius=4)${key}`));
    // an edit reaches the options' cache too
    col.set(1, ASPIRIN);
    const edited = (await chemGetFingerprints(col, Fingerprint.Morgan, false, {radius: 4})).map(onBits);
    expect(edited[1], edited[0], 'the edited row');
  });

  test('options: FCFP is computed by Crux in both engines, as RDKit\'s', async () => {
    // RDKit.js reads no useFeatures (it would give ECFP); Crux gives RDKit's feature-invariant Morgan
    await inBothEngines(async (engine) => {
      const col = molTable([ASPIRIN, CAFFEINE]);
      const fcfp = await chemGetFingerprints(col, Fingerprint.Morgan, false, {useFeatures: true});
      expect(onBits(fcfp[0]), RDKIT_BITS.aspirinFcfp4, `${engine}: aspirin`);
      expect(onBits(fcfp[1]), RDKIT_BITS.caffeineFcfp4, `${engine}: caffeine`);
      expect(cachedColumns(col), '~smiles.Morgan(features).Crux', `${engine}: Crux's cache`);
    });
  });

  test('options: chirality is computed by Crux in both engines, as RDKit 2026.03\'s', async () => {
    // Chem's RDKit.js 2024.09 under the new stereo perception gives L- and D-alanine the same Morgan bits, and its atom
    // pairs and torsions ignore chirality; Crux gives RDKit 2026.03's (it labels the molecule by the new CIP labeller)
    await inBothEngines(async (engine) => {
      const col = molTable([L_ALANINE, D_ALANINE]);
      const chiral = await chemGetFingerprints(col, Fingerprint.Morgan, false, {useChirality: true});
      expect(onBits(chiral[0]), RDKIT_BITS.lAlanineChiral, `${engine}: L-alanine`);
      expect(onBits(chiral[1]), RDKIT_BITS.dAlanineChiral, `${engine}: D-alanine`);
      // atom pairs and torsions: chirality takes part (Chem's RDKit.js ignores it); torsions tell the enantiomers
      // apart, atom pairs do not (as in RDKit: RDKit.js 2026.03.6 and Python give L- and D-alanine the same chiral
      // atom pairs)
      for (const type of [Fingerprint.AtomPair, Fingerprint.TopologicalTorsion]) {
        const fps = (await chemGetFingerprints(col, type, false, {useChirality: true})).map(onBits);
        const flat = (await chemGetFingerprints(col, type, false)).map(onBits);
        expect(fps[0] !== flat[0], true, `${engine}: ${type}: chirality takes part`);
        expect(fps[0] !== fps[1], type === Fingerprint.TopologicalTorsion, `${engine}: ${type}: L- against D-alanine`);
        expect(flat[0], flat[1], `${engine}: ${type} without chirality`);
      }
    });
  });

  test('options: through the getFingerprints function', async () => {
    const col = molTable([ASPIRIN, CAFFEINE]);
    const func = DG.Func.find({package: 'Chem', name: 'getFingerprints'})[0];
    const call = async (params: object) => ((await func.apply({col, _metric: BitArrayMetricsNames.Tanimoto,
      fingerprintType: Fingerprint.Morgan, ...params})).entries as (BitArray | null)[]).map(onBits);
    expect((await call({radius: 3}))[0], RDKIT_BITS.aspirinEcfp6, 'radius 3');
    expect((await call({useFeatures: true}))[1], RDKIT_BITS.caffeineFcfp4, 'FCFP');
    expect((await call({}))[0], (await freshBits([ASPIRIN], Fingerprint.Morgan))[0], 'no options: ECFP4');
  });

  test('fp_withoutDf', async () => {
    const col = DG.Column.fromStrings('smiles', [
      'CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4',
      'Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5',
    ]);
    const res = await chemGetFingerprints(col, Fingerprint.Morgan, false)!;
    expect(res.length, 2);
  }, {benchmark: true});
});
