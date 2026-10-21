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

const ASPIRIN = 'CC(=O)Oc1ccccc1C(=O)O';
const CAFFEINE = 'Cn1cnc2c1c(=O)n(C)c(=O)n2C';

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

/** The fingerprints of molecules in a fresh table (no caches). */
async function freshBits(smiles: string[], type: Fingerprint): Promise<string[]> {
  return (await chemGetFingerprints(molTable(smiles), type, false)).map(onBits);
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

  test('fp_withoutDf', async () => {
    const col = DG.Column.fromStrings('smiles', [
      'CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4',
      'Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5',
    ]);
    const res = await chemGetFingerprints(col, Fingerprint.Morgan, false)!;
    expect(res.length, 2);
  }, {benchmark: true});
});
