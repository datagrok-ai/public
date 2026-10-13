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

  test('fp_withoutDf', async () => {
    const col = DG.Column.fromStrings('smiles', [
      'CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4',
      'Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5',
    ]);
    const res = await chemGetFingerprints(col, Fingerprint.Morgan, false)!;
    expect(res.length, 2);
  }, {benchmark: true});
});
