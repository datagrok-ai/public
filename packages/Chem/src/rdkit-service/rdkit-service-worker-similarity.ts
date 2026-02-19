import {RdKitServiceWorkerBase} from './rdkit-service-worker-base';
import {defaultMorganFpLength, defaultMorganFpRadius, Fingerprint} from '../utils/chem-common';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getMolSafe, IMolContext} from '../utils/mol-creation_rdkit';

export interface IFpResult{
  fps: Array<Uint8Array | null>;
  smiles: Array<string | null> | null;
}

export class RdKitServiceWorkerSimilarity extends RdKitServiceWorkerBase {
  readonly _fpLength: number = defaultMorganFpLength;
  readonly _fpRadius: number = defaultMorganFpRadius;
  constructor(module: RDModule, webRoot: string) {
    super(module, webRoot);
  }

  /**
   * Calculates fingerprints either on pre-created array of RDMols or creating RDMOls on the fly.
   * If you want to use pre-created array of RDMols you should first create it by using initMoleculesStructures
   * web-worker method.
   *
   * @param {Fingerprint} fingerprintType Type of Fingerprint
   * @param {string[]} molecules List of molecule strings to calculate fingerprints on. If passed,
   * fps mols are created on the fly
   * @param {boolean} getCanonicalSmiles If passed canonical smiles are also calculated and returned
   * In case it is passed to function RDMols will be created on the fly.
   */

  async getFingerprints(fingerprintType: Fingerprint, molecules?: string[], getCanonicalSmiles?: boolean): Promise<IFpResult> {
    if (!molecules || this._requestTerminated)
      return {fps: [], smiles: null};
    let addedToCache = false;
    const fpLength = molecules.length;
    const fps = new Array<Uint8Array | null>(fpLength).fill(null);
    const morganFpParams = fingerprintType === Fingerprint.Morgan ?
      JSON.stringify({radius: this._fpRadius, nBits: this._fpLength}) : null;
    const canonicalSmilesArr = getCanonicalSmiles ? new Array<string | null>(fpLength).fill(null) : null;
    for (let i = 0; i < fpLength; ++i) {

      if (i % this._terminationCheckDelay === 0) //every N molecules check for termination flag
        await new Promise((r) => setTimeout(r, 0));
      if (this._requestTerminated)
        return { fps: fps, smiles: canonicalSmilesArr };

      const item = molecules[i];
      if (!item || item === '')
        continue;
      let rdMol = this._molsCache?.get(molecules[i]);
      if (!rdMol) {
        const mol: IMolContext = getMolSafe(item, {}, this._rdKitModule);
        rdMol = mol?.mol;
        if (rdMol)
          rdMol.is_qmol = mol?.isQMol;
      }
      if (rdMol) {
        try {
          if (canonicalSmilesArr)
            canonicalSmilesArr[i] = rdMol.get_smiles();
          switch (fingerprintType) {
            case Fingerprint.Pattern:
              fps[i] = rdMol.get_pattern_fp_as_uint8array();
              break;
            case Fingerprint.Morgan:
              if (!rdMol.is_qmol) {
                rdMol.remove_hs_in_place(); // hydrogens can cause identical molecules to have different fingerprints
                fps[i] = rdMol.get_morgan_fp_as_uint8array(morganFpParams!);
              }
              break;
            case Fingerprint.AtomPair:
              fps[i] = rdMol.get_atom_pair_fp_as_uint8array();
              break;
            case Fingerprint.MACCS:
              fps[i] = rdMol.get_maccs_fp_as_uint8array();
              break;
            case Fingerprint.RDKit:
              fps[i] = rdMol.get_rdkit_fp_as_uint8array();
              break;
            case Fingerprint.TopologicalTorsion:
              fps[i] = rdMol.get_topological_torsion_fp_as_uint8array();
              break;
            default:
              rdMol?.delete();
              throw Error('Unknown fingerprint type: ' + fingerprintType);
          }
          //only first N (MAX_MOL_CACHE_SIZE) molecules of the dataset are added to cache
          addedToCache = this.addToCache(rdMol);
        } catch {
          // nothing to do, fp is already null
        } finally {
          if (!addedToCache) { //do not delete mol in case it is in cache
            rdMol?.delete();
          }
        }
      }
    }
    return {fps: fps, smiles: canonicalSmilesArr};
  }
}
