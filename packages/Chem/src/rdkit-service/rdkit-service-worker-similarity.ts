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
   * @param {string[]} molecules List of molecule strings to calculate fingerprints on. If passed, fps mols are created on the fly
   * @param {boolean} getCanonicalSmiles If passed canonical smiles are also calculated and returned
   * In case it is passed to function RDMols will be created on the fly.
   */

  getFingerprints(fingerprintType: Fingerprint, molecules?: string[], getCanonicalSmiles?: boolean): IFpResult {
    if (this._rdKitMols === null && !molecules)
      return {fps: [], smiles: null};

    const fpLength = molecules ? molecules.length : this._rdKitMols!.length;
    const fps = new Array<Uint8Array | null>(fpLength).fill(null);
    const morganFpParams = fingerprintType === Fingerprint.Morgan ?
      JSON.stringify({radius: this._fpRadius, nBits: this._fpLength}) : null;
    const canonicalSmiles = getCanonicalSmiles ? new Array<string | null>(fpLength).fill(null) : null;
    for (let i = 0; i < fpLength; ++i) {
      let mol: IMolContext | null = null;
      if (molecules) {
        const item = molecules[i];
        if (item && item !== '')
          mol = getMolSafe(item, {}, this._rdKitModule);
      }
      const rdMol = molecules ? mol?.mol : this._rdKitMols![i];
      const isQMol = molecules ? mol?.isQMol : this._rdKitMols![i]?.is_qmol;
      try {
        if (canonicalSmiles && rdMol)
          canonicalSmiles[i] = rdMol.get_smiles();
        switch (fingerprintType) {
        case Fingerprint.Pattern:
          try {
            if (rdMol)
              fps[i] = rdMol.get_pattern_fp_as_uint8array();
          } catch {
            //do nothing, fp is already null
          }
          break;
        case Fingerprint.Morgan:
          try {
            if (rdMol && !isQMol)
              fps[i] = rdMol.get_morgan_fp_as_uint8array(morganFpParams!);
          } catch (error) {
            //do nothing, fp is already null
          }
          break;
        default:
          if (molecules)
            mol?.mol?.delete();
          throw Error('Unknown fingerprint type: ' + fingerprintType);
        }
      } catch {
        // nothing to do, fp is already null
      } finally {
      if (molecules)
        mol?.mol?.delete();
      }
    }
    return {fps: fps, smiles: canonicalSmiles};
  }
}
