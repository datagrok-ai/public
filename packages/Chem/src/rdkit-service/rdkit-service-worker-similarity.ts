import { RdKitServiceWorkerBase } from './rdkit-service-worker-base';
import { defaultMorganFpLength, defaultMorganFpRadius, Fingerprint } from '../utils/chem-common';
import { RDModule } from '@datagrok-libraries/chem-meta/src/rdkit-api';
import { getMolSafe } from '../utils/mol-creation_rdkit';

export interface IFingerprint {
  data: Uint8Array;
  length: number;
}

export class RdKitServiceWorkerSimilarity extends RdKitServiceWorkerBase {
  readonly _fpLength: number = defaultMorganFpLength;
  readonly _fpRadius: number = defaultMorganFpRadius;

  constructor(module: RDModule, webRoot: string) {
    super(module, webRoot);
  }

  getFingerprints(fingerprintType: Fingerprint, dict?: string[]): IFingerprint[] {
    if (this._rdKitMols === null && !dict)
      return [];

    const fpLength = dict ? dict.length : this._rdKitMols!.length;
    const fps = new Array<Uint8Array | null>(fpLength).fill(null);
    for (let i = 0; i < fpLength; ++i) {
      let mol;
      if (dict) {
        const item = dict[i];
        if (!item || item === '') {
          fps[i] = new Uint8Array();
          continue;
        }
        mol = getMolSafe(item, {}, this._rdKitModule);
      }
      const rdMol = dict ? mol?.mol : this._rdKitMols![i];
      const isQMol = dict ? mol?.isQMol : this._rdKitMols![i]?.is_qmol;
      try {
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
                fps[i] = rdMol.get_morgan_fp_as_uint8array(JSON.stringify({
                  radius: this._fpRadius,
                  nBits: this._fpLength,
                }));
            } catch (error) {
              //do nothing, fp is already null
            }
            break;
          default:
            if (dict)
              mol?.mol?.delete();
            throw Error('Unknown fingerprint type: ' + fingerprintType);
        }
      } catch {
        // nothing to do, fp is already null
      }
      if (dict)
        mol?.mol?.delete();
    }
    return fps!.map((el: any) => {
      return { data: el, length: el ? el.length : 0 };
    });
  }

}
