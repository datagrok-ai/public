import { RdKitServiceWorkerBase } from './rdkit-service-worker-base';
import { defaultMorganFpLength, defaultMorganFpRadius, Fingerprint } from '../utils/chem-common';
import { RDModule } from '@datagrok-libraries/chem-meta/src/rdkit-api';
import { getMolSafe } from '../utils/mol-creation_rdkit';

export class RdKitServiceWorkerSimilarity extends RdKitServiceWorkerBase {
  readonly _fpLength: number = defaultMorganFpLength;
  readonly _fpRadius: number = defaultMorganFpRadius;

  constructor(module: RDModule, webRoot: string) {
    super(module, webRoot);
  }

  getFingerprints(fingerprintType: Fingerprint) {
    if (this._rdKitMols === null)
      return [];

    const fps: (Uint8Array | null)[] = [];
    try {
      switch (fingerprintType) {
        case Fingerprint.Pattern:
          for (let i = 0; i < this._rdKitMols.length; ++i) {
            try {
              if (!this._rdKitMols[i])
                fps.push(null);
              else
                fps.push(this._rdKitMols[i]!.get_pattern_fp_as_uint8array());
            } catch {
              fps.push(null);
              continue;
            }
          }
          break;
        case Fingerprint.Morgan:
          for (let i = 0; i < this._rdKitMols.length; ++i) {
            try {
              if (!this._rdKitMols[i] || this._rdKitMols[i] && this._rdKitMols[i]?.is_qmol)
                fps.push(null);
              else
                fps.push(this._rdKitMols[i]!.get_morgan_fp_as_uint8array(JSON.stringify({
                  radius: this._fpRadius,
                  nBits: this._fpLength,
                })));
            } catch (error) {
              fps.push(null);
              continue;
            }
          }
          break;
        default:
          throw Error('Unknown fingerprint type: ' + fingerprintType);
      }
    } catch (e) {
      // nothing to do, bit is already 0
    }
    return fps!.map((el: any) => {
      return { data: el, length: el ? el.length : 0 };
    });
  }


  getFingerprintsWithMolsOnFly(dict: string[], fingerprintType: Fingerprint) {
    const fps: (Uint8Array | null)[] = [];
    for (let i = 0; i < dict.length; ++i) {
      const item = dict[i];
      let mol;
      if (!item || item === '')
        fps.push(new Uint8Array());
      else {
        const molSafe = getMolSafe(item, {}, this._rdKitModule);
        mol = molSafe.mol;
        if (mol === null) {
          fps.push(null);
        } else {
          switch (fingerprintType) {
            case Fingerprint.Pattern:
              try {
                fps.push(mol!.get_pattern_fp_as_uint8array());
              } catch {
                fps.push(null);
              }
              break;
            case Fingerprint.Morgan:
              try {
                if (mol.is_qmol)
                  fps.push(null);
                else
                  fps.push(mol.get_morgan_fp_as_uint8array(JSON.stringify({
                    radius: this._fpRadius,
                    nBits: this._fpLength,
                  })));
              } catch (error) {
                fps.push(null);
              }
              break;
            default:
              mol.delete();
              throw Error('Unknown fingerprint type: ' + fingerprintType);
          }
          mol.delete();
        }
      }
    }
    return fps!.map((el: any) => {
      return { data: el, length: el ? el.length : 0 };
    });
  }
}
