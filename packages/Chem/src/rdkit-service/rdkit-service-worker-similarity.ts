import {RdKitServiceWorkerBase} from './rdkit-service-worker-base';
import {defaultMorganFpLength, defaultMorganFpRadius, Fingerprint} from '../utils/chem-common';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

export class RdKitServiceWorkerSimilarity extends RdKitServiceWorkerBase {
  readonly _fpLength: number = defaultMorganFpLength;
  readonly _fpRadius: number = defaultMorganFpRadius;

  constructor(module: RDModule, webRoot: string) {
    super(module, webRoot);
  }

  getFingerprints(fingerprintType: Fingerprint) {
    if (this._rdKitMols === null)
      return [];

    const fps: Uint8Array[] = [];
    try {
      switch (fingerprintType) {
      case Fingerprint.Pattern:
        for (let i = 0; i < this._rdKitMols.length; ++i)
          fps.push(this._rdKitMols[i].get_pattern_fp_as_uint8array());
        break;
      case Fingerprint.Morgan:
        for (let i = 0; i < this._rdKitMols.length; ++i) {
          fps.push(this._rdKitMols[i].get_morgan_fp_as_uint8array(JSON.stringify({
            radius: this._fpRadius,
            nBits: this._fpLength,
          })));
        }
        break;
      default:
        throw Error('Unknown fingerprint type: ' + fingerprintType);
      }
    } catch (e) {
      // nothing to do, bit is already 0
    }
    return fps!.map((el: any) => {
      return {data: el, length: el.length};
    });
  }
}
