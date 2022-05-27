import {RdKitServiceWorkerBase} from './rdkit-service-worker-base';
import {defaultMorganFpLength, defaultMorganFpRadius, Fingerprint} from '../utils/chem-common';
import {RDModule} from "../rdkit-api";

export class RdKitServiceWorkerSimilarity extends RdKitServiceWorkerBase {
  readonly _fpLength: number = defaultMorganFpLength;
  readonly _fpRadius: number = defaultMorganFpRadius;

  constructor(module: RDModule, webRoot: string) {
    super(module, webRoot);
  }

  getFingerprints(fingerprintType: Fingerprint) {
    if (this._rdKitMols === null)
      return;

    const fps: Uint8Array[] = [];
    let fp: Uint8Array;
    for (let i = 0; i < this._rdKitMols.length; ++i) {
      try {
        // try map maybe?
        switch (fingerprintType) {
        case Fingerprint.Morgan:
          fp = this._rdKitMols[i].get_morgan_fp_as_uint8array(this._fpRadius, this._fpLength);
          break;
        case Fingerprint.Pattern:
          fp = this._rdKitMols[i].get_pattern_fp_as_uint8array();
          break;
        default:
          throw Error('Unknown fingerprint type: ' + fingerprintType);
        }
        fps.push(fp);
      } catch (e) {
        // nothing to do, bit is already 0
      }
    }
    return fps!.map((e: any) => {
      return {data: e, length: e.length};
    });
  }
}
