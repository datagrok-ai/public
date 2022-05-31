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

    let fps: Uint8Array[];
      try {
        switch (fingerprintType) {
          case Fingerprint.Morgan:
            fps = this._rdKitMols.map(el => el.get_morgan_fp_as_uint8array(this._fpRadius, this._fpLength));
            break;
          case Fingerprint.Pattern:
            fps = this._rdKitMols.map(el => el.get_pattern_fp_as_uint8array());
            break;
        default:
          throw Error('Unknown fingerprint type: ' + fingerprintType);
        }
    } catch (e) {
      // nothing to do, bit is already 0
    }
    return fps!.map((e: any) => {
      return {data: e, length: e.length};
    });
  }
}
