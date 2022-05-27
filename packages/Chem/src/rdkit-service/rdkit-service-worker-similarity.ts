import {RdKitServiceWorkerBase} from './rdkit-service-worker-base';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {rdKitFingerprintToBitArray} from '../utils/chem-common';
import {defaultMorganFpLength, defaultMorganFpRadius} from '../utils/chem-common';
import {RDModule} from "../rdkit-api";

export class RdKitServiceWorkerSimilarity extends RdKitServiceWorkerBase {
  readonly _fpLength: number = defaultMorganFpLength;
  readonly _fpRadius: number = defaultMorganFpRadius;

  constructor(module: RDModule, webRoot: string) {
    super(module, webRoot);
  }

  getMorganFingerprints() {
    if (this._rdKitMols === null)
      return;

    const morganFps: Uint8Array[] = [];
    for (let i = 0; i < this._rdKitMols.length; ++i) {
      // let fp = new BitArray(this._fpLength);
      try {
        const fp = this._rdKitMols[i].get_morgan_fp_as_uint8array(this._fpRadius, this._fpLength);
        morganFps.push(fp);
        // arr = rdKitFingerprintToBitArray(fp);
      } catch (e) {
        // nothing to do, bit is already 0
      }
    }

    return morganFps!.map((e: any) => {
      return {data: e, length: e.length};
    });
  }

  getPatternFingerprints() {
    if (this._rdKitMols === null)
      return;

    const patternFps: Uint8Array[] = [];
    for (let i = 0; i < this._rdKitMols.length; ++i) {
      // let fp = new BitArray(this._fpLength);
      try {
        const fp = this._rdKitMols[i].get_pattern_fp_as_uint8array();
        patternFps.push(fp);
        // arr = rdKitFingerprintToBitArray(fp);
      } catch (e) {
        // nothing to do, bit is already 0
      }
    }

    return patternFps!.map((e: any) => {
      return {data: e, length: e.length};
    });
  }
}
