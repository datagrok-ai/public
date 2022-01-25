import {RdKitServiceWorkerBase} from './rdkit-service-worker-base';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {rdKitFingerprintToBitArray} from '../utils/chem-common';
import {defaultMorganFpLength, defaultMorganFpRadius} from '../utils/chem-common';
import {RDModule} from "../rdkit-api";

export class RdKitServiceWorkerSimilarity extends RdKitServiceWorkerBase {
  _morganFps: BitArray[] | null = null;
  readonly _fpLength: number = defaultMorganFpLength;
  readonly _fpRadius: number = defaultMorganFpRadius;

  constructor(module: RDModule, webRoot: string) {
    super(module, webRoot);
  }

  initMorganFingerprints() {
    this.freeMorganFingerprints();
    if (this._rdKitMols === null)
      return;

    this._morganFps = [];
    for (let i = 0; i < this._rdKitMols.length; ++i) {
      const item = this._rdKitMols[i];
      let arr = new BitArray(this._fpLength);
      try {
        const fp = this._rdKitMols[i].get_morgan_fp_as_uint8array(this._fpRadius, this._fpLength);
        arr = rdKitFingerprintToBitArray(fp);
      } catch (e) {
        // nothing to do, bit is already 0
      }
      this._morganFps.push(arr);
    }
  }

  getMorganFingerprints() {
    if (this._morganFps === null)
      return [];

    return this._morganFps!.map((e: any) => {
      return {data: e.getRawData().buffer, length: e.length};
    });
  }

  freeMorganFingerprints() {
    this._morganFps = null;
  }
}
