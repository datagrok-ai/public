import {RdKitServiceWorkerBase} from './rdkit-service-worker-base';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {rdKitFingerprintToBitArray, tanimoto} from './chem-common';
import {defaultMorganFpLength, defaultMorganFpRadius} from './chem-common';

export class RdKitServiceWorkerSimilarity extends RdKitServiceWorkerBase {
  _rdKitMols: any[] | null = null;
  _tanimotoFps: BitArray[] | null = null;
  _sample: BitArray;
  readonly _fpLength: number = defaultMorganFpLength;
  readonly _fpRadius: number = defaultMorganFpRadius;

  constructor(module: Object, webRoot: string) {
    super(module, webRoot);
    this._sample = new BitArray(this._fpLength);
  }

  initMorganFingerprints() {
    this.freeMorganFingerprints();
    if (this._rdKitMols === null)
      return;

    this._tanimotoFps = [];
    for (let i = 0; i < this._rdKitMols.length; ++i) {
      const item = this._rdKitMols[i];
      let arr = new BitArray(this._fpLength);
      try {
        const fp = this._rdKitMols[i].get_morgan_fp(this._fpRadius, this._fpLength);
        arr = rdKitFingerprintToBitArray(fp);
      } catch (e) {
        // nothing to do, bit is already 0
      }
      this._tanimotoFps.push(arr);
    }
  }

  getMorganFingerprints() {
    if (this._tanimotoFps === null)
      return [];

    return this._tanimotoFps!.map((e: any) => {
      return {data: e.getRawData().buffer, length: e.length};
    });
  }

  freeMorganFingerprints() {
    this._tanimotoFps = null;
    this._sample = new BitArray(this._fpLength);
  }
}
