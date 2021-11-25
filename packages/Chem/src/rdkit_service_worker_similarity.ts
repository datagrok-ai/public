import {RdKitServiceWorkerBase} from './rdkit_service_worker_base';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {rdKitFingerprintToBitArray, tanimoto} from './chem_common';

/*
function rdKitFingerprintToBitArray(fp: string, fpLength: number) {
  let arr = new BitArray(fpLength);
  for (let j = 0; j < fpLength; ++j) {
    if (fp[j] === '1')
      arr.setTrue(j);
  }
  return arr;
}
 */

/*
function tanimoto(x: BitArray, y: BitArray) {
  const total = x.trueCount() + y.trueCount();
  if (total == 0)
    return 1.0;
  const common = x.andWithCountBits(y, true);
  return common / (total - common);
}

 */

export class RdKitServiceWorkerSimilarity extends RdKitServiceWorkerBase {

  _rdKitMols: any[] | null = null;
  _tanimotoFps: BitArray[] | null = null;
  _sample: BitArray;
  readonly _fpLength: number = 128;
  readonly _fpRadius: number = 2;

  constructor(module: Object, webRoot: string) {
    super(module, webRoot);
    this._sample = new BitArray(this._fpLength);
  }

  initMorganFingerprints() {
    this.freeMorganFingerprints();
    if (this._rdKitMols === null) {
      return;
    }
    this._tanimotoFps = [];
    for (let i = 0; i < this._rdKitMols.length; ++i) {
      let item = this._rdKitMols[i];
      let arr = new BitArray(this._fpLength);
      try {
        const fp = this._rdKitMols[i].get_morgan_fp(this._fpRadius, this._fpLength);
        arr = rdKitFingerprintToBitArray(fp, this._fpLength);
      } catch (e) {
        // nothing to do
      }
      this._tanimotoFps.push(arr);
    }
  }

  getSimilarities(sampleMol: string) {
    let distances = new Array(this._rdKitMols!.length).fill(0.0);
    try {
      const mol = this._rdKitModule.get_mol(sampleMol);
      const fp = mol.get_morgan_fp(this._fpRadius, this._fpLength);
      const sample = rdKitFingerprintToBitArray(fp, this._fpLength);
      for (let i = 0; i < this._rdKitMols!.length; ++i) {
        distances[i] = tanimoto(this._tanimotoFps![i], sample);
      }
    } finally {
      return distances;
    }
  }

  getMorganFingerprints() {
    return this._tanimotoFps!.map((e: any) => {return {data: e.getRawData().buffer, length: e.length}});
  }

  freeMorganFingerprints() {
    this._tanimotoFps = null;
    this._sample = new BitArray(this._fpLength);
  }

}