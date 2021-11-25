import {RdkitServiceSubstructure} from './rdkit_service_worker_substructure';
// import {BitSetFixedArray} from "./bitset-fixed-array";
import BitArray from '@datagrok-libraries/utils/src/bit-array';
export class RdKitServiceWorkerSimilarity extends RdkitServiceSubstructure {

  _tanimotoFps: BitArray[] | null = null;
  _sample: BitArray;
  readonly _fpLength: number = 128;
  readonly _fpRadius: number = 2;

  constructor(module: Object, webRoot: string) {
    super(module, webRoot);
    this._sample = new BitArray(this._fpLength);
  }

  initMorganFingerprints(/* the structures are already passed */) {
    this.freeMorganFingerprints();
    if (this._rdKitMols === null) {
      return;
    }
    this._tanimotoFps = []; // new  new BitArray(this._fpLength);
    for (let i = 0; i < this._rdKitMols.length; ++i) {
      let item = this._rdKitMols[i];
      let arr = new BitArray(this._fpLength);
      try {
        const fp = this._rdKitMols[i].get_morgan_fp(this._fpRadius, this._fpLength);
        arr = this._stringFpToArrBits(fp, arr, this._fpLength);
      } catch (e) {
        // nothing to do
      }
      this._tanimotoFps.push(arr);
    }
  }

  _tanimoto(arr: BitArray[], i: number, sample: BitArray): number {
    const total = arr[i].trueCount() + sample.trueCount();
    if (total == 0)
      return 1.0;
    const common = arr[i].andWithCountBits(sample, true);
    return common / (total - common);
  }

  getSimilarities(sampleMol: string) {
    let distances = new Array(this._rdKitMols!.length).fill(0.0);
    try {
      const mol = this._rdKitModule.get_mol(sampleMol);
      const fp = mol.get_morgan_fp(this._fpRadius, this._fpLength);
      this._sample = this._stringFpToArrBits(fp, this._sample, this._fpLength);
      for (let i = 0; i < this._rdKitMols!.length; ++i) {
        distances[i] = this._tanimoto(this._tanimotoFps!, i, this._sample);
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