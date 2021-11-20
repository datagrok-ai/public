import {RdkitServiceSubstructure} from './rdkit_service_worker_substructure';
import {BitSetFixedArray} from "./bitset-fixed-array";
import {BitSet} from "datagrok-api/src/dataframe";

export class RdKitServiceWorkerSimilarity extends RdkitServiceSubstructure {

  _tanimotoFps: BitSetFixedArray | null = null;
  _sample: BitSetFixedArray;
  readonly _fpLength: number = 128;
  readonly _fpRadius: number = 2;

  constructor(module: Object, webRoot: string) {
    super(module, webRoot);
    this._sample = new BitSetFixedArray(this._fpLength, 1);
  }

  _stringFpToArrBits(item: number, fp: string, arr: BitSetFixedArray) {
    for (let j = 0; j < this._fpLength; ++j) {
      if (fp[j] === '1')
        arr.setTrue(item, j);
      else if (fp[j] === '0')
        arr.setFalse(item, j);
    }
  }

  initTanimotoFingerprints(/* the structures are already passed */) {
    this.freeTanimotoFingerprints();
    if (this._rdKitMols === null) {
      return;
    }
    this._tanimotoFps = new BitSetFixedArray(this._fpLength, this._rdKitMols.length);
    for (let i = 0; i < this._rdKitMols.length; ++i) {
      let item = this._rdKitMols[i];
      try {
        const fp = this._rdKitMols[i].get_morgan_fp(this._fpRadius, this._fpLength);
        this._stringFpToArrBits(i, fp, this._tanimotoFps);
      } catch (e) {
        // nothing to do
      }
    }
  }

  _tanimoto(arr: BitSetFixedArray, i: number, sample: BitSetFixedArray): number {
    const total = arr.count(i) + sample.count(0);
    if (total == 0)
      return 1.0;
    const common = arr.andWithCountBits(i, sample);
    return common / (total - common);
  }

  getSimilarities(sampleMol: string) {
    let distances = new Array(this._rdKitMols!.length).fill(0.0);
    try {
      const mol = this._rdKitModule.get_mol(sampleMol);
      const fp = mol.get_morgan_fp(this._fpRadius, this._fpLength);
      this._stringFpToArrBits(0, fp, this._sample);
      for (let i = 0; i < this._rdKitMols!.length; ++i) {
        distances[i] = this._tanimoto(this._tanimotoFps!, i, this._sample);
      }
    } finally {
      return distances;
    }
  }

  freeTanimotoFingerprints() {
    this._tanimotoFps = null;
    this._sample = new BitSetFixedArray(this._fpLength, 1);
  }

}