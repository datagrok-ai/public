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

    console.log("001");
    this.freeTanimotoFingerprints();
    console.log("002");
    if (this._rdKitMols === null) {
      return 0;
    }
    console.log("003");
    this._tanimotoFps = new BitSetFixedArray(this._fpLength, this._rdKitMols.length);
    console.log("004");
    for (let i = 0; i < this._rdKitMols.length; ++i) {
      let item = this._rdKitMols[i];
      try {
        console.log("009");
        const fp = this._rdKitMols[i].get_morgan_fp(this._fpRadius, this._fpLength);
        this._stringFpToArrBits(i, fp, this._tanimotoFps);
        console.log("010");
      } catch (e) {
        console.log("oops");
        return 0;
        // nothing to do
      }
    }
    console.log("005");
    return 1;
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