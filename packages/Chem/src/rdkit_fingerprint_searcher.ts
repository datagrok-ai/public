import {BitSetFixedArray} from './bitset-fixed-array';

export class RdKitFingerprintSearcher {
  rdKitModule: any | null;
  arr: BitSetFixedArray | null;
  mols: any[] | null;
  readonly fp_length = 64;
  readonly fp_radius = 2;
  sample: BitSetFixedArray | null;

  constructor(module: Object) {
    this.rdKitModule = module;
    this.arr = null;
    this.mols = null;
    this.sample = new BitSetFixedArray(this.fp_length, 1);
  }

  _stringFpToArrBits(item: number, fp: string) {
    for (let j = 0; j < this.fp_length; ++j) {
      if (fp[j] === '1') {
this.arr!.setTrue(item, j);
      } else if (fp[j] === '0') {
this.arr!.setFalse(item, j);
      }
    }
  }

  _stringFpToSampleBits(fp: string) {
    for (let j = 0; j < this.fp_length; ++j) {
      if (fp[j] === '1') {
this.sample!.setTrue(0, j);
      } else if (fp[j] === '0') {
this.sample!.setFalse(0, j);
      }
    }
  }

  init(dict: string[]) {
    this.deinit();
    if (dict.length === 0) {
      this.arr = null;
      return;
    }
    this.arr = new BitSetFixedArray(this.fp_length, dict.length);
    this.mols = [];
    const hashToMolblock: {[_:string] : any} = {};
    const molIdxToHash = [];
    for (let i = 0; i < dict.length; ++i) {
      let item = dict[i];
      let mol;
      try {
        mol = this.rdKitModule.get_mol(item);
        if (item.includes('M  END')) {
          item = mol.normalize_2d_molblock();
          mol.straighten_2d_layout();
          if (!hashToMolblock[item]) {
            hashToMolblock[item] = mol.get_molblock();
          }
        }
      } catch (e) {
        console.error(
          'Possibly a malformed molString: `' + item + '`');
        // preserving indices with a placeholder
        mol = this.rdKitModule.get_mol('');
        // Won't rethrow
      }
      if (mol) {
        const fp = mol.get_morgan_fp(this.fp_radius, this.fp_length);
        this._stringFpToArrBits(i, fp);
        // mol.delete();
      }
      this.mols![i] = mol;
      molIdxToHash.push(item);
    }
    return {molIdxToHash, hashToMolblock};
  }

  search(queryMolString: string, querySmarts: string) {
    const matches: number[] = [];
    if (this.arr) {
      try {
        let queryMol = null;
        try {
          queryMol = this.rdKitModule.get_mol(queryMolString, '{"mergeQueryHs":true}');
        } catch (e2) {
          if (querySmarts !== null && querySmarts !== '') {
            console.log('Cannot parse a MolBlock. Switching to SMARTS');
            queryMol = this.rdKitModule.get_qmol(querySmarts);
          } else {
            throw 'SMARTS not set';
          }
        }
        if (queryMol) {
          if (queryMol.is_valid()) {
            const fp = queryMol.get_morgan_fp(this.fp_radius, this.fp_length);
            this._stringFpToSampleBits(fp);
            for (let i = 0; i < this.arr!.length; ++i) {
              if (this.arr!.entails(i, this.sample!)) {
                if (this.mols![i]!.get_substruct_match(queryMol) !== '{}') {
                  matches.push(i);
                }
              }
            }
          }
          queryMol.delete();
        }
      } catch (e) {
        console.error(
          'Possibly a malformed query: `' + queryMolString + '`');
        // Won't rethrow
      }
    }
    return '[' + matches.join(', ') + ']';
  }

  deinit() {
    this.arr = null;
  }
}
