import {RdKitServiceWorkerBase} from './rdkit_service_worker_base';
import {BitSetFixedArray} from "./bitset-fixed-array";

export class RdkitServiceSubstructure extends RdKitServiceWorkerBase {

  _rdKitMols: any[] | null = null;

  constructor(module: Object, webRoot: string) {
    super(module, webRoot);
  }

  initMoleculesStructures(dict: string[]) {
    this.freeMoleculesStructures();
    if (dict.length === 0) {
      return;
    }
    this._rdKitMols = [];
    let hashToMolblock: {[_:string] : any} = {};
    let molIdxToHash = [];
    for (let i = 0; i < dict.length; ++i) {
      let item = dict[i];
      let mol;
      try {
        mol = this._rdKitModule.get_mol(item);
        if (item.includes('M  END')) {
          item = mol.normalize_2d_molblock();
          mol.straighten_2d_layout();
          if (!hashToMolblock[item]) {
            hashToMolblock[item] = mol.get_molblock();
          }
        }
      } catch (e) {
        console.error(
          "Possibly a malformed molString: `" + item + "`");
        // preserving indices with a placeholder
        mol = this._rdKitModule.get_mol('');
        // Won't rethrow
      }
      this._rdKitMols!.push(mol);
      molIdxToHash.push(item);
    }
    return { molIdxToHash, hashToMolblock };
  }

  searchSubstructure(queryMolString: string, querySmarts: string) {
    let matches: number[] = [];
    if (this._rdKitMols) {
      try {
        let queryMol = null;
        try {
          queryMol = this._rdKitModule.get_mol(queryMolString, "{\"mergeQueryHs\":true}");
        } catch (e2) {
          if (querySmarts !== null && querySmarts !== '') {
            console.log("Cannot parse a MolBlock. Switching to SMARTS");
            queryMol = this._rdKitModule.get_qmol(querySmarts);
          } else {
            throw 'SMARTS not set';
          }
        }
        if (queryMol) {
          if (queryMol.is_valid()) {
            for (let i = 0; i < this._rdKitMols!.length; ++i)
              if (this._rdKitMols![i]!.get_substruct_match(queryMol) !== "{}")
                  matches.push(i);
          }
          queryMol.delete();
        }
      } catch (e) {
        console.error(
          "Possibly a malformed query: `" + queryMolString + "`");
        // Won't rethrow
      }
    }
    return '[' + matches.join(', ') + ']';
  }

  freeMoleculesStructures() {
    if (this._rdKitMols !== null) {
      for (let mol of this._rdKitMols!) {
        mol.delete();
        mol = null;
      }
      this._rdKitMols = null;
    }
  }

}