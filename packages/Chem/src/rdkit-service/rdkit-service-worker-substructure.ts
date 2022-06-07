import {RdKitServiceWorkerSimilarity} from './rdkit-service-worker-similarity';
import {isMolBlock} from '../utils/chem-utils';
import {RDModule} from "../rdkit-api";

interface InitMoleculesStructuresResult {
  molIdxToHash: string[];
  hashToMolblock: {
      [_: string]: any;
  }
}

export class RdKitServiceWorkerSubstructure extends RdKitServiceWorkerSimilarity {
  readonly _patternFpLength = 2048;
  readonly _patternFpUint8Length = 256;

  constructor(module: RDModule, webRoot: string) {
    super(module, webRoot);
  }

  initMoleculesStructures(dict: string[], normalizeCoordinates: boolean): InitMoleculesStructuresResult {
    this.freeMoleculesStructures();
    if (dict.length === 0)
      return {molIdxToHash: [], hashToMolblock: {}};
    this._rdKitMols = [];

    const hashToMolblock: {[_:string] : any} = {};
    const molIdxToHash = [];
    for (let i = 0; i < dict.length; ++i) {
      let item = dict[i];
      let mol = null;

      try {
        mol = this._rdKitModule.get_mol(item);

        if (normalizeCoordinates) {
          if (isMolBlock(item)) {
            mol.normalize_depiction();
            item = mol.compute_hash();
            mol.straighten_depiction();
            if (!hashToMolblock[item])
              hashToMolblock[item] = mol.get_molblock();
          }
        }
      } catch (e) {
        console.error('Chem | Possibly a malformed molString: `' + item + '`');
        mol?.delete();
        mol = this._rdKitModule.get_mol('');
      }
      this._rdKitMols.push(mol);
      if (normalizeCoordinates)
        molIdxToHash.push(item);
    }
    return {molIdxToHash, hashToMolblock};
  }

  searchSubstructure(queryMolString: string, querySmarts: string, patternFps?: Uint8Array[]): string {
    const matches: number[] = [];
    if (this._rdKitMols) {
      try {
        let queryMol = null;
        try {
          try {
            queryMol = this._rdKitModule.get_mol(queryMolString, '{"mergeQueryHs":true}');
          } catch (e2) {
            queryMol?.delete();
            queryMol = null;
            if (querySmarts !== null && querySmarts !== '') {
              console.log('Chem | Cannot parse a MolBlock. Switching to SMARTS');
              queryMol = this._rdKitModule.get_qmol(querySmarts);
            } else
              throw new Error('Chem | SMARTS not set');
          }
          if (queryMol && queryMol.is_valid()) {
            if (patternFps) {
              const fpRdKit = queryMol.get_pattern_fp_as_uint8array(this._patternFpLength);
              checkEl:
              for (let i = 0; i < patternFps.length; ++i) {
                for (let j = 0; j < this._patternFpUint8Length; ++j)
                  if ((patternFps[i][j] & fpRdKit[j]) != fpRdKit[j])
                    continue checkEl;

                if (this._rdKitMols[i]!.get_substruct_match(queryMol) !== '{}') // Is patternFP iff?
                  matches.push(i);
              }
            } else {
              for (let i = 0; i < this._rdKitMols!.length; ++i)
                if (this._rdKitMols[i]!.get_substruct_match(queryMol) !== '{}')
                  matches.push(i);
            }
          }
        } finally {
          queryMol?.delete();
        }
      } catch (e) {
        console.error(
          e, 'Possibly a malformed query: `' + queryMolString + '`');
        // Won't rethrow
      }
    }
    return '[' + matches.join(', ') + ']';
  }

  freeMoleculesStructures(): void {
    if (this._rdKitMols !== null) {
      for (let mol of this._rdKitMols!)
        mol.delete();
      this._rdKitMols = null;
    }
  }
}
