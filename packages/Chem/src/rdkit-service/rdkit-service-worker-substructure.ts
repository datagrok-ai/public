import {RdKitServiceWorkerSimilarity} from './rdkit-service-worker-similarity';
import {isMolBlock} from '../utils/chem-utils';
import {RDModule, RDMol} from '@datagrok-libraries/bio/src/rdkit-api';

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

  getMol(molString: string) {
    return this._rdKitModule.get_mol(molString, '{"mergeQueryHs":true}');
  }

  getQMol(molString: string) {
    return this._rdKitModule.get_qmol(molString);
  }

  searchSubstructure(queryMolString: string, queryMolBlockFailover: string, bitset?: boolean[]): string {
    const matches: number[] = [];
    if (this._rdKitMols) {
      let queryMol: RDMol | null = null;
      try {
        if (isMolBlock(queryMolString))
          queryMol = this.getMol(queryMolString);
        else {
          let createdQMol = false;
          try {
            queryMol = this.getQMol(queryMolString);
            createdQMol = queryMol.is_valid();
          } catch (e) {
            queryMol?.delete();
            queryMol = null;
          }
          if (createdQMol) {
            let mol = null;
            let createdMol = false;
            try {
              mol = this.getMol(queryMolString);
              createdMol = mol.is_valid();
            } catch (e) {
              mol?.delete();
            }
            if (createdMol) { // check the qmol is proper
              const match = mol!.get_substruct_match(queryMol!);
              if (match === '{}') {
                queryMol!.delete();
                queryMol = mol;
              }
            } // else, this looks to be a real SMARTS
          } else { // failover to queryMolBlockFailover
            queryMol = this.getMol(queryMolBlockFailover);
          }
        }
        if (queryMol && queryMol.is_valid()) {
            if (bitset) {
              for (let i = 0; i < bitset.length; ++i) {
                if (bitset[i] && this._rdKitMols[i]!.get_substruct_match(queryMol) !== '{}') // Is patternFP iff?
                  matches.push(i);
              }
            } else {
              for (let i = 0; i < this._rdKitMols!.length; ++i)
                if (this._rdKitMols[i]!.get_substruct_match(queryMol) !== '{}')
                  matches.push(i);
            }
          } else
          throw new Error('Chem | Search pattern cannot be set');
      } catch (e) {
        console.error(
          'Chem | Possibly a malformed query: `' + queryMolString + '`');
        // Won't rethrow
      } finally {
        queryMol?.delete();
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
