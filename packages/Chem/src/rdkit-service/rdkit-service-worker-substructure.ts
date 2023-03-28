import {RdKitServiceWorkerSimilarity} from './rdkit-service-worker-similarity';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {syncQueryAromatics} from '../utils/aromatic-utils';
import {getMolSafe} from '../utils/mol-creation_rdkit';
import {isMolBlock} from '../utils/chem-common';

export enum MolNotation {
  Smiles = 'smiles',
  Smarts = 'smarts',
  MolBlock = 'molblock', // molblock V2000
  V3KMolBlock = 'v3Kmolblock', // molblock V3000
  Unknown = 'unknown',
}

const MALFORMED_MOL_V2000 = `
Malformed

  0  0  0  0  0  0  0  0  0  0999 V2000
M  END`;
const MALFORMED_MOL_V3000 = `
Malformed

  0  0  0  0  0  0            999 V3000
M  END`;

function validateMol(mol: RDMol | null, molString: string) : void {
  if (mol === null)
    throw new Error('FATAL RDKit Error: Created a null molecule with no exception ' + molString);
  if (!mol.is_valid())
    throw new Error('FATAL RDKit Error: Created a not valid molecule with no exception ' + molString);
}

export class RdKitServiceWorkerSubstructure extends RdKitServiceWorkerSimilarity {
  readonly _patternFpLength = 2048;
  readonly _patternFpUint8Length = 256;

  constructor(module: RDModule, webRoot: string) {
    super(module, webRoot);
  }

  initMoleculesStructures(dict: string[]) : void {
    this.freeMoleculesStructures();
    this._rdKitMols = [];
    let logged = false;
    for (let i = 0; i < dict.length; ++i) {
      const item = dict[i];
      let mol;
      if (!item || item === '')
        mol = this._rdKitModule.get_mol('');
      else {
        mol = getMolSafe(item, {}, this._rdKitModule).mol;
        if (mol === null) {
          if (!logged) {
            const errorMessage = 'Chem | Possibly a malformed molString at init: `' + item + '`';
            logged = true;
          }
        }
      }
      this._rdKitMols.push(mol);
    }
  }

  getQMol(molString: string) : RDMol | null {
    let mol = null;
    try { mol = this._rdKitModule.get_qmol(molString); }
    catch(e) {
      if (mol !== null && mol.is_valid())
        mol.delete();
      return null;
    }
    validateMol(mol, molString);
    return mol;
  }

  searchSubstructure(queryMolString: string, queryMolBlockFailover: string, bitset?: boolean[]): string {
    const matches: number[] = [];
    if (this._rdKitMols === null)
      return '[' + matches.join(', ') + ']';

    let queryMol: RDMol | null = null;

    if (isMolBlock(queryMolString)) {
      if (queryMolString.includes(' H ') || queryMolString.includes('V3000'))
        queryMol = getMolSafe(queryMolString, {mergeQueryHs: true}, this._rdKitModule).mol;
      else {
        const molTmp = getMolSafe(queryMolString, {"mergeQueryHs":true, "kekulize": true}, this._rdKitModule).mol;
        if (molTmp !== null) {
          let molBlockAroma = null;
          try { molBlockAroma = molTmp!.get_aromatic_form(); }
          catch(e) { // looks like we get here when the molecule is already aromatic, so we just re-assign the block
            molBlockAroma = queryMolString;
          }

          molTmp.delete();
          const newQueryMolString = syncQueryAromatics(molBlockAroma, queryMolString);
          queryMolString = newQueryMolString;
        }
        queryMol = this.getQMol(queryMolString);
      }
    } else { // not a molblock
      queryMol = this.getQMol(queryMolString);
      if (queryMol !== null) {
        const mol = getMolSafe(queryMolString, {mergeQueryHs: true}, this._rdKitModule).mol;
        if (mol !== null) { // check the qmol is proper
          const match = mol.get_substruct_match(queryMol);
          if (match === '{}') {
            queryMol = mol;
          } else mol.delete();
        } // else, this looks to be a real SMARTS
      } else { // failover to queryMolBlockFailover
        queryMol = getMolSafe(queryMolBlockFailover, {mergeQueryHs: true}, this._rdKitModule).mol; // possibly get rid of fall-over in future
      }
    }

    if (queryMol !== null) {
        if (bitset) {
          for (let i = 0; i < bitset.length; ++i) {
            if (bitset[i] && this._rdKitMols[i] && this._rdKitMols[i]!.get_substruct_match(queryMol) !== '{}') // Is patternFP iff?
              matches.push(i);
          }
        } else {
          for (let i = 0; i < this._rdKitMols!.length; ++i) {
            if (this._rdKitMols[i] && this._rdKitMols[i]!.get_substruct_match(queryMol) !== '{}')
              matches.push(i);
          }
        }
        queryMol.delete();
    } else
      throw new Error('Chem | Search pattern cannot be set');

    return '[' + matches.join(', ') + ']';
  }

  freeMoleculesStructures(): void {
    if (this._rdKitMols !== null) {
      for (const mol of this._rdKitMols!)
        mol?.delete();
      this._rdKitMols = null;
    }
  }

  convertMolNotation(targetNotation: string): string[] {
    let result = (targetNotation === MolNotation.MolBlock) ? MALFORMED_MOL_V2000 :
      (targetNotation === MolNotation.V3KMolBlock) ? MALFORMED_MOL_V3000 : 'MALFORMED_INPUT_VALUE';

    const results = new Array(this._rdKitMols!.length);
    let mol: RDMol | null = null;
    for (let i = 0; i < this._rdKitMols!.length; ++i) {
      mol = this._rdKitMols![i];
      if (mol) {
        if (targetNotation === MolNotation.MolBlock) {
          if (!mol.has_coords())
            mol.set_new_coords();
         result = mol.get_molblock();
        }
        else if (targetNotation === MolNotation.Smiles)
          result = mol.get_smiles();
        else if (targetNotation === MolNotation.V3KMolBlock)
          result = mol.get_v3Kmolblock();
        else if (targetNotation === MolNotation.Smarts)
          result = mol.get_smarts();
      }
     results[i] = result;
    }

    console.log('Finished Worker ' + results.length);
    return results;
  }
}
