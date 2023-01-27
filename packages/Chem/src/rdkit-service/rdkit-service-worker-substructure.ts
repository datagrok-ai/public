import {RdKitServiceWorkerSimilarity} from './rdkit-service-worker-similarity';
import {isMolBlock} from '../utils/convert-notation-utils';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
//import {aromatizeMolBlock} from "../utils/aromatic-utils";

function syncQueryAromatics_1(molBlock: string,  bonds2Change: Array<number> | null = null) : string | Array<number> {
  let curPos = 0;
  curPos = molBlock.indexOf('\n', curPos) + 1;
  curPos = molBlock.indexOf('\n', curPos) + 1;
  curPos = molBlock.indexOf('\n', curPos) + 1;
  const atomCounts = parseInt(molBlock.substring(curPos, curPos + 3));
  const bondCounts = parseInt(molBlock.substring(curPos + 3, curPos + 6));

  for (let atomRowI = 0; atomRowI < atomCounts; atomRowI++) {
    curPos = molBlock.indexOf('\n', curPos) + 1;
  }

  const read = bonds2Change === null;
  bonds2Change ??= [];

  let bondOrder = -1;
  for (let bondRowI = 0; bondRowI < bondCounts; bondRowI++) {
    curPos = molBlock.indexOf('\n', curPos) + 1;
    if (read) {
      bondOrder = parseInt(molBlock.substring(curPos + 8, curPos + 9));
      if (bondOrder === 4)
       bonds2Change.push(bondRowI);
    }
    else {
      if (bonds2Change.includes(bondRowI))
        molBlock = molBlock.slice(0, curPos + 8) + '4' + molBlock.slice(curPos + 9);
    }
  }

  return read ? bonds2Change : molBlock;
}

function syncQueryAromatics_2(molBlockAroma: string, molBlock : string) : string {
  const bonds2Change = syncQueryAromatics_1(molBlock);
  const molModified = syncQueryAromatics_1(molBlockAroma, bonds2Change as Array<number>);
  return molModified as string;
}

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

    for (let i = 0; i < dict.length; ++i) {
      const item = dict[i];
      let mol = this.getMol(item, null, true);
      if (mol === null) {
        console.error('Chem | Possibly a malformed molString at init: `' + item + '`');
        mol = this._rdKitModule.get_mol('');
      }
      this._rdKitMols.push(mol);
    }
  }

  getMol(molString: string, options: string | null = null, qmolFallback: boolean = false) : RDMol | null {
    options ??= '{"mergeQueryHs":true}';
    let mol = null;
    try { mol = this._rdKitModule.get_mol(molString, options); }
    catch (e) {
      if (mol !== null && mol.is_valid())
        mol.delete();

      if (qmolFallback) {
        mol = this.getQMol(molString);
        if (mol === null) {
          console.error('Chem | Failed fallback to qmol. Possibly a malformed query: `' + molString + '`');
          return null;
        }
        validateMol(mol, molString);
        //console.error('Successfully performed fallback on qmol ' + molString);
        return mol;
      } else return null;
    }
    validateMol(mol, molString);
    return mol;
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
        queryMol = this.getMol(queryMolString);
      else {
        const molTmp = this.getMol(queryMolString, '{"mergeQueryHs":true, "kekulize": true}', true);
        if (molTmp !== null) {
          let molBlockAroma = null;
          try { molBlockAroma = molTmp!.get_aromatic_form(); }
          catch(e) { // looks like we get here when the molecule is already aromatic, so we just re-assign the block
            molBlockAroma = queryMolString;
          }

          molTmp.delete();
          const newQueryMolString = syncQueryAromatics_2(molBlockAroma, queryMolString);
          queryMolString = newQueryMolString;
          //const newQueryMolString = aromatizeMolBlock(queryMolString);
        }
        queryMol = this.getQMol(queryMolString);
      }
    } else { // not a molblock
      queryMol = this.getQMol(queryMolString);
      if (queryMol !== null) {
        const mol = this.getMol(queryMolString);
        if (mol !== null) { // check the qmol is proper
          const match = mol.get_substruct_match(queryMol);
          if (match === '{}') {
            queryMol = mol;
          } else mol.delete();
        } // else, this looks to be a real SMARTS
      } else { // failover to queryMolBlockFailover
        queryMol = this.getMol(queryMolBlockFailover); // possibly get rid of fall-over in future
      }
    }

    if (queryMol !== null) {
        if (bitset) {
          for (let i = 0; i < bitset.length; ++i) {
            if (bitset[i] && this._rdKitMols[i]!.get_substruct_match(queryMol) !== '{}') // Is patternFP iff?
              matches.push(i);
          }
        } else {
          for (let i = 0; i < this._rdKitMols!.length; ++i) {
            if (this._rdKitMols[i]!.get_substruct_match(queryMol) !== '{}')
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
        mol.delete();
      this._rdKitMols = null;
    }
  }
}
