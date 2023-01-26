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

export class RdKitServiceWorkerSubstructure extends RdKitServiceWorkerSimilarity {
  readonly _patternFpLength = 2048;
  readonly _patternFpUint8Length = 256;

  constructor(module: RDModule, webRoot: string) {
    super(module, webRoot);
  }

  initMoleculesStructures(dict: string[]) {
    this.freeMoleculesStructures();
    this._rdKitMols = [];

    for (let i = 0; i < dict.length; ++i) {
      const item = dict[i];
      let mol = null;

      try {
        mol = this._rdKitModule.get_mol(item);
      } catch (e) {
        console.error('Chem | Possibly a malformed molString: `' + item + '`');
        mol?.delete();
        mol = this._rdKitModule.get_mol('');
      }
      this._rdKitMols.push(mol);
    }
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
        if (isMolBlock(queryMolString)) {
          if (queryMolString.includes(' H '))
            queryMol = this.getMol(queryMolString);
          else {
            const molTmp = this._rdKitModule.get_mol(queryMolString, '{"mergeQueryHs":true, "kekulize": true}');
            const molBlockAroma = molTmp!.get_aromatic_form();
            molTmp!.delete();
            const newQueryMolString = syncQueryAromatics_2(molBlockAroma, queryMolString);
            //const newQueryMolString = aromatizeMolBlock(queryMolString);
            queryMol = this.getQMol(newQueryMolString)
          }
        }
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
            for (let i = 0; i < this._rdKitMols!.length; ++i) {
              if (this._rdKitMols[i]!.get_substruct_match(queryMol) !== '{}')
                matches.push(i);
            }
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
      for (const mol of this._rdKitMols!)
        mol.delete();
      this._rdKitMols = null;
    }
  }
}
