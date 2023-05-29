import {RdKitServiceWorkerSimilarity} from './rdkit-service-worker-similarity';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {syncQueryAromatics} from '../utils/aromatic-utils';
import {getMolSafe} from '../utils/mol-creation_rdkit';
import {isMolBlock} from '../utils/chem-common';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {RuleId} from '../panels/structural-alerts';

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
  if (!mol.is_valid()) {
    mol.delete();
    throw new Error('FATAL RDKit Error: Created a not valid molecule with no exception ' + molString);
  }
}

export class RdKitServiceWorkerSubstructure extends RdKitServiceWorkerSimilarity {
  constructor(module: RDModule, webRoot: string) {
    super(module, webRoot);
  }

  initMoleculesStructures(molecules: string[]): number {
    this.freeMoleculesStructures();
    this._rdKitMols = new Array<RDMol | null>(molecules.length).fill(null);
    let malformed = 0;
    for (let i = 0; i < molecules.length; ++i) {
      const item = molecules[i];
      if (item && item !== '') {
        const molSafe = getMolSafe(item, {}, this._rdKitModule);
        const mol = molSafe.mol;
        if (mol) {
          mol.is_qmol = molSafe.isQMol;
          this._rdKitMols[i] = mol;
        } else
          malformed++;
      }
    }
    return malformed;
  }

  getQMol(molString: string) : RDMol | null {
    let mol = null;
    try {mol = this._rdKitModule.get_qmol(molString);} catch (e) {
      if (mol !== null)
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

    let queryMol: RDMol | null;

    if (isMolBlock(queryMolString)) {
      if (queryMolString.includes(' H ') || queryMolString.includes('V3000'))
        queryMol = getMolSafe(queryMolString, {mergeQueryHs: true}, this._rdKitModule).mol;
      else {
        const molTmp = getMolSafe(queryMolString, {'mergeQueryHs': true, 'kekulize': true}, this._rdKitModule).mol;
        if (molTmp !== null) {
          let molBlockAroma: string | null;
          try {molBlockAroma = molTmp!.get_aromatic_form();} catch (e) {
            // looks like we get here when the molecule is already aromatic, so we just re-assign the block
            molBlockAroma = queryMolString;
          }

          molTmp.delete();
          queryMolString = syncQueryAromatics(molBlockAroma, queryMolString);
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
            queryMol.delete(); //remove mol object previously stored in queryMol
            queryMol = mol;
          } else
            mol.delete();
        } // else, this looks to be a real SMARTS
      } else { // failover to queryMolBlockFailover
        // possibly get rid of fall-over in future
        queryMol = getMolSafe(queryMolBlockFailover, {mergeQueryHs: true}, this._rdKitModule).mol;
      }
    }

    if (queryMol !== null) {
      if (bitset) {
        for (let i = 0; i < bitset.length; ++i) {
          try { //get_substruct_match can potentially throw an exception, need to catch
            // Is patternFP iff?
            if (bitset[i] && this._rdKitMols[i] && this._rdKitMols[i]!.get_substruct_match(queryMol) !== '{}')
              matches.push(i);
          } catch {
            continue;
          }
        }
      } else {
        for (let i = 0; i < this._rdKitMols!.length; ++i) {
          try {
            if (this._rdKitMols[i] && this._rdKitMols[i]!.get_substruct_match(queryMol) !== '{}')
              matches.push(i);
          } catch {
            continue;
          }
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
        } else if (targetNotation === MolNotation.Smiles)
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

  mostCommonStructure(molecules: string[], exactAtomSearch: boolean, exactBondSearch: boolean): string {
    const mols: RDMol[] = [];

    for (let i = 0; i < molecules.length; i++) {
      const molString = molecules[i];
      const molSafe = getMolSafe(molString!, {}, this._rdKitModule);
      if (molSafe.mol !== null && !molSafe.isQMol && molSafe.mol.is_valid())
        mols.push(molSafe.mol);
      else
        molSafe.mol?.delete();
    }
    if (mols.length > 0) {
      const arr = new Uint32Array(mols.length);

      for (let i = 0; i < mols.length; i++) {
      //@ts-ignore
        arr[i] = mols[i].$$.ptr;
      }

      //as wasm works with 8 bit and 32 bit are used
      const buff = this._rdKitModule._malloc(mols.length * 4);

      // >> 2 is the reduction of element number of 32 bit vs 8 bit for offset
      //@ts-ignore
      this._rdKitModule.HEAPU32.set(arr, buff >> 2);

      const smarts: string = this._rdKitModule.get_mcs(buff, mols.length, exactAtomSearch, exactBondSearch);

      this._rdKitModule._free(buff);

      for (let j = 0; j < mols.length; j++)
        mols[j].delete();
      return smarts;
    }
    return '';
  }

  getStructuralAlerts(alerts: {[rule in RuleId]?: string[]}, molecules?: string[]): {[rule in RuleId]?: boolean[]} {
    if (this._rdKitMols === null && typeof molecules === 'undefined') {
      console.debug(`getStructuralAlerts: No molecules to process`);
      return {};
    }

    const ruleSmartsMap: {[rule in RuleId]?: (RDMol | null)[]} = {};
    const rules = Object.keys(alerts) as RuleId[];
    for (const rule of rules) {
      const ruleLength = alerts[rule]!.length;
      ruleSmartsMap[rule] = new Array(ruleLength);
      for (let smartsIdx = 0; smartsIdx < ruleLength; smartsIdx++)
        ruleSmartsMap[rule]![smartsIdx] = this.getQMol(alerts[rule]![smartsIdx]);
    }

    const molsCount = molecules?.length ?? this._rdKitMols!.length;
    // Prepare the result storage
    const resultValues: {[ruleId in RuleId]?: BitArray} = {};
    for (const rule of rules)
      resultValues[rule] = new BitArray(molsCount, false);

    // Run the structural alerts detection
    for (let molIdx = 0; molIdx < molsCount; molIdx++) {
      const mol = typeof molecules !== 'undefined' ? getMolSafe(molecules[molIdx], {}, this._rdKitModule).mol :
        this._rdKitMols![molIdx];
      if (mol === null) {
        console.debug(`Molecule ${molIdx} is null`);
        continue;
      }

      for (const rule of rules) {
        const ruleSmarts = ruleSmartsMap[rule]!;
        for (let alertIdx = 0; alertIdx < ruleSmarts.length; alertIdx++) {
          const smarts = ruleSmarts[alertIdx];
          if (smarts === null) {
            console.debug(`Smarts ${alertIdx} for rule ${rule} is null at molecule ${molIdx}`);
            continue;
          }

          const matches = mol.get_substruct_match(smarts);
          if (matches !== '{}') {
            resultValues[rule]!.setTrue(molIdx);
            break;
          }
        }
      }
      mol.delete();
    }

    for (const smartsList of Object.values(ruleSmartsMap)) {
      for (const smarts of smartsList)
        smarts?.delete();
    }

    this._rdKitMols = null;

    return Object.fromEntries(Object.entries(resultValues).map(([k, val]) => [k, val.getRangeAsList(0, val.length)]));
  }
}
