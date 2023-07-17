import {RdKitServiceWorkerSimilarity} from './rdkit-service-worker-similarity';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {IMolContext, getMolSafe, isSmarts} from '../utils/mol-creation_rdkit';
import {isMolBlock} from '../utils/chem-common';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import { RuleId } from '../panels/structural-alerts';

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

  searchSubstructure(queryMolString: string, queryMolBlockFailover: string, bitset?: boolean[]): string {
    const matches: number[] = [];
    if (this._rdKitMols === null)
      return '[' + matches.join(', ') + ']';

    let queryMol = null;

    if (isMolBlock(queryMolString)) {
      if (queryMolString.includes(' H ') || queryMolString.includes('V3000'))
        queryMol = getMolSafe(queryMolString, {mergeQueryHs: true}, this._rdKitModule).mol;
      else {
        try {
          queryMol = this._rdKitModule.get_qmol(queryMolString);
          queryMol.set_aromatic_form();
        } catch (e) {
          if (queryMol) {
            queryMol.delete();
            queryMol = null;
          }
        }
      }
    } else { // not a molblock
      try {
        queryMol = this._rdKitModule.get_qmol(queryMolString);
      } catch (e) {}
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

  getStructuralAlerts(alerts: {[rule in RuleId]?: string[]}, molecules?: string[]): {[rule in RuleId]?: boolean[]} {
    if (this._rdKitMols === null && typeof molecules === 'undefined') {
      console.debug(`getStructuralAlerts: No molecules to process`);
      return {};
    }

    const ruleSmartsMap: {[rule in RuleId]?: (RDMol | null)[]} = {};
    const rules = Object.keys(alerts) as RuleId[];
    for (let rule of rules) {
      const ruleLength = alerts[rule]!.length;
      ruleSmartsMap[rule] = new Array(ruleLength);
      for (let smartsIdx = 0; smartsIdx < ruleLength; smartsIdx++) {
        let qmol = null;
        try {
          qmol = this._rdKitModule.get_qmol(alerts[rule]![smartsIdx]);
        } catch {}
        ruleSmartsMap[rule]![smartsIdx] = qmol;
      }
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

  getMolSafe(molecules: string[], details: object = {}, warnOff: boolean = true, checkIfSmarts: boolean = true): IMolContext[] {
    const mols: IMolContext[] = [];
    for (let i = 0; i < molecules.length; ++i)
      mols[i] = getMolSafe(molecules[i], details, this._rdKitModule, warnOff, checkIfSmarts);
    return mols;
  }
}
