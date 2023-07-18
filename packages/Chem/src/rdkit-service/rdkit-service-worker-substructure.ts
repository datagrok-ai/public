import {RdKitServiceWorkerSimilarity} from './rdkit-service-worker-similarity';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getMolSafe, getQueryMolSafe} from '../utils/mol-creation_rdkit';
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

export class RdKitServiceWorkerSubstructure extends RdKitServiceWorkerSimilarity {
  constructor(module: RDModule, webRoot: string) {
    super(module, webRoot);
  }

  initMoleculesStructures(molecules: string[], useSubstructLib?: boolean): number {
    console.log(`initMoleculesStructures`);
    this.freeMoleculesStructures();
    if (useSubstructLib) {
      this._substructLibrary = new this._rdKitModule.SubstructLibrary();
      this._malformedIdxs = new BitArray(molecules.length);
    } else
      this._rdKitMols = new Array<RDMol | null>(molecules.length).fill(null);
    let numMalformed = 0;
    for (let i = 0; i < molecules.length; ++i) {
      const item = molecules[i];
      if (item && item !== '') {
        const molSafe = getMolSafe(item, {}, this._rdKitModule);
        const mol = molSafe.mol;
        if (mol) {
          mol.is_qmol = molSafe.isQMol;
          if (useSubstructLib) {
            if (mol.is_qmol) {
              this._malformedIdxs!.setFast(i, true);
              numMalformed++;
            } else
              this._substructLibrary?.add_trusted_smiles(mol.get_smiles());
            mol.delete();
          } else
            this._rdKitMols![i] = mol;
        } else {
          if (useSubstructLib)
            this._malformedIdxs!.setFast(i, true);
          numMalformed++;
        }
      }
    }
    return numMalformed;
  }

  postTerminationFlag(flag: boolean) {
    this._terminateFlag = flag;
  }

  async searchSubstructure(queryMolString: string, queryMolBlockFailover: string, molecules?: string[],
    useSubstructLib?: boolean): Promise<Uint32Array> {
    if (this._rdKitMols === null && !useSubstructLib && !molecules)
      throw new Error('Chem | Molecules for substructure serach haven\'t been provided');

    const queryMol = getQueryMolSafe(queryMolString, queryMolBlockFailover, this._rdKitModule);

    if (queryMol !== null) {
      const matches = useSubstructLib ? this.searchWithSubstructLib(queryMol) :
        molecules ? await this.searchWithPatternFps(queryMol, molecules) : this.searchWithoutPatternFps(queryMol);
      queryMol.delete();
      return matches;
    } else
      throw new Error('Chem | Search pattern cannot be set');
  }

  searchWithSubstructLib(queryMol: RDMol): Uint32Array {
    const matches = new BitArray(this._substructLibrary!.size() + this._malformedIdxs!.trueCount());
    const str = this._substructLibrary!.get_matches(queryMol, undefined, undefined, this._substructLibrary!.size());
    const matchesIdxs = JSON.parse(str);
    // re-calculate matches idxs considering malformed data
    let nonMalformedCounter = 0;
    let matchesCounter = 0;

    for (let i = -1; (i = this._malformedIdxs!.findNext(i, false)) !== -1;) {
      if (nonMalformedCounter === matchesIdxs[matchesCounter]) {
        matchesCounter++;
        matches.setFast(i, true);
      }
      nonMalformedCounter++;
    }
    return matches.buffer;
  }

  async searchWithPatternFps(queryMol: RDMol, molecules: string[], startIdx: number = 0): Promise<Uint32Array> {
    console.log('**************in searchWithPatternFps');
    const matches = new BitArray(molecules.length);
    const calc = (start: number): Promise<void> => {
      return new Promise((resolve) => {
        setTimeout(async () => {
          this._terminateFlag && console.log(this._terminateFlag);
          if (start < molecules.length && !this._terminateFlag) {
            this.searchWithPatternFpsBatch(queryMol,
              molecules.slice(start, Math.min(start + 10, molecules.length)), matches, start);
            await calc(Math.min(start + 10, molecules.length));
          }
          if (this._terminateFlag) {
            this._terminateFlag = false;
            console.log('Terminated flag set');
          }
          resolve();
        }, 0);
      });
    };
    await calc(0);
    return matches.buffer;
  }

  searchWithPatternFpsBatch(queryMol: RDMol, molecules: string[], matches: BitArray, startIdx: number) {
    let mol: RDMol | null = null;
    const details = JSON.stringify({sanitize: false, removeHs: false, assignStereo: false});
    for (let i = 0; i < molecules.length; ++i) {
      try {
        mol = this._rdKitModule.get_mol(molecules[i], details);
        if (mol) {
          if (mol.get_substruct_match(queryMol) !== '{}')
            matches.setFast(i + startIdx, true);
        }
      } catch {
        continue;
      } finally {
        mol?.delete();
      }
    }
  }

  searchWithoutPatternFps(queryMol: RDMol): Uint32Array {
    const matches = new BitArray(this._rdKitMols!.length);
    for (let i = 0; i < this._rdKitMols!.length; ++i) {
      try {
        if (this._rdKitMols![i] && this._rdKitMols![i]!.get_substruct_match(queryMol) !== '{}')
          matches.setFast(i, true);
      } catch {
        continue;
      }
    }
    return matches.buffer;
  }


  freeMoleculesStructures(): void {
    if (this._substructLibrary) {
      this._substructLibrary.delete();
      this._substructLibrary = null;
      this._malformedIdxs = null;
    }
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
    for (const rule of rules) {
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
}
