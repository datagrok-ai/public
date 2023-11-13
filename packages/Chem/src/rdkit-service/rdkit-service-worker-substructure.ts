import {RdKitServiceWorkerSimilarity} from './rdkit-service-worker-similarity';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getMolSafe, getQueryMolSafe} from '../utils/mol-creation_rdkit';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {RuleId} from '../panels/structural-alerts';
import { SubstructureSearchType } from '../constants';

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
    console.log(`initMoleculesStructures`);
    this.freeMoleculesStructures();
    this._rdKitMols = new Array<RDMol | null>(molecules.length).fill(null);
    let numMalformed = 0;
    for (let i = 0; i < molecules.length; ++i) {
      const item = molecules[i];
      if (item && item !== '') {
        const molSafe = getMolSafe(item, {}, this._rdKitModule);
        const mol = molSafe.mol;
        if (mol)
            this._rdKitMols![i] = mol;
        else
          numMalformed++;
      }
    }
    return numMalformed;
  }

  async searchSubstructure(queryMolString: string, queryMolBlockFailover: string, molecules?: string[],
    searchType?: SubstructureSearchType): Promise<Uint32Array> {
    if (!molecules)
      throw new Error('Chem | Molecules for substructure serach haven\'t been provided');

    const queryMol = getQueryMolSafe(queryMolString, queryMolBlockFailover, this._rdKitModule); 
    let queryCanonicalSmiles = '';
    if (queryMol !== null) {
      if (searchType === SubstructureSearchType.EXACT_MATCH) {
        try {
          queryCanonicalSmiles = queryMol.get_smiles();
          //need to get canonical smiles from mol (not qmol) since qmol implicitly merges query hydrogens
          queryCanonicalSmiles = this._rdKitModule.get_mol(queryMolString).get_smiles();
        } catch {}
      }
      const matches = await this.searchWithPatternFps(queryMol, molecules,
        searchType ?? SubstructureSearchType.CONTAINS, queryCanonicalSmiles);
      queryMol.delete();
      return matches;
    } else
      throw new Error('Chem | Search pattern cannot be set');
  }


  async searchWithPatternFps(queryMol: RDMol, molecules: string[], searchType: SubstructureSearchType,
    queryCanonicalSmiles: string): Promise<Uint32Array> {
    const matches = new BitArray(molecules.length);
    if (this._requestTerminated)
      return matches.buffer;
    const details = JSON.stringify({sanitize: false, removeHs: false, assignStereo: false});
    for (let i = 0; i < molecules.length; ++i) {
      const terminationCheckDelay = queryCanonicalSmiles ? this._terminationCheckDelay * 10 : this._terminationCheckDelay;

      if (i % terminationCheckDelay === 0) //every N molecules check for termination flag
        await new Promise((r) => setTimeout(r, 0));
      if (this._requestTerminated)
        return matches.buffer;

      if (queryCanonicalSmiles) {
        matches.setFast(i, molecules[i] === queryCanonicalSmiles);
      } else {
        let mol: RDMol | null = null;
        let isCached = false;
        try {
          const cachedMol = this._molsCache?.get(molecules[i]);
          mol = cachedMol ?? this._rdKitModule.get_mol(molecules[i], details);
          if (cachedMol || this.addToCache(mol))
            isCached = true;
          if (mol) {
            if (this.searchBySearchType(mol, queryMol, searchType))
              matches.setFast(i, true);
          }
        } catch {
          continue;
        } finally {
          !isCached && mol?.delete();
        }
      }
    }
    return matches.buffer;
  }

  searchBySearchType(mol: RDMol, queryMol: RDMol, searchType: SubstructureSearchType): boolean {
    switch (searchType) {
      case SubstructureSearchType.CONTAINS:
        return mol.get_substruct_match(queryMol) !== '{}';
      case SubstructureSearchType.INCLUDED_IN:
        return queryMol.get_substruct_match(mol) !== '{}';
      case SubstructureSearchType.NOT_CONTAINS:
          return mol.get_substruct_match(queryMol) == '{}';
      case SubstructureSearchType.NOT_INCLUDED_IN:
          return queryMol.get_substruct_match(mol) == '{}';
      case SubstructureSearchType.EXACT_MATCH:
        const match1 = mol.get_substruct_match(queryMol);
        const match2 = queryMol.get_substruct_match(mol);
        return match1 !== '{}' && match2 !== '{}';
      default:
        throw Error('Unknown search type: ' + searchType);
    }
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

  invalidateCache() {
    this._cacheCounter = 0;
    if (this._molsCache) {
      this._molsCache.forEach((it) => it?.delete());
      this._molsCache.clear();
    }
  }

  mmpGetFragments(molecules: string[]): [string, string][][] {
    const frags: [string, string][][] = new Array<[string, string][]>(molecules.length);
    for (let i = 0; i < molecules.length; i++) {
      let mol;
      try {
        mol = this._rdKitModule.get_mol(molecules[i]);
        if (mol) {
          const res = this._rdKitModule.get_matched_fragments(mol, 1, 1, 20);
          const length = res.second.size();
          frags[i] = new Array<[string, string]>(length);

          for (let j = 0; j < length; j++) {
            try {
              const fffSplit = res.second.next().get_smiles().split('.');
              const firstIsFirst = fffSplit[0].length >= fffSplit[1].length;
              frags[i][j] = [firstIsFirst ? fffSplit[0] : fffSplit[1], firstIsFirst ? fffSplit[1] : fffSplit[0]];
            } catch (e: any) {
              frags[i][j] = ['', ''];
            } finally {

            }
          }

          res.first.delete();
          res.second.delete();
        } else
          frags[i] = new Array<[string, string]>(0);
      } catch (e: any) {
        frags[i] = new Array<[string, string]>(0);
      } finally {
        mol?.delete();
      }
    }

    return frags;
  }

  mmpGetMcs(molecules: [string, string][]): string[] {
    const res: string[] = new Array<string>(molecules.length);
    for (let i = 0; i < molecules.length; i++) {
      let mol1;
      let mol2;
      let mols;
      try {
        mol1 = getMolSafe(molecules[i][0], {}, this._rdKitModule);
        mol2 = getMolSafe(molecules[i][1], {}, this._rdKitModule);
        mols = new this._rdKitModule.MolList();
        if (mol1.mol && mol2.mol) {
          mols.append(mol1.mol!);
          mols.append(mol2.mol!);
          res[i] = this._rdKitModule.get_mcs_as_smarts(mols, JSON.stringify({
            AtomCompare: 'Elements',
            BondCompare: 'OrderExact',
            RingMatchesRingOnly: true
          }));
        } else
          res[i] = '';
      } catch (e: any) {
        res[i] = '';
      } finally {
        mol1?.mol?.delete();
        mol2?.mol?.delete();
      }
    }

    return res;
  }

  setTerminateFlag(flag: boolean) {
    this._requestTerminated = flag;
  }

  mostCommonStructure(molecules: string[], exactAtomSearch: boolean, exactBondSearch: boolean): string {
    let mols;
    try {
      mols = new this._rdKitModule.MolList();
      for (let i = 0; i < molecules.length; i++) {
        const molString = molecules[i];
        if (!molString)
          continue;
        let molSafe;
        try {
          molSafe = getMolSafe(molString!, {}, this._rdKitModule);
          if (molSafe.mol !== null && !molSafe.isQMol)
            mols.append(molSafe.mol);
        } finally {
          molSafe?.mol?.delete();
        }
      }
      let mcsSmarts: string|null = null;
      if (mols.size() > 1) {
        mcsSmarts = this._rdKitModule.get_mcs_as_smarts(mols, JSON.stringify({
          AtomCompare: exactAtomSearch ? 'Elements' : 'Any',
          BondCompare: exactBondSearch ? 'OrderExact' : 'Order',
          //RingMatchesRingOnly: true
        }));
      }
      return mcsSmarts ?? '';
    } finally {
      mols?.delete();
    }
  }
}
