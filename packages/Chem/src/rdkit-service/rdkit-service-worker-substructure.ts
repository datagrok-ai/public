import {RdKitServiceWorkerSimilarity} from './rdkit-service-worker-similarity';
import {MolList, RDModule, RDMol, RGroupDecomp} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {IMolContext, getMolSafe, getQueryMolSafe} from '../utils/mol-creation_rdkit';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {RuleId} from '../panels/structural-alerts';
import {SubstructureSearchType} from '../constants';
import {stringArrayToMolList, getUncommonAtomsAndBonds} from '../utils/chem-common';
import {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';

export enum MolNotation {
  Smiles = 'smiles',
  Smarts = 'smarts',
  MolBlock = 'molblock', // molblock V2000
  V3KMolBlock = 'v3Kmolblock', // molblock V3000
  Unknown = 'unknown',
}

export interface IRGroupAnalysisResult {
  colNames: string[];
  smiles: Array<string>[];
  atomsToHighLight: Array<Array<Uint32Array>>;
  bondsToHighLight: Array<Array<Uint32Array>>;
}

export interface IMmpFragmentsResult {
  frags: [string, string][][];
  smiles: string[];
}

export type InverseSubstructureRes = {
  inverse1: (ISubstruct | null)[],
  inverse2: (ISubstruct | null)[],
  fromAligned: string[],
  toAligned: string[]
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
      const terminationCheckDelay = queryCanonicalSmiles ?
        this._terminationCheckDelay * 10 :
        this._terminationCheckDelay;

      if (i % terminationCheckDelay === 0) //every N molecules check for termination flag
        await new Promise((r) => setTimeout(r, 0));
      if (this._requestTerminated)
        return matches.buffer;

      if (queryCanonicalSmiles)
        matches.setFast(i, molecules[i] === queryCanonicalSmiles);
      else {
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

  async convertMolNotation(molecules: string[], targetNotation: string): Promise<string[]> {
    if (!molecules || this._requestTerminated)
      return [];
    let addedToCache = false;
    const result = (targetNotation === MolNotation.MolBlock) ? MALFORMED_MOL_V2000 :
      (targetNotation === MolNotation.V3KMolBlock) ? MALFORMED_MOL_V3000 : 'MALFORMED_INPUT_VALUE';
    const results = new Array<string>(molecules.length).fill(result);

    for (let i = 0; i < molecules!.length; ++i) {
      //every N molecules check for termination flag
      if (i % this._terminationCheckDelay === 0)
        await new Promise((r) => setTimeout(r, 0));
      if (this._requestTerminated)
        return results;
      const item = molecules[i];
      if (!item || item === '') {
        results[i] = '';
        continue;
      }
      let rdMol = this._molsCache?.get(molecules[i]);
      if (!rdMol) {
        const mol: IMolContext = getMolSafe(item, {}, this._rdKitModule);
        rdMol = mol?.mol;
        if (rdMol)
          rdMol.is_qmol = mol?.isQMol;
      }
      if (rdMol) {
        try {
          switch (targetNotation) {
          case MolNotation.MolBlock:
            if (!rdMol.has_coords())
              rdMol.set_new_coords();
            results[i] = rdMol.get_molblock();
            break;
          case MolNotation.Smiles:
            if (!rdMol.is_qmol)
              results[i] = rdMol.get_smiles();
            break;
          case MolNotation.V3KMolBlock:
            results[i] = rdMol.get_v3Kmolblock();
            break;
          case MolNotation.Smarts:
            results[i] = rdMol.get_smarts();
            break;
          default:
            results[i] = 'Unknown notation: ' + targetNotation;
            throw Error('Unknown notation: ' + targetNotation);
          }
          addedToCache = this.addToCache(rdMol);
        } catch {
          // nothing to do, fp is already null
        } finally {
          if (!addedToCache) {
            //do not delete mol in case it is in cache
            rdMol?.delete();
          }
        }
      }
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
    if (molsCount === 0)
      return {};
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

  rGroupAnalysis(molecules: string[], coreMolecule: string, coreIsQMol?: boolean, options?: string): IRGroupAnalysisResult {
    let mols: MolList | null = null;
    let res: RGroupDecomp | null = null;
    let core: RDMol | null = null;
    let cols: {[colName: string] : MolList | null} = {};
    let totalColsNum = 0;
    let colNames: string [] = [];
    const rgroupTargetAtomsPropName = '_rgroupTargetAtoms';
    const rgroupTargetBondsPropName = '_rgroupTargetBonds';
    const numOfNonRGroupCols = 1; //Core column
    const coreColName = 'Core';
    const molColName = 'Mol';
    const resCols = [];
    try {
      mols = stringArrayToMolList(molecules, this._rdKitModule);
      try {
        core = coreIsQMol ? this._rdKitModule.get_qmol(coreMolecule) :
          this._rdKitModule.get_mol(coreMolecule, JSON.stringify({
            makeDummiesQueries: true,
            mappedDummiesAreRGroups: true,
          }));
        if (!core)
          throw new Error(`Core is possibly malformed`);
      } catch (e) {
        throw new Error(`Core is possibly malformed`);
      }

      //res = this._rdKitModule.rgroups(core!, mols!, options ? options : '');
      res = this._rdKitModule.get_rgd(core!, options ? options : '');
      const unmatches: number[] = [];
      for (let i = 0; i < molecules.length; i ++) {
        const match = res!.add(mols!.at(i));
        if (match == -1)
          unmatches.push(i);
      }

      res!.process();

      cols = res!.get_rgroups_as_columns();
      colNames = Object.keys(cols).filter((it) => it !== molColName); //exclude Mol column from result since we do not need it

      totalColsNum = colNames.length;

      let counter = 0;
      if (totalColsNum > 0) {
        const atomsToHighlight = Array<Array<Uint32Array>>(totalColsNum - numOfNonRGroupCols)
          .fill([]).map((u) => [] as Uint32Array[]);
        const bondsToHighlight = Array<Array<Uint32Array>>(totalColsNum - numOfNonRGroupCols)
          .fill([]).map((u) => [] as Uint32Array[]);
        for (let i = 0; i < totalColsNum; i++) {
          const isRGroupCol = colNames[i] !== coreColName;
          const col = Array<string>(molecules.length);
          for (let j = 0; j < molecules.length; j++) {
            if (unmatches[counter] !== j) {
              const rgroup = cols[colNames[i]]!.at(j - counter);
              if (isRGroupCol) {
                if (rgroup.has_prop(rgroupTargetAtomsPropName))
                  atomsToHighlight[i - numOfNonRGroupCols][j] = new Uint32Array(JSON.parse(rgroup.get_prop(rgroupTargetAtomsPropName)));
                if (rgroup.has_prop(rgroupTargetBondsPropName))
                  bondsToHighlight[i - numOfNonRGroupCols][j] = new Uint32Array(JSON.parse(rgroup.get_prop(rgroupTargetBondsPropName)));
              }
              col[j] = rgroup.get_smiles();
            } else {
              counter++;
              col[j] = '';
            }
          }
          resCols.push(col);
          counter = 0;
        }
        return {colNames: colNames, smiles: resCols, atomsToHighLight: atomsToHighlight, bondsToHighLight: bondsToHighlight};
      }
      return {colNames: [], smiles: [], atomsToHighLight: [], bondsToHighLight: []};
    } catch (e: any) {
      throw new Error(e.message);
    } finally {
      core?.delete();
      mols?.delete();
      res?.delete();
      for (let i = 0; i < totalColsNum; i ++)
        cols[colNames[i]]?.delete();
    }
  }

  invalidateCache() {
    this._cacheCounter = 0;
    if (this._molsCache) {
      this._molsCache.forEach((it) => it?.delete());
      this._molsCache.clear();
    }
  }

  mmpGetFragments(molecules: string[]): IMmpFragmentsResult {
    const size = molecules.length;
    const frags: [string, string][][] = new Array<[string, string][]>(size);
    const smiles = new Array<string>(size);
    for (let i = 0; i < size; i++) {
      let mol;
      try {
        mol = this._rdKitModule.get_mol(molecules[i]);
        smiles[i] = mol.get_smiles();
        if (mol) {
          const res = mol.get_mmpa_frags(1, 1, 20);
          const length = res.sidechains.size();
          frags[i] = new Array<[string, string]>(length);

          for (let j = 0; j < length; j++) {
            let frag = null;
            try {
              frag = res.sidechains.next();
              const split = frag.get_smiles().split('.');

              //the following logic is for case when additional entities, like salts, present
              let firstFragment = '';
              let secondFragment = '';
              let additionalFragments: string[];

              if (split.length == 2) {
                firstFragment = split[0];
                secondFragment = split[1];
                additionalFragments = [];
              } else {
                let oneFragmentReady = false;
                additionalFragments = new Array<string>(split.length - 2);
                let counter = 0;

                for (let k = 0; k < split.length; k++) {
                  if (split[k].includes('[*')) {
                    if (oneFragmentReady)
                      secondFragment = split[k];
                    else {
                      firstFragment = split[k];
                      oneFragmentReady = true;
                    }
                  } else {
                    additionalFragments[counter] = split[k];
                    counter++;
                  }
                }
              }

              const firstIsFirst = firstFragment.length >= secondFragment.length;

              //swap
              if (!firstIsFirst) {
                const temp = firstFragment;
                firstFragment = secondFragment;
                secondFragment = temp;
              }

              //add additional entities to smallest fragment
              for (let k = 0; k < additionalFragments.length; k++)
                secondFragment += '.' + additionalFragments[k];

              frags[i][j] = [firstFragment, secondFragment];
            } catch (e: any) {
              frags[i][j] = ['', ''];
            } finally {
              frag?.delete();
            }
          }

          res.cores.delete();
          res.sidechains.delete();
        } else
          frags[i] = new Array<[string, string]>(0);
      } catch (e: any) {
        frags[i] = new Array<[string, string]>(0);
        smiles[i] = '';
      } finally {
        mol?.delete();
      }
    }

    return {frags, smiles};
  }

  mmpLinkFragments(cores: string[], fragments: string[]): string[] {
    const size = cores.length;
    const smiles = new Array<string>(size);
    for (let i = 0; i < size; i++) {
      let mol;
      let smilesGen = '';
      try {
        const smi = `${cores[i]}.${fragments[i]}`.replaceAll('([*:1])', '9').replaceAll('[*:1]', '9');
        mol = getMolSafe(smi, {}, this._rdKitModule);
        smilesGen = mol.mol!.get_smiles();
        smiles[i] = smilesGen;
      } catch (e: any) {
        smiles[i] = '';
      } finally {
        mol?.mol?.delete();
      }
    }

    return smiles;
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
            RingMatchesRingOnly: true,
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


  getInverseSubstructuresAndAlign(cores: string[],
    from: string[], to: string[]):
    InverseSubstructureRes {
    const fromAligned = new Array<string>(from.length);
    const toAligned = new Array<string>(from.length);
    const res1 = new Array<(ISubstruct | null)>(from.length);
    const res2 = new Array<(ISubstruct | null)>(from.length);

    for (let i = 0; i < from.length; i++) {
      //aligning molecules
      let mol1 = null;
      let mol2 = null;
      let mcsMol = null;

      const opts = JSON.stringify({
        useCoordGen: true,
        allowRGroups: true,
        acceptFailure: false,
        alignOnly: true,
      });

      try {
        const core = cores[i].replace('[*:1]', '[H]');
        mcsMol = this._rdKitModule.get_mol(core);
        mol1 = this._rdKitModule.get_mol(from[i]);
        mol2 = this._rdKitModule.get_mol(to[i]);
        mcsMol.set_new_coords();
        mol1.generate_aligned_coords(mcsMol, opts);
        mol2.generate_aligned_coords(mcsMol, opts);
        fromAligned[i] = mol1.get_molblock();
        toAligned[i] = mol2.get_molblock();
        res1[i] = getUncommonAtomsAndBonds(from[i], mcsMol, this._rdKitModule, '#bc131f');
        //@ts-ignore
        res1[i]['highlightBondWidthMultiplier'] = 40;
        res2[i] = getUncommonAtomsAndBonds(to[i], mcsMol, this._rdKitModule, '#49bead');
        //@ts-ignore
        res2[i]['highlightBondWidthMultiplier'] = 40;
      } catch (e: any) {
        fromAligned[i] = '';
        toAligned[i] = '';
      } finally {
        mol1?.delete();
        mol2?.delete();
        mcsMol?.delete();
      }
    }
    return {inverse1: res1, inverse2: res2, fromAligned: fromAligned, toAligned: toAligned};
  }
}
