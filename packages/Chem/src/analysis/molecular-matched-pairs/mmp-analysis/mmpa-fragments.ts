import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {mmpFragmentsByEngine} from '../../../engine/analysis-jobs';
import {webGPUMMP} from '@datagrok-libraries/math/src/webGPU/mmp/webGPU-mmp';
import {IMmpFragmentsResult} from '../../../rdkit-service/rdkit-service-worker-substructure';
import {MMP_CONSTRICTIONS, MMP_ERRORS, MmpFragments, MmpRules, MolecularPair} from './mmpa-misc';

/**
* Runs parallel fragmentation for molecules (on the Crux workers in Crux mode, "Chemical engine")
* @param {DG.Column} molecules column with molecules
*/
export async function getMmpFrags(molecules: string[]): Promise<[MmpFragments, string[]]> {
  const res = await mmpFragmentsByEngine(molecules);
  return [encodeFragments(res), res.smiles];
}

export async function getMmpRules(
  fragsOut: MmpFragments, fragmentCutoff: number, gpu: boolean, strictCPU = false):
  Promise<[MmpRules, number, boolean]> {
  let rules: MmpRules | null = null;
  let allCaseesNumber = 0;
  let useGpu = false;

  try {
    if (fragsOut.fragCodes.length < 10 || !gpu || strictCPU) {
      if (fragsOut.fragCodes.length > MMP_CONSTRICTIONS.CPU)
        throw new Error(MMP_ERRORS.FRAGMENTS_CPU);

      [rules, allCaseesNumber, useGpu] = getMmpRulesCPU(fragsOut, fragmentCutoff);
    } else {
      useGpu = true;
      [rules, allCaseesNumber, useGpu] = await getMmpRulesGPU(fragsOut, fragmentCutoff);
    }
  } catch (e: any) {
    const eMsg: string = e instanceof Error ? e.message : e.toString();
    if (eMsg === MMP_ERRORS.FRAGMENTS_CPU) {
      grok.shell.warning(MMP_ERRORS.GPU_ABORTED);
      grok.shell.error(MMP_ERRORS.FRAGMENTS_CPU);
      throw new Error(MMP_ERRORS.FRAGMENTS_CPU);
    }
    if (useGpu) {
      grok.shell.warning(MMP_ERRORS.GPU_ABORTED);
      [rules, allCaseesNumber, useGpu] = await getMmpRules(fragsOut, fragmentCutoff, gpu, true);
    } else {
      grok.shell.error(MMP_ERRORS.PAIRS);
      throw new Error(MMP_ERRORS.PAIRS);
    }
  }

  return [rules!, allCaseesNumber, useGpu];
}

function getBestFragmentPair(
  dimFirst: number, idxFirst: number,
  dimSecond: number, idxSecond: number,
  fragsOut: MmpFragments): [number, number, number] | null {
  let core: number | null = null;
  let r1: number | null = null; // molecule minus core for first molecule in pair
  let r2: number | null = null; // molecule minus core for second molecule in pair

  //here we get the best possible fragment pair
  for (let p1 = 0; p1 < dimFirst; p1++) {
    for (let p2 = 0; p2 < dimSecond; p2++) {
      if (fragsOut.fragCodes[idxFirst][p1][0] === fragsOut.fragCodes[idxSecond][p2][0]) {
        const newCore = fragsOut.fragCodes[idxFirst][p1][0];
        // core is a fragment index (0 is valid) — compare against null, not falsiness
        if (core === null || fragsOut.sizes[newCore] > fragsOut.sizes[core]) {
          core = newCore;
          r1 = fragsOut.fragCodes[idxFirst][p1][1];
          r2 = fragsOut.fragCodes[idxSecond][p2][1];
        }
      }
    }
  }

  // null when the two molecules share no common core (or one has no fragmentations)
  return core === null ? null : [core, r1!, r2!];
}

/**
 * Where {@link fillRules} finds a fragment in `smilesFrags` (its first index) and a rule by its fragments (`sr1`, then
 * `sr2`: the last rule with them), kept as the rules are built instead of searching the arrays for every pair.
 */
export interface MmpRulesIndex {
  frags: Map<number, number>;
  rules: Map<number, Map<number, number>>;
}

export function newMmpRulesIndex(): MmpRulesIndex {
  return {frags: new Map(), rules: new Map()};
}

function addRule(mmpRules: MmpRules, index: MmpRulesIndex, sr1: number, sr2: number,
  pair: {fs: number, ss: number, core: number}): void {
  mmpRules.rules.push({sr1, sr2, pairs: [pair]});
  let bySr2 = index.rules.get(sr1);
  if (!bySr2)
    index.rules.set(sr1, bySr2 = new Map());
  bySr2.set(sr2, mmpRules.rules.length - 1);
}

function addFrag(mmpRules: MmpRules, index: MmpRulesIndex, frag: number): number {
  mmpRules.smilesFrags.push(frag);
  const at = mmpRules.smilesFrags.length - 1;
  if (!index.frags.has(frag))
    index.frags.set(frag, at);
  return at;
}

function fillRules(
  mmpRules: MmpRules, index: MmpRulesIndex,
  rFirst:number, idxFirst: number,
  rSecond: number, idxSecond: number,
  core: number, ruleCounter: number, allCasesCounter: number) {
  let ruleSmiles1 = index.frags.get(rFirst) ?? -1;
  let ruleSmiles2 = index.frags.get(rSecond) ?? -1;
  const ruleIndexStraight = index.rules.get(ruleSmiles1)?.get(ruleSmiles2) ?? -1;
  const ruleIndexInverse = index.rules.get(ruleSmiles2)?.get(ruleSmiles1) ?? -1;

  if (ruleSmiles1 === -1)
    ruleSmiles1 = addFrag(mmpRules, index, rFirst);
  if (ruleSmiles2 === -1)
    ruleSmiles2 = addFrag(mmpRules, index, rSecond);

  const indxFirst = ruleSmiles1 < ruleSmiles2;
  if (ruleIndexStraight === -1) {
    addRule(mmpRules, index, indxFirst ? ruleSmiles1: ruleSmiles2, indxFirst ? ruleSmiles2: ruleSmiles1,
      {fs: indxFirst ? idxFirst : idxSecond,
        ss: indxFirst ? idxSecond : idxFirst,
        core});
    ruleCounter++;
    addRule(mmpRules, index, indxFirst ? ruleSmiles2: ruleSmiles1, indxFirst ? ruleSmiles1: ruleSmiles2,
      {fs: indxFirst ? idxSecond : idxFirst,
        ss: indxFirst ? idxFirst : idxSecond,
        core});
    ruleCounter++;
    allCasesCounter += 2;
  } else {
    mmpRules.rules[ruleIndexStraight].pairs
      .push({fs: idxFirst, ss: idxSecond, core});
    mmpRules.rules[ruleIndexInverse].pairs
      .push({fs: idxSecond, ss: idxFirst, core});
    allCasesCounter += 2;
  }

  return [ruleCounter, allCasesCounter];
}

function encodeFragments(fragsOut: IMmpFragmentsResult): MmpFragments {
  const res: [number, number][][] = new Array(fragsOut.frags.length)
    .fill(null).map((_, i) => new Array(fragsOut.frags[i].length).fill(null).map((_) => [0, 0]));

  let fragmetnCounter = 0;
  const fragmentMap: Record<string, number> = {};
  for (let i = 0; i < fragsOut.frags.length; i++) {
    for (let j = 0; j < fragsOut.frags[i].length; j++) {
      if (fragmentMap[fragsOut.frags[i][j][0]] == undefined) {
        fragmentMap[fragsOut.frags[i][j][0]] = fragmetnCounter;
        fragmetnCounter++;
      }
      if (fragmentMap[fragsOut.frags[i][j][1]] == undefined) {
        fragmentMap[fragsOut.frags[i][j][1]] = fragmetnCounter;
        fragmetnCounter++;
      }
      res[i][j][0] = fragmentMap[fragsOut.frags[i][j][0]];
      res[i][j][1] = fragmentMap[fragsOut.frags[i][j][1]];
    }
  }

  const maxFragmentIndex = fragmetnCounter;
  const fragIdToFragName = new Array<string>(maxFragmentIndex);
  Object.entries(fragmentMap).forEach(([key, val]) => {
    fragIdToFragName[val] = key;
  });

  const fragSizes = new Uint32Array(maxFragmentIndex).map((_, i) => fragIdToFragName[i].length);

  return {fragCodes: res, idToName: fragIdToFragName, sizes: fragSizes};
}

export function getMmpRulesCPU(fragsOut: MmpFragments, fragmentCutoff: number): [MmpRules, number, boolean] {
  const mmpRules: MmpRules = {rules: [], smilesFrags: []};
  const index = newMmpRulesIndex();
  const dim = fragsOut.fragCodes.length;

  const pairs: MolecularPair[] = [];

  for (let i = 0; i < dim; i++) {
    const dimFirstMolecule = fragsOut.fragCodes[i].length;
    for (let j = i + 1; j < dim; j++) {
      const dimSecondMolecule = fragsOut.fragCodes[j].length;
      const best = getBestFragmentPair(dimFirstMolecule, i, dimSecondMolecule, j, fragsOut);
      if (best === null)
        continue;
      const [core, firstR, secondR] = best;
      if (fragsOut.idToName[core] === '' ||
        fragsOut.sizes[firstR] / fragsOut.sizes[core] > fragmentCutoff ||
        fragsOut.sizes[secondR] / fragsOut.sizes[core] > fragmentCutoff)
        continue;

      pairs.push({first: i, second: j, core, firstR, secondR});
    }
  }

  let ruleCounter = 0;
  let allCasesCounter = 0;

  for (let i = 0; i < pairs.length; i++) {
    [ruleCounter, allCasesCounter] = fillRules(mmpRules, index,
      pairs[i].firstR, pairs[i].first, pairs[i].secondR, pairs[i].second, pairs[i].core,
      ruleCounter, allCasesCounter);
  }

  return [mmpRules, allCasesCounter, false];
}

async function getMmpRulesGPU(fragsOut: MmpFragments, fragmentCutoff: number):
  Promise<[MmpRules, number, boolean]> {
  //const [encodedFrags, fragIdToFragName, fragSizes] = encodeFragments(fragsOut);
  const pairs = await webGPUMMP(fragsOut.fragCodes, fragsOut.sizes, fragmentCutoff);
  const [mmpRules, allCasesCounter] = mmpRulesOfPairs(pairs!.frag1Idx, pairs!.mol1Idx, pairs!.frag2Idx,
    pairs!.mol2Idx, pairs!.coreIdx);
  return [mmpRules, allCasesCounter, true];
}

/** The rules of the pairs in order, pair `i` the fragments `frag1[i]` of molecule `mol1[i]` and `frag2[i]` of
 * `mol2[i]` around the fragment `core[i]` (the GPU's pairs), and the number of cases. */
export function mmpRulesOfPairs(frag1: ArrayLike<number>, mol1: ArrayLike<number>, frag2: ArrayLike<number>,
  mol2: ArrayLike<number>, core: ArrayLike<number>): [MmpRules, number] {
  const mmpRules: MmpRules = {rules: [], smilesFrags: []};
  const index = newMmpRulesIndex();
  let ruleCounter = 0;
  let allCasesCounter = 0;
  for (let i = 0; i < core.length; i++) {
    [ruleCounter, allCasesCounter] = fillRules(mmpRules, index, frag1[i], mol1[i], frag2[i], mol2[i], core[i],
      ruleCounter, allCasesCounter);
  }
  return [mmpRules, allCasesCounter];
}
