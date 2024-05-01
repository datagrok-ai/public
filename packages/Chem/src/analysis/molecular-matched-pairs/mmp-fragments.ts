import * as DG from 'datagrok-api/dg';
import {getRdKitService} from '../../utils/chem-common-rdkit';

export type MmpRules = {
  rules: {
    smilesRule1: number,
    smilesRule2: number,
    pairs: {firstStructure: number, secondStructure: number}[]
  } [],
  smilesFrags: string[]
};

type Fragments = [string, string][][];

/**
* Runs paralled fragmentation for molecules
* @param {DG.Column} molecules column with molecules
*/
export async function getMmpFrags(molecules: DG.Column): Promise<Fragments> {
  const service = await getRdKitService();
  const res = await service.mmpGetFragments(molecules.toList());
  return res;
}

function getBestFragmentPair(
  dimFirst: number, idxFirst: number,
  dimSecond: number, idxSecond: number,
  frags: Fragments): [string, string, string] {
  let core = '';
  let r1 = ''; // molecule minus core for first molecule in pair
  let r2 = ''; // molecule minus core for second molecule in pair

  //here we get the best possible fragment pair
  for (let p1 = 0; p1 < dimFirst; p1++) {
    for (let p2 = 0; p2 < dimSecond; p2++) {
      if (frags[idxFirst][p1][0] === frags[idxSecond][p2][0]) {
        const newCore = frags[idxFirst][p1][0];
        if (newCore.length > core.length) {
          core = newCore;
          r1 = frags[idxFirst][p1][1];
          r2 = frags[idxSecond][p2][1];
        }
      }
    }
  }

  return [core, r1, r2];
}

function fillRules(
  mmpRules: MmpRules,
  rFirst:string, idxFirst: number,
  rSecond: string, idxSecond: number,
  ruleCounter: number, allCasesCounter: number) {
  let ruleSmiles1 = mmpRules.smilesFrags.indexOf(rFirst);
  let ruleSmiles2 = mmpRules.smilesFrags.indexOf(rSecond);
  let ruleIndexStraight = -1;
  let ruleIndexInverse = -1;

  for (let ind = 0; ind < mmpRules.rules.length; ind++) {
    if (mmpRules.rules[ind].smilesRule1 == ruleSmiles1 && mmpRules.rules[ind].smilesRule2 == ruleSmiles2)
      ruleIndexStraight = ind;
    if (mmpRules.rules[ind].smilesRule1 == ruleSmiles2 && mmpRules.rules[ind].smilesRule2 == ruleSmiles1)
      ruleIndexInverse = ind;
  }

  if (ruleSmiles1 == -1) {
    mmpRules.smilesFrags.push(rFirst);
    ruleSmiles1 = mmpRules.smilesFrags.length -1;
  }
  if (ruleSmiles2 == -1) {
    mmpRules.smilesFrags.push(rSecond);
    ruleSmiles2 = mmpRules.smilesFrags.length -1;
  }

  const indxFirst = ruleSmiles1 < ruleSmiles2;
  if (ruleIndexStraight == -1) {
    mmpRules.rules.push({
      smilesRule1: indxFirst ? ruleSmiles1: ruleSmiles2,
      smilesRule2: indxFirst ? ruleSmiles2: ruleSmiles1,
      pairs: [],
    });
    mmpRules.rules[ruleCounter].pairs
      .push({firstStructure: indxFirst ? idxFirst : idxSecond, secondStructure: indxFirst ? idxSecond : idxFirst});
    ruleCounter++;
    mmpRules.rules.push({
      smilesRule1: indxFirst ? ruleSmiles2: ruleSmiles1,
      smilesRule2: indxFirst ? ruleSmiles1: ruleSmiles2,
      pairs: [],
    });
    mmpRules.rules[ruleCounter].pairs
      .push({firstStructure: indxFirst ? idxSecond : idxFirst, secondStructure: indxFirst ? idxFirst : idxSecond});
    ruleCounter++;
    allCasesCounter += 2;
  } else {
    mmpRules.rules[ruleIndexStraight].pairs
      .push({firstStructure: idxFirst, secondStructure: idxSecond});
    mmpRules.rules[ruleIndexInverse].pairs
      .push({firstStructure: idxSecond, secondStructure: idxFirst});
    allCasesCounter += 2;
  }

  return [ruleCounter, allCasesCounter];
}

export function getMmpRules(frags: Fragments, fragmentCutoff: number): [MmpRules, number] {
  const mmpRules: MmpRules = {rules: [], smilesFrags: []};
  const dim = frags.length;

  const nelements = dim*(dim - 1) / 2; //number of elements in upper triangle
  const cores = new Array<string>(nelements);
  const firstRs = new Array<string>(nelements);
  const secondRs = new Array<string>(nelements);

  for (let i = 0; i < dim; i++) {
    const dimFirstMolecule = frags[i].length;
    const rowPart = i*dim- i*(i + 1) / 2 - i - 1; //number of elements passed in upper triangele
    for (let j = i + 1; j < dim; j++) {
      const idx = rowPart + j;

      const dimSecondMolecule = frags[j].length;
      [cores[idx], firstRs[idx], secondRs[idx]] = getBestFragmentPair(dimFirstMolecule, i, dimSecondMolecule, j, frags);
    }
  }

  let ruleCounter = 0;
  let allCasesCounter = 0;

  let idx = 0;
  for (let i = 0; i < dim; i++) {
    for (let j = i + 1; j < dim; j++) {
      if (cores[idx] === '' ||
          firstRs[idx].length / cores[idx].length > fragmentCutoff ||
          secondRs[idx].length / cores[idx].length > fragmentCutoff) {idx++; continue;}

      [ruleCounter, allCasesCounter] =
        fillRules(mmpRules, firstRs[idx], i, secondRs[idx], j, ruleCounter, allCasesCounter);
      idx++;
    }
  }

  return [mmpRules, allCasesCounter];
}
