import * as DG from 'datagrok-api/dg';
import {getRdKitService} from '../../utils/chem-common-rdkit';
import {webGPUMMP} from '@datagrok-libraries/math/src/webGPU/mmp/webGPU-mmp';

export type MmpRules = {
  rules: {
    smilesRule1: number,
    smilesRule2: number,
    pairs: {firstStructure: number, secondStructure: number}[]
  } [],
  smilesFrags: string[]
};

type Fragments = [string, string][][];

type MolecularPair = {
  first: number,
  second: number,
  core: string,
  firstR: string,
  secondR: string
}

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

  //const nelements = dim*(dim - 1) / 2; //number of elements in upper triangle

  const pairs: MolecularPair[] = [];

  for (let i = 0; i < dim; i++) {
    const dimFirstMolecule = frags[i].length;
    //const rowPart = i*dim- i*(i + 1) / 2 - i - 1; //number of elements passed in upper triangle
    for (let j = i + 1; j < dim; j++) {
      //const idx = rowPart + j;
      const dimSecondMolecule = frags[j].length;
      //[cores[idx], firstRs[idx], secondRs[idx]] =
      //getBestFragmentPair(dimFirstMolecule, i, dimSecondMolecule, j, frags);
      const [core, firstR, secondR] = getBestFragmentPair(dimFirstMolecule, i, dimSecondMolecule, j, frags);
      if (core === '' ||
      firstR.length / core.length > fragmentCutoff ||
      secondR.length / core.length > fragmentCutoff) {/*idx++;*/ continue;}

      pairs.push({first: i, second: j, core, firstR, secondR});
    }
  }

  let ruleCounter = 0;
  let allCasesCounter = 0;

  for (let i = 0; i < pairs.length; i++) {
    [ruleCounter, allCasesCounter] = fillRules(mmpRules,
      pairs[i].firstR, pairs[i].first, pairs[i].secondR, pairs[i].second,
      ruleCounter, allCasesCounter);
  }

  return [mmpRules, allCasesCounter];
}

//From full matrix i and j adresses to k in linear upper triangle k = i*dim + j - i*(i + 1) / 2 - i - 1;


export async function getMmpRulesP(frags: Fragments, fragmentCutoff: number): Promise<[MmpRules, number]> {
  const mmpRules: MmpRules = {rules: [], smilesFrags: []};
  const dim = frags.length;

  const nelements = dim*(dim - 1) / 2; //number of elements in upper triangle
  const pairIndexes = new Array<[number, number]>(nelements);

  let k = 0;
  for (let i = 0; i < dim; i++) {
    for (let j = i + 1; j < dim; j++) {
      pairIndexes[k] = [i, j]; k++;
    }
  }

  const service = await getRdKitService();
  //const pairs: [number, number, string, string, string][] = [];
  const pairs = await service.mmpGetPairs(pairIndexes, frags, fragmentCutoff);

  // for (let i = 0; i < nelements; i++) {
  //   const dimFirstMolecule = frags[pairIndexes[i][0]].length;
  //   const dimSecondMolecule = frags[pairIndexes[i][1]].length;
  //   const [core, firstR, secondR] =
  //     getBestFragmentPair(dimFirstMolecule, pairIndexes[i][0], dimSecondMolecule, pairIndexes[i][1], frags);
  //   if (core === '' ||
  //   firstR.length / core.length > fragmentCutoff ||
  //   secondR.length / core.length > fragmentCutoff)
  //     continue;

  //   pairs.push([pairIndexes[i][0], pairIndexes[i][1], core, firstR, secondR]);
  // }

  let ruleCounter = 0;
  let allCasesCounter = 0;

  for (let i = 0; i < pairs.length; i++) {
    [ruleCounter, allCasesCounter] = fillRules(mmpRules,
      pairs[i][3], pairs[i][0], pairs[i][4], pairs[i][1],
      ruleCounter, allCasesCounter);
  }

  return [mmpRules, allCasesCounter];
}

export async function getMmpRulesPP(frags: Fragments, fragmentCutoff: number): Promise<[MmpRules, number]> {
  const mmpRules: MmpRules = {rules: [], smilesFrags: []};

  const dim = frags.length;

  //const nelements = dim*(dim - 1) / 2; //number of elements in upper triangle

  const pairs: MolecularPair[] = [];

  // for (let i = 0; i < dim; i++) {
  //   const dimFirstMolecule = frags[i].length;
  //   //const rowPart = i*dim- i*(i + 1) / 2 - i - 1; //number of elements passed in upper triangle
  //   for (let j = i + 1; j < dim; j++) {
  //     //const idx = rowPart + j;
  //     const dimSecondMolecule = frags[j].length;
  //     //[cores[idx], firstRs[idx], secondRs[idx]] =
  //     //getBestFragmentPair(dimFirstMolecule, i, dimSecondMolecule, j, frags);
  //     const [core, firstR, secondR] = getBestFragmentPair(dimFirstMolecule, i, dimSecondMolecule, j, frags);
  //     if (core === '' ||
  //     firstR.length / core.length > fragmentCutoff ||
  //     secondR.length / core.length > fragmentCutoff) {/*idx++;*/ continue;}

  //     pairs.push({first: i, second: j, core, firstR, secondR});
  //   }
  // }


  const {res: encodedFrags, fragmentMap, fragIdToFragName, fragSizes} = encodeFragments(frags);
  const out = await webGPUMMP(encodedFrags, fragSizes, fragmentCutoff);
  return [mmpRules, 100];
}


function encodeFragments(frags: Fragments) {
  const res: [number, number][][] = new Array(frags.length)
    .fill(null).map((_, i) => new Array(frags[i].length).fill(null).map((_) => [0, 0]));

  let fragmetnCounter = 0;
  const fragmentMap: Record<string, number> = {};
  for (let i = 0; i < frags.length; i++) {
    for (let j = 0; j < frags[i].length; j++) {
      if (!fragmentMap[frags[i][j][0]]) {
        fragmentMap[frags[i][j][0]] = fragmetnCounter;
        fragmetnCounter++;
      }
      if (!fragmentMap[frags[i][j][1]]) {
        fragmentMap[frags[i][j][1]] = fragmetnCounter;
        fragmetnCounter++;
      }
      res[i][j][0] = fragmentMap[frags[i][j][0]];
      res[i][j][1] = fragmentMap[frags[i][j][1]];
    }
  }

  const maxFragmentIndex = fragmetnCounter;
  const fragIdToFragName = new Array<string>(maxFragmentIndex);
  Object.entries(fragmentMap).forEach(([key, val]) => {
    fragIdToFragName[val] = key;
  });

  const fragSizes = new Uint32Array(maxFragmentIndex).map((_, i) => fragIdToFragName[i].length);


  return {res, fragmentMap, fragIdToFragName, fragSizes};
}
