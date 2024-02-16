import * as DG from 'datagrok-api/dg';
import {getRdKitService} from '../../utils/chem-common-rdkit';

const CUTOFF_MOL_SIZE = 0.4;

export type MmpRules = {
  rules: {
    smilesRule1: number,
    smilesRule2: number,
    pairs: {firstStructure: number, secondStructure: number}[]
  } [],
  smilesFrags: string[]
};

/**
* Runs paralled fragmentation for molecules
* @param {DG.Column} molecules column with molecules
*/
export async function getMmpFrags(molecules: DG.Column): Promise<[string, string][][]> {
  const service = await getRdKitService();
  const res = await service.mmpGetFragments(molecules.toList());
  return res;
}

export function getMmpRules(frags: [string, string][][]): [MmpRules, number] {
  const mmpRules: MmpRules = {rules: [], smilesFrags: []};
  const dim = frags.length;
  let ruleCounter = 0;
  let allCasesCounter = 0;

  for (let i = 0; i < dim; i++) {
    const dimFirstMolecule = frags[i].length;
    for (let j = i + 1; j < dim; j++) {
      const dimSecondMolecule = frags[j].length;
      let core = '';
      let r1 = ''; // molecule minus core for first molecule in pair
      let r2 = ''; // molecule minus core for second molecule in pair

      //here we get the best possible fragment pair
      for (let p1 = 0; p1 < dimFirstMolecule; p1++) {
        for (let p2 = 0; p2 < dimSecondMolecule; p2++) {
          if (frags[i][p1][0] === frags[j][p2][0]) {
            const newCore = frags[i][p1][0];
            if (newCore.length > core.length) {
              core = newCore;
              r1 = frags[i][p1][1];
              r2 = frags[j][p2][1];
            }
          }
        }
      }

      if (core === '' || r1.length / core.length > CUTOFF_MOL_SIZE || r2.length / core.length > CUTOFF_MOL_SIZE)
        continue;

      let ruleSmiles1 = mmpRules.smilesFrags.indexOf(r1);
      let ruleSmiles2 = mmpRules.smilesFrags.indexOf(r2);
      let ruleIndexStraight = -1;
      let ruleIndexInverse = -1;

      for (let ind = 0; ind < mmpRules.rules.length; ind++) {
        if (mmpRules.rules[ind].smilesRule1 == ruleSmiles1 && mmpRules.rules[ind].smilesRule2 == ruleSmiles2)
          ruleIndexStraight = ind;
        if (mmpRules.rules[ind].smilesRule1 == ruleSmiles2 && mmpRules.rules[ind].smilesRule2 == ruleSmiles1)
          ruleIndexInverse = ind;
      }

      if (ruleSmiles1 == -1) {
        mmpRules.smilesFrags.push(r1);
        ruleSmiles1 = mmpRules.smilesFrags.length -1;
      }
      if (ruleSmiles2 == -1) {
        mmpRules.smilesFrags.push(r2);
        ruleSmiles2 = mmpRules.smilesFrags.length -1;
      }

      const indxFirst = ruleSmiles1 < ruleSmiles2;
      if (ruleIndexStraight == -1) {
        mmpRules.rules.push({
          smilesRule1: indxFirst ? ruleSmiles1: ruleSmiles2,
          smilesRule2: indxFirst ? ruleSmiles2: ruleSmiles1,
          pairs: [],
        });
        mmpRules.rules[ruleCounter].pairs.push({firstStructure: indxFirst ? i : j, secondStructure: indxFirst ? j : i});
        ruleCounter++;
        mmpRules.rules.push({
          smilesRule1: indxFirst ? ruleSmiles2: ruleSmiles1,
          smilesRule2: indxFirst ? ruleSmiles1: ruleSmiles2,
          pairs: [],
        });
        mmpRules.rules[ruleCounter].pairs.push({firstStructure: indxFirst ? j : i, secondStructure: indxFirst ? i : j});
        ruleCounter++;
        allCasesCounter += 2;
      } else {
        mmpRules.rules[ruleIndexStraight].pairs
          .push({firstStructure: i, secondStructure: j});
        mmpRules.rules[ruleIndexInverse].pairs
          .push({firstStructure: j, secondStructure: i});
        allCasesCounter += 2;
      }
    }
  }

  return [mmpRules, allCasesCounter];
}
