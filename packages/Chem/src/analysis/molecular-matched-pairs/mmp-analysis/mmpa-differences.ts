import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {MmpAllCasesBasedData, MmpFragments, MmpInitData, MmpRules, MmpRulesBasedData} from './mmpa-misc';
import {SortData} from '../mmp-viewer/mmp-viewer';

export function getPlainData(rules: MmpRules, frags: MmpFragments,
  initData: MmpInitData, allCasesNumber: number, fragSortingInfo: SortData): [MmpRulesBasedData, MmpAllCasesBasedData] {
  const [fromFrag, toFrag, occasions] = getAllRulesOcasions(rules, frags, fragSortingInfo); //rules n objects
  const [maxActs, meanDiffs, molFrom, molTo, pairNum,
    molNumFrom, molNumTo, pairsFromSmiles, pairsToSmiles,
    ruleNum, diffs, activityPairsIdxs] =
    calculateActivityDiffs(initData.molecules, initData.activities, rules, frags, allCasesNumber, occasions);

  const rulesBased: MmpRulesBasedData = {fromFrag, toFrag, occasions, meanDiffs};
  const allCasesBased: MmpAllCasesBasedData = {maxActs, molFrom, molTo, pairNum,
    molNumFrom, molNumTo, pairsFromSmiles, pairsToSmiles, ruleNum, diffs, activityPairsIdxs};

  return [rulesBased, allCasesBased];
}

function getAllRulesOcasions(mmpr: MmpRules, frags: MmpFragments, fragSortingInfo: SortData):
[string [], string [], Int32Array] {
  const allSize = mmpr.rules.length;
  const fromFrag = new Array<string>(allSize);
  const toFrag = new Array<string>(allSize);
  const occasions = new Int32Array(allSize);
  for (let i = 0; i < allSize; i++) {
    const fromFragment = mmpr.smilesFrags[mmpr.rules[i].smilesRule1];
    fromFrag[i] = frags.idToName[fromFragment];
    toFrag[i] = frags.idToName[mmpr.smilesFrags[mmpr.rules[i].smilesRule2]];
    occasions[i] = mmpr.rules[i].pairs.length;
    const fragIdxInSortingInfo = fragSortingInfo.fragmentIdxs.indexOf(fromFragment);
    if (fragIdxInSortingInfo === -1) {
      fragSortingInfo.fragmentIdxs.push(fromFragment);
      fragSortingInfo.frequencies.push(mmpr.rules[i].pairs.length);
    } else
      fragSortingInfo.frequencies[fragIdxInSortingInfo] += mmpr.rules[i].pairs.length;
  }
  return [fromFrag, toFrag, occasions];
}

function calculateActivityDiffs(
  molecules: string[],
  activities: Float32Array[],
  mmpr: MmpRules, frags: MmpFragments, allCasesNumber: number,
  occasions: Int32Array) :
  [ maxActs: number [],
    meanDiffs: Float32Array[],
    molFrom: string[],
    molTo: string[],
    pairNum: Int32Array,
    molNumFrom: Int32Array,
    molNumTo: Int32Array,
    pairsFromSmiles: string[],
    pairsToSmiles: string[],
    ruleNum: Int32Array,
    diffs: Float32Array[],
    activityPairsIdxs: BitArray[]
 ] {
  const variates = activities.length;
  const maxActs = new Array<number>(variates).fill(0);
  const allSize = mmpr.rules.length;

  const meanDiffs = new Array<Float32Array>(variates);
  for (let i = 0; i < variates; i++)
    meanDiffs[i] = new Float32Array(allSize);

  const molFrom = new Array<string>(allCasesNumber);
  const molTo = new Array<string>(allCasesNumber);
  const pairNum = new Int32Array(allCasesNumber);
  const molNumFrom = new Int32Array(allCasesNumber);
  const molNumTo = new Int32Array(allCasesNumber);
  const diffs = new Array<Float32Array>(variates);
  for (let i = 0; i < variates; i++)
    diffs[i] = new Float32Array(allCasesNumber);

  const pairsFromSmiles = new Array<string>(allCasesNumber);
  const pairsToSmiles = new Array<string>(allCasesNumber);
  const ruleNum = new Int32Array(allCasesNumber);

  const activityPairsIdxs = new Array<BitArray>(variates);
  for (let i = 0; i < variates; i++)
    activityPairsIdxs[i] = new BitArray(allCasesNumber);

  let pairIdx = 0;
  for (let i = 0; i < allSize; i++) {
    const mean = new Float32Array(variates);
    for (let j = 0; j < occasions[i]; j++) {
      const idx1 = mmpr.rules[i].pairs[j].firstStructure;
      const idx2 = mmpr.rules[i].pairs[j].secondStructure;

      molFrom[pairIdx] = molecules[idx1];
      molTo[pairIdx] = molecules[idx2];

      for (let k = 0; k < variates; k++) {
        //TODO: make more efficient
        const col = activities[k];
        const diff = col[idx2] - col[idx1];
        diffs[k][pairIdx] = diff;
        if (diff > 0) {
          if (diff > maxActs[k])
            maxActs[k] = diff;
          activityPairsIdxs[k].setBit(pairIdx, true, false);
        }

        mean[k] += diff;
      }

      molNumFrom[pairIdx] = idx1;
      molNumTo[pairIdx] = idx2;
      pairNum[pairIdx] = pairIdx;
      pairsFromSmiles[pairIdx] = frags.idToName[mmpr.smilesFrags[mmpr.rules[i].smilesRule1]];
      pairsToSmiles[pairIdx] = frags.idToName[mmpr.smilesFrags[mmpr.rules[i].smilesRule2]];
      ruleNum[pairIdx] = i;

      pairIdx++;
    }

    for (let k = 0; k < variates; k++) {
      mean[k] /= occasions[i];
      meanDiffs[k][i] = mean[k];
    }
  }

  return [maxActs, meanDiffs, molFrom, molTo,
    pairNum, molNumFrom, molNumTo, pairsFromSmiles, pairsToSmiles, ruleNum, diffs, activityPairsIdxs];
}
