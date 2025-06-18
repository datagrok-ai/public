import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {MmpAllCasesBasedData, MmpFragments, MmpInitData, MmpRules, MmpRulesBasedData} from './mmpa-misc';
import {SortData} from '../mmp-viewer/mmp-viewer';
import * as DG from 'datagrok-api/dg';
import {MmpDiffTypes} from '../mmp-function-editor';

export function getPlainData(rules: MmpRules, frags: MmpFragments,
  initData: MmpInitData, allCasesNumber: number, fragSortingInfo: SortData): [MmpRulesBasedData, MmpAllCasesBasedData] {
  const [fromFrag, toFrag, occasions] = getAllRulesOcasions(rules, frags, fragSortingInfo);
  const [maxActs, meanDiffs, molFrom, molTo, pairNum,
    molNumFrom, molNumTo, pairsFromSmiles, pairsToSmiles,
    ruleNum, diffs, activityPairsIdxs, coreNums, pairedActivities] =
    calculateActivityDiffs(initData.molecules, initData.activities, initData.diffTypes,
      rules, frags, allCasesNumber, occasions);

  // fragmentParentIndices - for each fragment pair there is an array of
  // parent molecule indexes where this fragment was identified
  // moleculePairParentIndices - indexes of 'from' parent molecules and 'to' parent molecules
  // moleculePairParentIndices - indexes of pairs from molecule pairs dataset for each fragment pair
  const fragmentParentIndices: Array<Uint32Array> = new Array<Uint32Array>(rules.rules.length);
  const moleculePairParentIndices: [Uint32Array, Uint32Array] = [
    new Uint32Array(allCasesNumber),
    new Uint32Array(allCasesNumber),
  ];
  const moleculePairIndices: Array<Uint32Array> = new Array<Uint32Array>(rules.rules.length);
  let pairIdx = 0;
  for (let i = 0; i < rules.rules.length; i++) {
    const rule = rules.rules[i];
    const indices = new Uint32Array(rule.pairs.length * 2);
    const pairsIndices = new Uint32Array(rule.pairs.length);
    for (let j = 0; j < rule.pairs.length; j++) {
      const pair = rule.pairs[j];
      // Set fragment parent indices
      indices[j * 2] = pair.fs;
      indices[j * 2 + 1] = pair.ss;
      // Set molecule pair parent indices
      moleculePairParentIndices[0][pairIdx] = pair.fs;
      moleculePairParentIndices[1][pairIdx] = pair.ss;
      pairsIndices[j] = pairIdx;
      pairIdx++;
    }
    fragmentParentIndices[i] = indices;
    moleculePairIndices[i] = pairsIndices;
  }

  const rulesBased: MmpRulesBasedData = {
    fromFrag,
    toFrag,
    occasions,
    meanDiffs,
    fragmentParentIndices,
    moleculePairParentIndices,
    moleculePairIndices,
  };
  const allCasesBased: MmpAllCasesBasedData = {
    maxActs,
    molFrom,
    molTo,
    pairNum,
    molNumFrom,
    molNumTo,
    pairsFromSmiles,
    pairsToSmiles,
    ruleNum,
    diffs,
    activityPairsIdxs,
    coreNums,
    pairedActivities,
  };

  return [rulesBased, allCasesBased];
}

function getAllRulesOcasions(mmpr: MmpRules, frags: MmpFragments, fragSortingInfo: SortData):
[string [], string [], Int32Array] {
  const allSize = mmpr.rules.length;
  const fromFrag = new Array<string>(allSize);
  const toFrag = new Array<string>(allSize);
  const occasions = new Int32Array(allSize);
  for (let i = 0; i < allSize; i++) {
    const fromFragment = mmpr.smilesFrags[mmpr.rules[i].sr1];
    fromFrag[i] = frags.idToName[fromFragment];
    toFrag[i] = frags.idToName[mmpr.smilesFrags[mmpr.rules[i].sr2]];
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

const delta = (x: number, y: number) => {return x - y;};
const ratio = (x: number, y: number) => {return x / y;};

function calculateActivityDiffs(
  molecules: string[],
  activities: Float32Array[],
  diffTypes: string[],
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
    activityPairsIdxs: BitArray[],
    coreNums: Int32Array,
    pairedActivities: Array<[Float32Array, Float32Array]>,
 ] {
  const variates = activities.length;
  const maxActs = new Array<number>(variates).fill(0);
  const allSize = mmpr.rules.length;

  const meanDiffs = new Array<Float32Array>(variates);
  const pairedActivities = new Array<[Float32Array, Float32Array]>(variates);
  for (let i = 0; i < variates; i++)
    meanDiffs[i] = new Float32Array(allSize);

  const molFrom = new Array<string>(allCasesNumber);
  const molTo = new Array<string>(allCasesNumber);
  const coresNums = new Int32Array(allCasesNumber);
  const pairNum = new Int32Array(allCasesNumber);
  const molNumFrom = new Int32Array(allCasesNumber);
  const molNumTo = new Int32Array(allCasesNumber);
  const diffs = new Array<Float32Array>(variates);
  for (let i = 0; i < variates; i++) {
    pairedActivities[i] = [new Float32Array(allCasesNumber), new Float32Array(allCasesNumber)];
    diffs[i] = new Float32Array(allCasesNumber);
  }

  const pairsFromSmiles = new Array<string>(allCasesNumber);
  const pairsToSmiles = new Array<string>(allCasesNumber);
  const ruleNum = new Int32Array(allCasesNumber);

  const activityPairsIdxs = new Array<BitArray>(variates);
  for (let i = 0; i < variates; i++)
    activityPairsIdxs[i] = new BitArray(allCasesNumber);

  let pairIdx = 0;
  for (let i = 0; i < allSize; i++) {
    const mean = new Float32Array(variates).fill(0);
    const validCounts = new Int32Array(variates).fill(0);
    for (let j = 0; j < occasions[i]; j++) {
      mmpr.rules[i].pairs[j].id = pairIdx;
      const idx1 = mmpr.rules[i].pairs[j].fs;
      const idx2 = mmpr.rules[i].pairs[j].ss;

      molFrom[pairIdx] = molecules[idx1];
      molTo[pairIdx] = molecules[idx2];
      coresNums[pairIdx] = mmpr.rules[i].pairs[j].core;

      for (let k = 0; k < variates; k++) {
        const deltaFunc = diffTypes[k] as MmpDiffTypes === MmpDiffTypes.delta ? delta :
          diffTypes[k] as MmpDiffTypes === MmpDiffTypes.ratio ? ratio : delta;
        //TODO: make more efficient
        const col = activities[k];
        pairedActivities[k][0][pairIdx] = col[idx1];
        pairedActivities[k][1][pairIdx] = col[idx2];
        if (col[idx2] == DG.FLOAT_NULL || col[idx1] == DG.FLOAT_NULL) {
          diffs[k][pairIdx] = DG.FLOAT_NULL;
          continue;
        }
        const diff = deltaFunc(col[idx2], col[idx1]);
        diffs[k][pairIdx] = diff;
        if (diff > 0) {
          if (diff > maxActs[k])
            maxActs[k] = diff;
          activityPairsIdxs[k].setBit(pairIdx, true, false);
        }

        mean[k] += diff;
        validCounts[k]++;
      }

      molNumFrom[pairIdx] = idx1;
      molNumTo[pairIdx] = idx2;
      pairNum[pairIdx] = pairIdx;
      pairsFromSmiles[pairIdx] = frags.idToName[mmpr.smilesFrags[mmpr.rules[i].sr1]];
      pairsToSmiles[pairIdx] = frags.idToName[mmpr.smilesFrags[mmpr.rules[i].sr2]];
      ruleNum[pairIdx] = i;

      pairIdx++;
    }

    for (let k = 0; k < variates; k++)
      meanDiffs[k][i] = validCounts[k] > 0 ?mean[k] / validCounts[k] : DG.FLOAT_NULL;
  }

  return [maxActs, meanDiffs, molFrom, molTo, pairNum, molNumFrom, molNumTo, pairsFromSmiles,
    pairsToSmiles, ruleNum, diffs, activityPairsIdxs, coresNums, pairedActivities];
}
