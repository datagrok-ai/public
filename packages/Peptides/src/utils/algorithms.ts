import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import * as C from './constants';
import * as type from './types';
import {getTypedArrayConstructor} from './misc';
import {CLUSTER_TYPE, ClusterStats, ClusterTypeStats, MonomerPositionStats, PositionStats, SummaryStats} from '../model';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {Stats, getStats} from './statistics';


type MutationCliffInfo = {pos: string, seq1monomer: string, seq2monomer: string, seq1Idx: number, seq2Idx: number};

export function findMutations(activityArray: type.RawData, monomerInfoArray: type.RawColumn[],
  settings: type.PeptidesSettings = {},
  targetOptions: {targetCol?: type.RawColumn | null, currentTarget?: string | null} = {}): type.MutationCliffs {
  const nCols = monomerInfoArray.length;
  if (nCols === 0)
    throw new Error(`PepAlgorithmError: Couldn't find any column of semType '${C.SEM_TYPES.MONOMER}'`);

  settings.minActivityDelta ??= 0;
  settings.maxMutations ??= 1;
  const currentTargetIdx = targetOptions.targetCol?.cat!.indexOf(targetOptions.currentTarget!) ?? -1;

  const substitutionsInfo: type.MutationCliffs = new Map();
  const nRows = activityArray.length;
  for (let seq1Idx = 0; seq1Idx < nRows - 1; seq1Idx++) {
    if (currentTargetIdx !== -1 && targetOptions.targetCol?.rawData[seq1Idx] !== currentTargetIdx)
      continue;

    for (let seq2Idx = seq1Idx + 1; seq2Idx < nRows; seq2Idx++) {
      if (currentTargetIdx !== -1 && targetOptions.targetCol?.rawData[seq2Idx] !== currentTargetIdx)
        continue;

      let substCounter = 0;
      const activityValSeq1 = activityArray[seq1Idx];
      const activityValSeq2 = activityArray[seq2Idx];
      const delta = activityValSeq1 - activityValSeq2;
      if (Math.abs(delta) < settings.minActivityDelta)
        continue;

      let substCounterFlag = false;
      const tempData: MutationCliffInfo[] = new Array(monomerInfoArray.length);
      let tempDataIdx = 0;
      for (const monomerInfo of monomerInfoArray) {
        const seq1categoryIdx = monomerInfo.rawData[seq1Idx];
        const seq2categoryIdx = monomerInfo.rawData[seq2Idx];
        if (seq1categoryIdx === seq2categoryIdx)
          continue;

        substCounter++;
        substCounterFlag = substCounter > settings.maxMutations;
        if (substCounterFlag)
          break;

        tempData[tempDataIdx++] ={
          pos: monomerInfo.name,
          seq1monomer: monomerInfo.cat![seq1categoryIdx],
          seq2monomer: monomerInfo.cat![seq2categoryIdx],
          seq1Idx: seq1Idx,
          seq2Idx: seq2Idx,
        };
      }

      if (substCounterFlag || substCounter === 0)
        continue;

      // Separate processing loop in case substCOunter is 0 or out of restricted range to prevent unnecessary computations
      for (let i = 0; i < tempDataIdx; i++) {
        const tempDataElement = tempData[i];
        //Working with seq1monomer
        const seq1monomer = tempDataElement.seq1monomer;
        if (!substitutionsInfo.has(seq1monomer))
          substitutionsInfo.set(seq1monomer, new Map());

        const position = tempDataElement.pos;

        let positionsMap = substitutionsInfo.get(seq1monomer)!;
        if (!positionsMap.has(position))
          positionsMap.set(position, new Map());

        let indexes = positionsMap.get(position)!;
        if (indexes.has(seq1Idx))
          (indexes.get(seq1Idx)! as number[]).push(seq2Idx);
        else
          indexes.set(seq1Idx, [seq2Idx]);

        //Working with seq2monomer
        const seq2monomer = tempDataElement.seq2monomer;
        if (!substitutionsInfo.has(seq2monomer))
          substitutionsInfo.set(seq2monomer, new Map());

        positionsMap = substitutionsInfo.get(seq2monomer)!;
        if (!positionsMap.has(position))
          positionsMap.set(position, new Map());

        indexes = positionsMap.get(position)!;
        if (indexes.has(seq2Idx))
          (indexes.get(seq2Idx)! as number[]).push(seq1Idx);
        else
          indexes.set(seq2Idx, [seq1Idx]);
      }
    }
  }

  const TypedArray = getTypedArrayConstructor(nRows);
  for (const positionMap of substitutionsInfo.values()) {
    for (const indexMap of positionMap.values()) {
      for (const [index, indexArray] of indexMap.entries())
        indexMap.set(index, new TypedArray(indexArray));
    }
  }

  return substitutionsInfo;
}

export function calculateMonomerPositionStatistics(df: DG.DataFrame, positionColumns: DG.Column<string>[],
  options: {isFiltered?: boolean, columns?: string[]} = {}): MonomerPositionStats {
  options.isFiltered ??= false;
  const monomerPositionObject = {general: {}} as MonomerPositionStats & {general: SummaryStats};
  let activityColData: Float64Array = df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED).getRawData() as Float64Array;
  let sourceDfLen = df.rowCount;

  if (options.isFiltered) {
    sourceDfLen = df.filter.trueCount;
    const tempActivityData = new Float64Array(sourceDfLen);
    const selectedIndexes = df.filter.getSelectedIndexes();
    for (let i = 0; i < sourceDfLen; ++i)
      tempActivityData[i] = activityColData[selectedIndexes[i]];
    activityColData = tempActivityData;
    positionColumns = DG.DataFrame.fromColumns(positionColumns).clone(df.filter).columns.toList();
  }
  options.columns ??= positionColumns.map((col) => col.name);

  for (const posCol of positionColumns) {
    if (!options.columns.includes(posCol.name))
      continue;
    const posColData = posCol.getRawData();
    const posColCateogries = posCol.categories;
    const currentPositionObject = {general: {}} as PositionStats & {general: SummaryStats};

    for (let categoryIndex = 0; categoryIndex < posColCateogries.length; ++categoryIndex) {
      const monomer = posColCateogries[categoryIndex];
      if (monomer === '')
        continue;

      const boolArray: boolean[] = new Array(sourceDfLen).fill(false);
      for (let i = 0; i < sourceDfLen; ++i) {
        if (posColData[i] === categoryIndex)
          boolArray[i] = true;
      }
      const bitArray = BitArray.fromValues(boolArray);
      const stats = bitArray.allFalse || bitArray.allTrue ?
        {count: sourceDfLen, meanDifference: 0, ratio: 1.0, pValue: null, mask: bitArray} :
        getStats(activityColData, bitArray);
      currentPositionObject[monomer] = stats;
      getSummaryStats(currentPositionObject.general, stats);
    }
    monomerPositionObject[posCol.name] = currentPositionObject;
    getSummaryStats(monomerPositionObject.general, null, currentPositionObject.general);
  }
  return monomerPositionObject;
}

export function getSummaryStats(genObj: SummaryStats, stats: Stats | null = null, summaryStats: SummaryStats | null = null): void {
  if (stats === null && summaryStats === null)
    throw new Error(`MonomerPositionStatsError: either stats or summaryStats must be present`);

  const possibleMaxCount = stats?.count ?? summaryStats!.maxCount;
  genObj.maxCount ??= possibleMaxCount;
  if (genObj.maxCount < possibleMaxCount)
    genObj.maxCount = possibleMaxCount;

  const possibleMinCount = stats?.count ?? summaryStats!.minCount;
  genObj.minCount ??= possibleMinCount;
  if (genObj.minCount > possibleMinCount)
    genObj.minCount = possibleMinCount;

  const possibleMaxMeanDifference = stats?.meanDifference ?? summaryStats!.maxMeanDifference;
  genObj.maxMeanDifference ??= possibleMaxMeanDifference;
  if (genObj.maxMeanDifference < possibleMaxMeanDifference)
    genObj.maxMeanDifference = possibleMaxMeanDifference;

  const possibleMinMeanDifference = stats?.meanDifference ?? summaryStats!.minMeanDifference;
  genObj.minMeanDifference ??= possibleMinMeanDifference;
  if (genObj.minMeanDifference > possibleMinMeanDifference)
    genObj.minMeanDifference = possibleMinMeanDifference;

  if (!isNaN(stats?.pValue ?? NaN)) {
    const possibleMaxPValue = stats?.pValue ?? summaryStats!.maxPValue;
    genObj.maxPValue ??= possibleMaxPValue;
    if (genObj.maxPValue < possibleMaxPValue)
      genObj.maxPValue = possibleMaxPValue;

    const possibleMinPValue = stats?.pValue ?? summaryStats!.minPValue;
    genObj.minPValue ??= possibleMinPValue;
    if (genObj.minPValue > possibleMinPValue)
      genObj.minPValue = possibleMinPValue;
  }

  const possibleMaxRatio = stats?.ratio ?? summaryStats!.maxRatio;
  genObj.maxRatio ??= possibleMaxRatio;
  if (genObj.maxRatio < possibleMaxRatio)
    genObj.maxRatio = possibleMaxRatio;

  const possibleMinRatio = stats?.ratio ?? summaryStats!.minRatio;
  genObj.minRatio ??= possibleMinRatio;
  if (genObj.minRatio > possibleMinRatio)
    genObj.minRatio = possibleMinRatio;
}

export function calculateClusterStatistics(df: DG.DataFrame, clustersColumnName: string,
  customClusters: DG.Column<boolean>[]): ClusterTypeStats {
  const rowCount = df.rowCount;
  const origClustCol = df.getCol(clustersColumnName);
  const origClustColData = origClustCol.getRawData();
  const origClustColCat = origClustCol.categories;
  const origClustMasks: BitArray[] = Array.from({length: origClustColCat.length},
    () => new BitArray(rowCount, false));
  for (let rowIdx = 0; rowIdx < rowCount; ++rowIdx)
    origClustMasks[origClustColData[rowIdx]].setTrue(rowIdx);

  const customClustMasks = customClusters.map(
    (v) => BitArray.fromUint32Array(rowCount, v.getRawData() as Uint32Array));
  const customClustColNamesList = customClusters.map((v) => v.name);

  const activityColData: type.RawData = df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED).getRawData();

  const origClustStats: ClusterStats = {};
  const customClustStats: ClusterStats = {};

  for (const clustType of Object.values(CLUSTER_TYPE)) {
    const masks = clustType === CLUSTER_TYPE.ORIGINAL ? origClustMasks : customClustMasks;
    const clustNames = clustType === CLUSTER_TYPE.ORIGINAL ? origClustColCat : customClustColNamesList;
    const resultStats = clustType === CLUSTER_TYPE.ORIGINAL ? origClustStats : customClustStats;
    for (let maskIdx = 0; maskIdx < masks.length; ++maskIdx) {
      const mask = masks[maskIdx];
      const stats = mask.allTrue || mask.allFalse ? {count: mask.length, meanDifference: 0, ratio: 1.0, pValue: null, mask: mask} :
        getStats(activityColData, mask);
      resultStats[clustNames[maskIdx]] = stats;
    }
  }

  const resultStats = {} as ClusterTypeStats;
  resultStats[CLUSTER_TYPE.ORIGINAL] = origClustStats;
  resultStats[CLUSTER_TYPE.CUSTOM] = customClustStats;
  return resultStats;
}
