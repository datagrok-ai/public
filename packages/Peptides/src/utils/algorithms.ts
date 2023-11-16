import * as DG from 'datagrok-api/dg';
import * as C from './constants';
import * as type from './types';
import {ParallelMutationCliffs} from './parallel-mutation-cliffs';
import {CLUSTER_TYPE} from '../model';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {ClusterStats, ClusterTypeStats, MonomerPositionStats, PositionStats, Stats, SummaryStats, getStats} from './statistics';

export type TargetOptions = {targetCol?: type.RawColumn | null, currentTarget?: string | null};

export async function findMutations(activityArray: type.RawData, monomerInfoArray: type.RawColumn[],
  settings: type.PeptidesSettings = {}, targetOptions: TargetOptions = {}): Promise<type.MutationCliffs> {
  const nCols = monomerInfoArray.length;
  if (nCols === 0)
    throw new Error(`PepAlgorithmError: Couldn't find any column of semType '${C.SEM_TYPES.MONOMER}'`);

  settings.minActivityDelta ??= 0;
  settings.maxMutations ??= 1;
  const substitutionsInfo = await new ParallelMutationCliffs().calc(activityArray, monomerInfoArray, settings, targetOptions);
  return substitutionsInfo;
}

export function calculateMonomerPositionStatistics(df: DG.DataFrame, positionColumns: DG.Column<string>[],
  options: {isFiltered?: boolean, columns?: string[]} = {}): MonomerPositionStats {
  options.isFiltered ??= false;
  const monomerPositionObject = {general: {}} as MonomerPositionStats & {general: SummaryStats};
  const activityCol = df.getCol(C.COLUMNS_NAMES.ACTIVITY);
  let activityColData: Float64Array = activityCol.getRawData() as Float64Array;
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
        {count: sourceDfLen, meanDifference: 0, ratio: 1.0, pValue: null, mask: bitArray, mean: activityCol.stats.avg} :
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

  const activityCol = df.getCol(C.COLUMNS_NAMES.ACTIVITY);
  const activityColData = activityCol.getRawData() as Float64Array;

  const origClustStats: ClusterStats = {};
  const customClustStats: ClusterStats = {};

  for (const clustType of Object.values(CLUSTER_TYPE)) {
    const masks = clustType === CLUSTER_TYPE.ORIGINAL ? origClustMasks : customClustMasks;
    const clustNames = clustType === CLUSTER_TYPE.ORIGINAL ? origClustColCat : customClustColNamesList;
    const resultStats = clustType === CLUSTER_TYPE.ORIGINAL ? origClustStats : customClustStats;
    for (let maskIdx = 0; maskIdx < masks.length; ++maskIdx) {
      const mask = masks[maskIdx];
      const stats = mask.allTrue || mask.allFalse ?
        {count: mask.length, meanDifference: 0, ratio: 1.0, pValue: null, mask: mask, mean: activityCol.stats.avg} :
        getStats(activityColData, mask);
      resultStats[clustNames[maskIdx]] = stats;
    }
  }

  const resultStats = {} as ClusterTypeStats;
  resultStats[CLUSTER_TYPE.ORIGINAL] = origClustStats;
  resultStats[CLUSTER_TYPE.CUSTOM] = customClustStats;
  return resultStats;
}
