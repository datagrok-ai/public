/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as C from './constants';
import * as type from './types';
import {ParallelMutationCliffs} from './parallel-mutation-cliffs';
import {CLUSTER_TYPE} from '../viewers/logo-summary';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {
  bitSetToBitArray,
  ClusterStats,
  ClusterTypeStats,
  getStats,
  MonomerPositionStats,
  PositionStats,
  StatsItem,
  SummaryStats,
} from './statistics';

export type MutationCliffsOptions = {
  maxMutations?: number,
  minActivityDelta?: number,
  filter?: Uint32Array,
  singlePosition?: {position: number}
};

/**
 * Finds mutation cliffs in the set of sequences.
 * @param activityArray - Activity column raw data.
 * @param monomerInfoArray - Split sequence raw columns.
 * @param options - Options for the mutation cliffs algorithm.
 * @return - Mutation cliffs map.
 */
export async function findMutations(activityArray: type.RawData, monomerInfoArray: type.RawColumn[],
  options: MutationCliffsOptions = {}): Promise<type.MutationCliffs> {
  const nCols = monomerInfoArray.length;
  if (nCols === 0)
    throw new Error(`PepAlgorithmError: Couldn't find any column of semType '${C.SEM_TYPES.MONOMER}'`);


  options.minActivityDelta ??= 0;
  options.maxMutations ??= 1;
  const mutationCliffsService = new ParallelMutationCliffs();
  const substitutionsInfo = await mutationCliffsService.calc(activityArray, monomerInfoArray, options);
  mutationCliffsService.terminate();
  return substitutionsInfo;
}

/**
 * Calculates statistics for mutation cliffs, used for mutation cliffst table (coloring, tooltips, distribution...)
 * @param cliffs - mutation cliffs data
 * @param activityArray - array of activities
 *
 */
export function calculateCliffsStatistics(
  cliffs: type.MutationCliffs, activityArray: type.RawData,
): type.MutationCliffStats {
  const res: type.MutationCliffStats['stats'] = new Map();
  let minDiff = 999999; let maxDiff = -999999; let minCount = 2; let maxCount = 2;
  for (const monomer of cliffs.keys()) {
    const monomerStatsMap: Map<string, StatsItem> = new Map();
    res.set(monomer, monomerStatsMap);
    // monomer substitutions map from mutations cliffs
    const monomerSubMap = cliffs.get(monomer)!;
    for (const position of monomerSubMap.keys()) {
      const subMap = monomerSubMap.get(position)!;
      if (subMap.size === 0)
        continue;
      // create two masks, one filtering the activities to only mutation cliffs (with given monomer at given position and its substitutions)
      // another one corresponding to only the given monomer at given position within the mutation cliffs
      const filterMask = new BitArray(activityArray.length, false);
      const maskMCWithMonomerAtPosition = new BitArray(activityArray.length, false);
      for (const index of subMap.keys()) {
        // set the filter mask to true for all sequences within the mutation cliff pairs
        filterMask.setFast(index, true);
        const toIndexes = subMap.get(index)!;
        toIndexes.forEach((i) => filterMask.setFast(i, true));
        // set the mask for sequences with the given monomer at the given position within the mutation cliffs
        maskMCWithMonomerAtPosition.setFast(index, true);
      }
      const stats = getStats(activityArray, maskMCWithMonomerAtPosition, filterMask);
      stats.mask = filterMask; // store the filter mask for later use in the viewer
      minDiff = Math.min(minDiff, stats.meanDifference);
      maxDiff = Math.max(maxDiff, stats.meanDifference);
      minCount = Math.min(minCount, stats.count);
      maxCount = Math.max(maxCount, stats.count);
      // here, stats will show the following
      // count - number of sequences with the given monomer at the given position within the mutation cliffs
      // mask.trueCount - number of unique sequences within the mutation cliffs (with given monomer at given position and its substitutions)
      // meanDifference - difference between mean activity of sequences with the given monomer at the given position within the mutation cliffs
      // and mean activity of other sequences within the mutation cliffs (with given monomer at given position substitutions)
      monomerStatsMap.set(position, stats);
    }
  }
  return {stats: res, minDiff, maxDiff, minCount, maxCount};
}


/**
 * Calculates statistics for each monomer position.
 * @param activityCol - Activity column.
 * @param filter - Dataframe filter to consider.
 * @param positionColumns - Position columns containing monomers.
 * @param [options] - Options for the algorithm.
 * @param [options.isFiltered] - Whether the dataframe is filtered.
 * @param [options.columns] - Columns to consider when calculating statistics.
 * @param [options.target] - Target column and category to consider.
 * @param [options.aggValue] - Column and aggregation type to consider instead of count.
 * @return - Statistics for each monomer position.
 */
export function calculateMonomerPositionStatistics(activityCol: DG.Column<number>, filter: DG.BitSet,
  positionColumns: DG.Column<string>[], options: {
    isFiltered?: boolean,
    columns?: string[],
    aggValue?: {
      col: DG.Column,
      type: DG.AGG
    }
  } = {}): MonomerPositionStats {
  options.isFiltered ??= false;
  const monomerPositionObject = {general: {}} as MonomerPositionStats & { general: SummaryStats };
  const activityColData: Float64Array = activityCol.getRawData() as Float64Array;
  const sourceDfLen = activityCol.length;
  options.columns ??= positionColumns.map((col) => col.name);
  // if (options.isFiltered) {
  //   sourceDfLen = filter.trueCount;
  //   const tempActivityData = new Float64Array(sourceDfLen);
  //   const selectedIndexes = filter.getSelectedIndexes();
  //   for (let i = 0; i < sourceDfLen; ++i)
  //     tempActivityData[i] = activityColData[selectedIndexes[i]];


  //   activityColData = tempActivityData;
  //   positionColumns = DG.DataFrame.fromColumns(positionColumns).clone(filter, options.columns).columns.toList();
  //   if (options.aggValue)
  //     options.aggValue.col = options.aggValue.col.clone(filter);
  // }
  const filterBitArray = options.isFiltered ? bitSetToBitArray(filter) : undefined;

  for (const posCol of positionColumns) {
    if (!options.columns.includes(posCol.name))
      continue;


    const posColData = posCol.getRawData();
    const posColCateogries = posCol.categories;
    const currentPositionObject = {general: {}} as PositionStats & { general: SummaryStats };

    for (let categoryIndex = 0; categoryIndex < posColCateogries.length; ++categoryIndex) {
      const monomer = posColCateogries[categoryIndex];
      if (monomer === '')
        continue;


      const boolArray: boolean[] = new Array(sourceDfLen).fill(false);
      for (let i = 0; i < sourceDfLen; ++i) {
        if (posColData[i] === categoryIndex && (!options.isFiltered || filter.get(i)))
          boolArray[i] = true;
      }
      const bitArray = BitArray.fromValues(boolArray);
      if (bitArray.allFalse)
        continue;
      const stats = getStats(activityColData, bitArray, filterBitArray, options.aggValue);
      currentPositionObject[monomer] = stats;
      getSummaryStats(currentPositionObject.general, stats);
    }
    monomerPositionObject[posCol.name] = currentPositionObject;
    getSummaryStats(monomerPositionObject.general, null, currentPositionObject.general);
  }
  return monomerPositionObject;
}

/**
 * Calculates summary statistics for the monomer position statistics such as maximum and minimum values for each
 * statistic in general and on each position.
 * @param genObj - Object to store the summary statistics to.
 * @param stats - Statistics for a single monomer position.
 * @param summaryStats - Summary statistics for all monomer positions.
 */
export function getSummaryStats(genObj: SummaryStats, stats: StatsItem | null = null,
  summaryStats: SummaryStats | null = null): void {
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

/**
 * Calculates statistics for each cluster type.
 * @param df - Dataframe containing the clusters column.
 * @param clustersColumnName - Name of the original clusters column.
 * @param customClusters - Array of custom clusters columns names.
 * @param activityCol - Activity column.
 * @return - Statistics for each cluster type.
 */
export function calculateClusterStatistics(df: DG.DataFrame, clustersColumnName: string,
  customClusters: DG.Column<boolean>[], activityCol: DG.Column<number>): ClusterTypeStats {
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

  const activityColData = activityCol.getRawData() as Float64Array;

  const origClustStats: ClusterStats = {};
  const customClustStats: ClusterStats = {};

  const filterMask = df.filter.anyFalse ? bitSetToBitArray(df.filter) : undefined;
  for (const clustType of Object.values(CLUSTER_TYPE)) {
    const masks = clustType === CLUSTER_TYPE.ORIGINAL ? origClustMasks : customClustMasks;
    const clustNames = clustType === CLUSTER_TYPE.ORIGINAL ? origClustColCat : customClustColNamesList;
    const resultStats = clustType === CLUSTER_TYPE.ORIGINAL ? origClustStats : customClustStats;
    for (let maskIdx = 0; maskIdx < masks.length; ++maskIdx) {
      const mask = masks[maskIdx];
      resultStats[clustNames[maskIdx]] = mask.allTrue || mask.allFalse ?
        {count: mask.length, meanDifference: 0, ratio: 1.0, pValue: null, mask: mask, mean: activityCol.stats.avg} :
        getStats(activityColData, mask, filterMask);
    }
  }

  const resultStats = {} as ClusterTypeStats;
  resultStats[CLUSTER_TYPE.ORIGINAL] = origClustStats;
  resultStats[CLUSTER_TYPE.CUSTOM] = customClustStats;
  return resultStats;
}
