import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import * as type from './types';
import * as C from '../utils/constants';

import {getActivityDistribution, getStatsTableMap} from '../widgets/distribution';
import {getDistributionPanel, getDistributionTable} from './misc';
import {getAggregatedColumnValues, MonomerPositionStats} from './statistics';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';

export type TooltipOptions = {
  fromViewer?: boolean, isMutationCliffs?: boolean, x: number, y: number, monomerPosition: type.SelectionItem,
  mpStats: MonomerPositionStats, aggrColValues?: StringDictionary,
  isMostPotentResidues?: boolean, cliffStats?: type.MutationCliffStats['stats'],
  postfixes?: StringDictionary, additionalStats?: StringDictionary
};

/**
 * Shows tooltip at the given coordinates.
 * @param df - Dataframe to show tooltip at.
 * @param activityCol - Activity column.
 * @param columns - Aggregation columns.
 * @param options - Tooltip options.
 * @return - Flag if the tooltip is shown.
 */
export function showTooltip(df: DG.DataFrame, activityCol: DG.Column<number>, columns: [string, DG.AggregationType][],
  options: TooltipOptions): boolean {
  options.fromViewer ??= false;
  options.isMutationCliffs ??= false;
  options.isMostPotentResidues ??= false;
  if (options.monomerPosition.positionOrClusterType !== C.COLUMNS_NAMES.MONOMER)
    showTooltipAt(df, activityCol, columns, options);


  return true;
}

//TODO: move out to viewer code
/**
 * Shows dataframe tooltip at the given coordinates.
 * @param df - Dataframe to show tooltip at.
 * @param activityCol - Activity column.
 * @param columns - Aggregation columns.
 * @param options - Tooltip options.
 * @return - Tooltip element.
 */
export function showTooltipAt(df: DG.DataFrame, activityCol: DG.Column<number>, columns: [string, DG.AggregationType][],
  options: TooltipOptions): HTMLDivElement | null {
  options.fromViewer ??= false;
  options.isMutationCliffs ??= false;
  options.isMostPotentResidues ??= false;
  options.additionalStats ??= {};
  if (!options.cliffStats || !options.isMutationCliffs) {
    const stats = options
      .mpStats[options.monomerPosition.positionOrClusterType]![options.monomerPosition.monomerOrCluster];
    if (!stats?.count)
      return null;


    const mask = DG.BitSet.fromBytes(stats.mask.buffer.buffer, activityCol.length);
    const hist = getActivityDistribution(getDistributionTable(activityCol, mask), true);
    const tableMap = getStatsTableMap(stats);
    if (options.fromViewer) {
      tableMap['Mean difference'] = `${tableMap['Mean difference']}${options.isMostPotentResidues ? ' (size)' : ''}`;
      if (tableMap['p-value'])
        tableMap['p-value'] = `${tableMap['p-value']}${options.isMostPotentResidues ? ' (color)' : ''}`;
    }
    const aggregatedColMap = options.aggrColValues ?? getAggregatedColumnValues(df, columns, {mask: mask});
    const resultMap = {...options.additionalStats, ...tableMap, ...aggregatedColMap};
    for (const [key, value] of Object.entries(options.postfixes ?? {})) {
      if (resultMap[key])
        resultMap[key] = `${resultMap[key]}${value}`;
    }
    const distroStatsElem = getDistributionPanel(hist, resultMap);
    ui.tooltip.show(distroStatsElem, options.x, options.y);
    return distroStatsElem;
  } else {
    const stats = options.cliffStats?.get(options.monomerPosition.monomerOrCluster)
      ?.get(options.monomerPosition.positionOrClusterType)
      ;
    if (!stats)
      return null;
    const mask = DG.BitSet.fromBytes(stats.mask.buffer.buffer, activityCol.length);
    const hist = getActivityDistribution(getDistributionTable(activityCol, mask), true);
    const tableMap = getStatsTableMap(stats, {countName: 'Unique count'});
    if (options.fromViewer) {
      tableMap['Mean difference'] = `${tableMap['Mean difference']}${' (Color)'}`;
      if (tableMap['Unique count'])
        tableMap['Unique count'] = `${tableMap['Unique count']}${' (Size)'}`;
    }
    const aggregatedColMap = options.aggrColValues ?? getAggregatedColumnValues(df, columns, {mask: mask});
    const resultMap = {...options.additionalStats, ...tableMap, ...aggregatedColMap};
    const distroStatsElem = getDistributionPanel(hist, resultMap);
    ui.tooltip.show(distroStatsElem, options.x, options.y);
    return distroStatsElem;
  }
}
