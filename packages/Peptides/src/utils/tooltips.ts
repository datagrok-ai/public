import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import * as type from './types';
import * as C from '../utils/constants';

import {getActivityDistribution, getStatsTableMap} from '../widgets/distribution';
import {getDistributionPanel, getDistributionTable} from './misc';
import {getMonomerWorksInstance} from '../package';
import {AggregationColumns, getAggregatedColumnValues, MonomerPositionStats} from './statistics';

export type TooltipOptions = {
  fromViewer?: boolean, isMutationCliffs?: boolean, x: number, y: number,
  monomerPosition: type.SelectionItem, mpStats: MonomerPositionStats
};

export function showMonomerTooltip(monomer: string, x: number, y: number): boolean {
  const tooltipElements: HTMLDivElement[] = [];
  const monomerName = monomer.toLowerCase();

  const mw = getMonomerWorksInstance();
  const mol = mw?.getCappedRotatedMonomer('PEPTIDE', monomer);

  if (mol) {
    tooltipElements.push(ui.div(monomerName));
    const options = {autoCrop: true, autoCropMargin: 0, suppressChiralText: true};
    tooltipElements.push(grok.chem.svgMol(mol, undefined, undefined, options));
  } else if (monomer !== '')
    tooltipElements.push(ui.div(monomer));
  else
    return true;


  ui.tooltip.show(ui.divV(tooltipElements), x, y);

  return true;
}

export function showTooltip(df: DG.DataFrame, activityCol: DG.Column<number>, columns: AggregationColumns,
  options: TooltipOptions): boolean {
  options.fromViewer ??= false;
  options.isMutationCliffs ??= false;
  if (options.monomerPosition.positionOrClusterType === C.COLUMNS_NAMES.MONOMER)
    showMonomerTooltip(options.monomerPosition.monomerOrCluster, options.x, options.y);
  else
    showTooltipAt(df, activityCol, columns, options);
  return true;
}

//TODO: move out to viewer code
export function showTooltipAt(df: DG.DataFrame, activityCol: DG.Column<number>, columns: AggregationColumns,
  options: TooltipOptions): HTMLDivElement | null {
  options.fromViewer ??= false;
  options.isMutationCliffs ??= false;
  const stats = options
    .mpStats[options.monomerPosition.positionOrClusterType]![options.monomerPosition.monomerOrCluster];
  if (!stats?.count)
    return null;

  const mask = DG.BitSet.fromBytes(stats.mask.buffer.buffer, activityCol.length);
  const hist = getActivityDistribution(getDistributionTable(activityCol, mask), true);

  const tableMap = getStatsTableMap(stats);
  if (options.fromViewer) {
    tableMap['Mean difference'] = `${tableMap['Mean difference']}${options.isMutationCliffs ? ' (size)' : ''}`;
    if (tableMap['p-value'])
      tableMap['p-value'] = `${tableMap['p-value']}${options.isMutationCliffs ? ' (color)' : ''}`;
  }
  const aggregatedColMap = getAggregatedColumnValues(df, Object.entries(columns), {mask: mask});
  const resultMap = {...tableMap, ...aggregatedColMap};

  const distroStatsElem = getDistributionPanel(hist, resultMap);

  ui.tooltip.show(distroStatsElem, options.x, options.y);

  return distroStatsElem;
}
