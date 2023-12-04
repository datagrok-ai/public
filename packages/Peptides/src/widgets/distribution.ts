import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';

import $ from 'cash-dom';

import * as C from '../utils/constants';
import {AggregationColumns, getAggregatedColumnValues, getStats, StatsItem} from '../utils/statistics';
import {getDistributionPanel, getDistributionTable} from '../utils/misc';
import {SARViewer} from '../viewers/sar-viewer';
import {LogoSummaryTable} from '../viewers/logo-summary';

export type DistributionTableOptions = { peptideSelection: DG.BitSet, columns: AggregationColumns };
export type PeptideViewer = SARViewer | LogoSummaryTable;

export function getDistributionWidget(table: DG.DataFrame, options: DistributionTableOptions): HTMLDivElement {
  if (!table.selection.anyTrue)
    return ui.divText('No distribution');

  const activityCol = table.getCol(C.COLUMNS_NAMES.ACTIVITY);
  const rowCount = activityCol.length;

  const distributionHost = ui.div([], 'd4-flex-wrap');

  const updateDistributionHost = (): void => {
    const res: HTMLDivElement[] = [];
    if (!table.selection.anyTrue)
      res.push(ui.divText('No distribution'));
    else {
      const hist = getActivityDistribution(getDistributionTable(activityCol, table.selection,
        options.peptideSelection));
      const bitArray = BitArray.fromString(table.selection.toBinaryString());
      const mask = DG.BitSet.create(rowCount,
        bitArray.allFalse || bitArray.allTrue ? (_): boolean => true : (i): boolean => bitArray.getBit(i));
      const aggregatedColMap = getAggregatedColumnValues(table, options.columns, {filterDf: true, mask});
      const stats = bitArray.allFalse || bitArray.allTrue ?
        {count: rowCount, pValue: null, meanDifference: 0, ratio: 1, mask: bitArray, mean: activityCol.stats.avg} :
        getStats(activityCol.getRawData(), bitArray);
      const tableMap = getStatsTableMap(stats);
      const resultMap: { [key: string]: any } = {...tableMap, ...aggregatedColMap};
      const distributionRoot = getDistributionPanel(hist, resultMap);
      $(distributionRoot).addClass('d4-flex-col');

      res.push(distributionRoot);
    }

    $(distributionHost).empty().append(res);
  };

  updateDistributionHost();
  return ui.divV([distributionHost]);
}

export function getActivityDistribution(table: DG.DataFrame, isTooltip: boolean = false,
): DG.Viewer<DG.IHistogramLookSettings> {
  const hist = table.plot.histogram({
    filteringEnabled: false,
    valueColumnName: C.COLUMNS_NAMES.ACTIVITY,
    splitColumnName: C.COLUMNS_NAMES.SPLIT_COL,
    legendVisibility: 'Never',
    showXAxis: true,
    showColumnSelector: false,
    showRangeSlider: false,
    showBinSelector: false,
    backColor: isTooltip ? '#fdffe5' : '#fffff',
  }) as DG.Viewer<DG.IHistogramLookSettings>;
  hist.root.style.width = 'auto';
  return hist;
}

export function getStatsTableMap(stats: StatsItem, options: { fractionDigits?: number } = {}): StringDictionary {
  options.fractionDigits ??= 3;
  const tableMap: StringDictionary = {
    'Count': `${stats.count} (${stats.ratio.toFixed(options.fractionDigits)}%)`,
    'Mean difference': stats.meanDifference.toFixed(options.fractionDigits),
    'Mean activity': stats.mean.toFixed(options.fractionDigits),
  };
  if (stats.pValue !== null)
    tableMap['p-value'] = stats.pValue < 0.01 ? '<0.01' : stats.pValue.toFixed(options.fractionDigits);
  return tableMap;
}
