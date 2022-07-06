import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as C from '../utils/constants';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {getStats, Stats} from '../utils/filtering-statistics';
import * as type from '../utils/types';

const allLabel = 'All';

export function getDistributionWidget(table: DG.DataFrame): DG.Widget {
  const splitCol = table.col(C.COLUMNS_NAMES.SPLIT_COL);
  const activityScaledCol = table.columns.bySemType(C.SEM_TYPES.ACTIVITY_SCALED)!;
  const selectionObject: type.SelectionObject = JSON.parse(table.tags[C.TAGS.SELECTION]);
  if (!splitCol || !selectionObject)
    return new DG.Widget(ui.divText('No distribution'));

  const positions = Object.keys(selectionObject);
  let aarStr = allLabel;
  let otherStr = '';

  if (positions.length) {
    aarStr = '';

    for (const position of positions) {
      aarStr += `${position}: {`;

      for (const aar of selectionObject[position])
        aarStr += `${aar}, `;

      aarStr = aarStr.slice(0, aarStr.length - 2);
      aarStr += '}; ';
    }
    aarStr = aarStr.slice(0, aarStr.length - 2);
    otherStr = 'Other';
  }

  const distributionTable = DG.DataFrame.fromColumns([activityScaledCol, splitCol]);
  const stats = getStats(activityScaledCol.toList(), table.selection);
  return new DG.Widget(getDistributionAndStats(distributionTable, stats, aarStr, otherStr));
}

export function getDistributionAndStats(
    table: DG.DataFrame, stats: Stats, thisLabel: string, otherLabel: string = '', isTooltip: boolean = false,
  ): HTMLDivElement {
  const labels = ui.divV([
    ui.label(thisLabel, {style: {color: DG.Color.toHtml(otherLabel == '' ? DG.Color.blue : DG.Color.orange)}}),
    ui.label(otherLabel, {style: {color: DG.Color.toHtml(DG.Color.blue)}})]);

  const histRoot = table.plot.histogram({
    filteringEnabled: false,
    valueColumnName: table.columns.bySemType(C.SEM_TYPES.ACTIVITY_SCALED)?.name,
    splitColumnName: C.COLUMNS_NAMES.SPLIT_COL,
    legendVisibility: 'Never',
    showXAxis: true,
    showColumnSelector: false,
    showRangeSlider: false,
    showBinSelector: !isTooltip,
    backColor: isTooltip ? '#fdffe5' : '#fffff',
  }).root;
  histRoot.style.width = 'auto';

  const tableMap: StringDictionary = {
    'Statistics:': '',
    'Count': stats.count.toString(),
    'Ratio': stats.ratio.toFixed(2),
    'p-value': stats.pValue < 0.01 ? '<0.01' : stats.pValue.toFixed(2),
    'Mean difference': stats.meanDifference.toFixed(2),
  };

  const result = ui.divV([labels, histRoot, ui.tableFromMap(tableMap)]);
  result.style.minWidth = '200px';
  if (isTooltip)
    histRoot.style.maxHeight = '150px';
  return result;
}
