// Which viewer charts a comparison result, with its options and minimum height; no rendering.

import * as DG from 'datagrok-api/dg';
import {RUN_COLUMN} from './types';
import {Comparison} from './comparison-builders';

export type ScalarChartType = 'radar' | 'pcplot';

export type ChartSpec =
  | {kind: 'hint', text: string}
  | {
    kind: 'viewer',
    type: string,
    options: Record<string, any>,
    minHeight: number,
    // radar vs PC plot switch above a multi-scalar chart
    scalarChartSwitch: boolean,
    // shared vs independent y scale switch above a multi-value line chart
    sharedScaleSwitch: boolean,
  };

const hint = (text: string): ChartSpec => ({kind: 'hint', text});

const viewer = (
  type: string, options: Record<string, any>, minHeight: number, scalarChartSwitch = false, sharedScaleSwitch = false,
): ChartSpec => ({kind: 'viewer', type, options, minHeight, scalarChartSwitch, sharedScaleSwitch});

export function chartSpec(
  comparison: Comparison | null,
  hasSelection: boolean,
  scalarChartType: ScalarChartType,
  radarAvailable: boolean,
  sharedScale = true,
): ChartSpec {
  if (!hasSelection)
    return hint('Select what to compare above');
  if (!comparison || comparison.chartDf.rowCount === 0)
    return hint('Nothing to show: no data points across the selected runs');

  if (comparison.kind === 'multi-scalar') {
    // a 2-axis radar degenerates to a line; the PC plot (a slopegraph) charts it until a 3rd value
    const radarEligible = radarAvailable && comparison.valueColumnNames.length >= 3;
    return radarEligible && scalarChartType === 'radar' ?
      viewer('Radar', {
        valuesColumnNames: comparison.valueColumnNames,
        colorColumnName: RUN_COLUMN,
        showCurrentRow: false,
      }, 300, true) :
      viewer(DG.VIEWER.PC_PLOT, {
        columnNames: comparison.valueColumnNames,
        colorColumnName: RUN_COLUMN,
      }, 300, radarEligible);
  }

  if (comparison.kind === 'scalar') {
    return viewer(DG.VIEWER.BAR_CHART, {
      valueColumnName: comparison.valueColumnName,
      valueAggrType: 'avg',
      splitColumnName: RUN_COLUMN,
    }, 250);
  }

  const {valueColumnNames} = comparison;
  if (comparison.isKeyIndex) {
    // multi-value suggestions need a line-chartable index, but the index may switch to a
    // string column after selection, and datetime mixed with numeric charts only as elapsed time
    if (valueColumnNames.length > 1) {
      return hint(comparison.mixedTimeIndex ?
        'Datetime and numeric indexes chart together only in the relative timeseries mode — set it ' +
          'on every table, or select a single value' :
        'Multiple values cannot chart over a categorical index — pick a numeric or datetime ' +
          'index, or select a single value');
    }
    return viewer(DG.VIEWER.BAR_CHART, {
      valueColumnName: comparison.valueColumnName,
      valueAggrType: 'avg',
      splitColumnName: comparison.indexColumnName,
      stackColumnName: RUN_COLUMN,
    }, 250);
  }

  if (comparison.isScatter) {
    // color: value name (multi-value melt) > split > Run; runs move to marker shapes when color is taken
    const {melted, splitColumnName} = comparison;
    return viewer(DG.VIEWER.SCATTER_PLOT, melted ? {
      xColumnName: comparison.indexColumnName,
      yColumnName: melted.valueColumnName,
      colorColumnName: melted.seriesColumnName,
      markersColumnName: RUN_COLUMN,
    } : {
      xColumnName: comparison.indexColumnName,
      yColumnName: comparison.valueColumnName,
      colorColumnName: splitColumnName ?? RUN_COLUMN,
      ...splitColumnName ? {markersColumnName: RUN_COLUMN} : {},
    }, 300);
  }

  const shareable = valueColumnNames.length > 1 && comparison.sameUnits;
  const shared = shareable && sharedScale;
  return viewer(DG.VIEWER.LINE_CHART, {
    xColumnName: comparison.indexColumnName,
    yColumnNames: valueColumnNames,
    splitColumnNames: [RUN_COLUMN, comparison.stepColumnName, comparison.splitColumnName]
      .filter((name): name is string => !!name),
    multiAxis: shared,
    yGlobalScale: shared,
  }, shared ? 250 : 250 * Math.max(1, valueColumnNames.length), false, shareable);
}
