/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {
  getSeriesFitFunction,
  getChartBounds,
} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {LogOptions} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {FitChartCellRenderer} from './fit-renderer';
import {
  getOrCreateParsedChartData, getColumnChartOptions, getDataFrameChartOptions, mergeProperties,
  substituteZeroes, CHART_OPTIONS, SERIES_OPTIONS,
} from './fit-chart-data';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {ColorType, SeriesColorType, getSeriesColor, isColorValid} from './render-utils';
import {calculateSeriesFit, getChartDataAggrStats, aggregatedStatisticsProperties} from './fit-statistics';
import {fitChartDataProperties} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {fitSeriesProperties} from '@datagrok-libraries/statistics/src/fit/fit-engine';
import {
  changeCurvesOptions, chartPropertiesFor, normalizeStatisticNames, seriesPropertiesFor, ManipulationLevel,
} from './fit-options';

type CurveStatisticParams = {propName: string, seriesNumber: number} | {propName: string, aggrType: string};

/** Runs a curve statistic so its `join(table)` output makes the platform add it as a calculated column. */
export async function addStatisticColumn(gridCell: DG.GridCell, funcName: string,
  params: CurveStatisticParams): Promise<void> {
  const func = DG.Func.find({package: 'Curves', name: funcName})[0];
  if (!func) {
    grok.shell.error(`Curves: ${funcName} is not registered`);
    return;
  }
  const df = gridCell.cell.dataFrame;
  const curveColumnName = gridCell.cell.column.name;
  const before = new Set(df.columns.names());
  try {
    await func.prepare({table: df, curveColumn: gridCell.cell.column, ...params})
      .call(false, undefined, {processed: false});
  } catch (e) {
    // the caller is an icon click handler, so an unhandled rejection just looks like a dead button
    grok.shell.error(`Could not add the statistic column: ${e}`);
    return;
  }
  placeAfterCurveColumn(gridCell, df.columns.names().filter((name) => !before.has(name)), curveColumnName);
}

/** `join(table)` appends, so the position is restored afterwards - in the grid too, which takes a
 * column's position from the dataframe only when it first builds a GridColumn for it. */
function placeAfterCurveColumn(gridCell: DG.GridCell, added: string[], curveColumnName: string): void {
  if (added.length === 0)
    return;
  const df = gridCell.cell.dataFrame;
  df.columns.setOrder(insertedAfter(df.columns.names(), added, curveColumnName));

  const grid = gridCell.grid;
  if (!grid)
    return;
  const gridNames: string[] = [];
  for (let i = 0; i < grid.columns.length; i++) {
    const name = grid.columns.byIndex(i)?.name;
    if (name)
      gridNames.push(name);
  }
  grid.columns.setOrder(insertedAfter(gridNames, added, curveColumnName));
}

function insertedAfter(names: string[], added: string[], anchor: string): string[] {
  const rest = names.filter((name) => !added.includes(name));
  const at = rest.indexOf(anchor) + 1;
  return [...rest.slice(0, at), ...added, ...rest.slice(at)];
}

/** A multi-choice is capped at 155px and always draws its scrollbar, which cut the statistics off and
 * left a second scrollbar inside a panel that already scrolls. These lists are short enough to read. */
function showWholeList(input: DG.InputBase): DG.InputBase {
  const checks = input.root.querySelector('.ui-input-multi-choice-checks') as HTMLElement | null;
  if (checks) {
    checks.style.maxHeight = 'none';
    checks.style.overflowY = 'visible';
  }
  return input;
}

export class FitGridCellHandler extends DG.ObjectHandler {
  get type(): string {
    return 'GridCell';
  }

  isApplicable(x: any): boolean {
    return x instanceof DG.GridCell && x.cellType === FitConstants.FIT_CELL_TYPE;
  }

  // TODO: add aspect ratio for the cell
  // TODO: add legend
  // TODO: add the table for the values on the cell or don't render it at all
  // TODO: fix the curves demo app

  renderProperties(gridCell: DG.GridCell, _context: any = null): HTMLElement {
    const acc = ui.accordion('Curves property panel');
    // TODO: make just the base ui.input.choice after nullable option is added
    const switchProperty = DG.Property.js('level', DG.TYPE.STRING, {description: 'Controls the level at which properties will be switched',
      defaultValue: 'Column', choices: ['Dataframe', 'Column', 'Cell'], nullable: false});
    const switchLevelInput = ui.input.forProperty(switchProperty);

    // temporarily because input doesn't show the tooltip
    ui.tooltip.bind(switchLevelInput.captionLabel, 'Controls the level at which properties will be switched');

    const chartData = getOrCreateParsedChartData(gridCell.cell);
    const columnChartOptions = getColumnChartOptions(gridCell.cell.column);
    const dfChartOptions = getDataFrameChartOptions(gridCell.cell.dataFrame);

    const seriesOptionsRefresh = {onValueChanged: (v: any, inputBase: DG.InputBase) =>
      changeCurvesOptions(gridCell, inputBase, SERIES_OPTIONS, switchLevelInput.value as ManipulationLevel)};
    const chartOptionsRefresh = {onValueChanged: (v: any, inputBase: DG.InputBase) =>
      changeCurvesOptions(gridCell, inputBase, CHART_OPTIONS, switchLevelInput.value as ManipulationLevel)};

    const setValidColors = (colorFieldName: SeriesColorType) => {
      if (dfChartOptions.seriesOptions && !isColorValid(dfChartOptions.seriesOptions[colorFieldName]) &&
        columnChartOptions.seriesOptions && !isColorValid(columnChartOptions.seriesOptions[colorFieldName])) {
        if (chartData.seriesOptions) {
          if (!isColorValid(chartData.seriesOptions[colorFieldName])) {
            chartData.seriesOptions[colorFieldName] = DG.Color.toHtml(colorFieldName === 'outlierColor' ?
              DG.Color.red : DG.Color.getCategoricalColor(0));
          }
        } else {
          if (!isColorValid(chartData.series ? chartData.series[0][colorFieldName] : '')) {
chartData.series![0][colorFieldName] = DG.Color.toHtml(colorFieldName === 'outlierColor' ?
  DG.Color.red : DG.Color.getCategoricalColor(0));
          }
        }
      }
    };
    const tableCell = gridCell.cell;
    setValidColors('pointColor');
    setValidColors('fitLineColor');
    setValidColors('outlierColor');

    mergeProperties(fitSeriesProperties, columnChartOptions.seriesOptions, chartData.seriesOptions ? chartData.seriesOptions :
      chartData.series ? chartData.series[0] ?? {} : {});
    mergeProperties(fitSeriesProperties, dfChartOptions.seriesOptions, chartData.seriesOptions ? chartData.seriesOptions :
      chartData.series ? chartData.series[0] ?? {} : {});
    mergeProperties(fitChartDataProperties, columnChartOptions.chartOptions,
      chartData.chartOptions ? chartData.chartOptions : {});
    mergeProperties(fitChartDataProperties, dfChartOptions.chartOptions,
      chartData.chartOptions ? chartData.chartOptions : {});

    if (chartData.chartOptions?.allowXZeroes && chartData.chartOptions?.logX &&
      chartData.series?.some((series) => series.points.some((p) => p.x === 0)))
      substituteZeroes(chartData);

    const form = ui.form([switchLevelInput]);
    const seriesSource = chartData.seriesOptions ? chartData.seriesOptions : chartData.series![0];
    // a custom JS fit function is stored as a description object, which a choice input cannot display
    const customFitFunction = seriesSource && typeof seriesSource.fitFunction === 'object' ?
      seriesSource.fitFunction : null;
    const fitSeriesChildren = seriesPropertiesFor(chartData, customFitFunction).map((p) => {
      if (p.name !== 'fitFunction' || !customFitFunction)
        return showWholeList(ui.input.forProperty(p, seriesSource, seriesOptionsRefresh));
      // the input is bound to the name, so a pick has to be written back to the series
      const named = {...seriesSource, fitFunction: customFitFunction.name};
      return ui.input.forProperty(p, named, {onValueChanged: (v: any, inputBase: DG.InputBase) => {
        if (v !== customFitFunction.name)
          seriesSource.fitFunction = v;
        seriesOptionsRefresh.onValueChanged(v, inputBase);
      }});
    });
    ui.forms.addGroup(form, 'Series options', fitSeriesChildren);
    // normalize the stored legacy names in place, so the input binds to the real options object
    if (chartData.chartOptions?.showStatistics)
      chartData.chartOptions.showStatistics = normalizeStatisticNames(chartData, chartData.chartOptions.showStatistics);
    ui.forms.addGroup(form, 'Chart options', chartPropertiesFor(chartData).map((p) =>
      showWholeList(ui.input.forProperty(p, chartData.chartOptions, chartOptionsRefresh))));
    acc.addPane('Options', () => form);

    const choices = (chartData.series?.length ?? 0) > 1 ? ['all', 'aggregated'] : ['all'];
    const seriesStatsProperty = DG.Property.js('series', DG.TYPE.STRING,
      {description: 'Controls whether to show series statistics or aggregated statistics',
        defaultValue: 'all', choices: choices, nullable: false});
    const seriesStatsInput = ui.input.forProperty(seriesStatsProperty, null, {onValueChanged: () => {
      acc.getPane('Fit').root.lastElementChild!.replaceChildren(createFitPane());
    }});
    const aggrTypeProperty = DG.Property.js('aggregation type', DG.TYPE.STRING,
      {description: 'Controls which aggregation to use on the series statistics',
        defaultValue: 'med', choices: Object.values(DG.STATS), nullable: false});
    const aggrTypeInput = ui.input.forProperty(aggrTypeProperty, null, {onValueChanged: () => {
      acc.getPane('Fit').root.lastElementChild!.replaceChildren(createFitPane());
    }});

    function createFitPane(): HTMLElement {
      const hostItems = (chartData.series?.length ?? 0) > 1 ? seriesStatsInput.stringValue === 'aggregated' ?
        [seriesStatsInput.root, aggrTypeInput.root] : [seriesStatsInput.root] : null;
      const host = ui.divV(hostItems!);
      const dataBounds = getChartBounds(chartData);
      if (dataBounds.x <= 0 && chartData.chartOptions) chartData.chartOptions.logX = false;
      if (dataBounds.y <= 0 && chartData.chartOptions) chartData.chartOptions.logY = false;

      if (seriesStatsInput.stringValue === 'all') {
        const chartLogOptions: LogOptions = {logX: chartData.chartOptions?.logX, logY: chartData.chartOptions?.logY};
        for (let i = 0; i < chartData.series!.length; i++) {
          const series = chartData.series![i];
          const seriesFit = calculateSeriesFit(series, i, chartLogOptions, tableCell);

          const color = getSeriesColor(series, i, ColorType.FIT_LINE);
          const seriesName = series.name ?? 'series ' + i;
          host.appendChild(ui.panel([
            ui.h1(seriesName, {style: {color: color}}),
            // descriptors come from the series' own fit function, so only applicable statistics are listed
            ui.input.form(seriesFit, getSeriesFitFunction(series).statisticsProperties, {
              onCreated: (input) => input.root.appendChild(ui.iconFA('plus', () =>
                addStatisticColumn(gridCell, 'curveStatistic',
                  {propName: input.property.name, seriesNumber: i}),
              `Calculate ${input.property.name} for the whole column`))
            })
          ]));
        }
      } else {
        const seriesStatistics = getChartDataAggrStats(chartData, aggrTypeInput.stringValue, tableCell);
        host.appendChild(ui.panel([
          ui.h1(`series ${aggrTypeInput.stringValue}`),
          // same per-fit-function descriptors as the single-series pane, deduplicated across series
          ui.input.form(seriesStatistics, aggregatedStatisticsProperties(chartData, aggrTypeInput.stringValue), {
            onCreated: (input) => input.root.appendChild(ui.iconFA('plus', () =>
              addStatisticColumn(gridCell, 'curveAggrStatistic',
                {propName: input.property.name, aggrType: aggrTypeInput.stringValue}),
            `Calculate ${input.property.name} ${aggrTypeInput.stringValue} for the whole column`))
          })
        ]));
      }

      return host;
    }
    const gcw = DG.GridCellWidget.fromGridCell(gridCell);
    gcw.canvas.addEventListener('click', (e) => gcw.render());
    const chartPane = acc.addPane('Chart', () => gcw.root);
    const screenBounds = FitChartCellRenderer.inflateScreenBounds(gridCell.bounds);
    if (screenBounds.width < FitConstants.MIN_POINTS_AND_STATS_VISIBILITY_PX_WIDTH ||
      screenBounds.height < FitConstants.MIN_POINTS_AND_STATS_VISIBILITY_PX_HEIGHT)
      chartPane.expanded = true;

    acc.addPane('Fit', () => createFitPane());

    return acc.root;
  }
}
