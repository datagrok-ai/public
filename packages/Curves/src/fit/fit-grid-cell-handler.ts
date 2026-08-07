/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {
  getSeriesFitFunction,
  getChartBounds,
} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {
  fitChartDataProperties,
  IFitChartData,
  IFitSeries,
  IFitChartOptions,
  IFitSeriesOptions,
  IFitFunctionDescription,
  LogOptions,
} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {isColorValid, FitChartCellRenderer} from './fit-renderer';
import {
  getOrCreateParsedChartData, getColumnChartOptions, getDataFrameChartOptions, mergeProperties,
  substituteZeroes, CHART_OPTIONS, SERIES_OPTIONS,
} from './fit-chart-data';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {parseCellValue, isNativeFormat} from './curve-converter';
import {ColorType, SeriesColorType, getSeriesColor} from './render-utils';
import {calculateSeriesFit, getChartDataAggrStats, aggregatedStatisticsProperties} from './fit-statistics';
import {fitFunctions, fitSeriesProperties, getStatisticProperty} from '@datagrok-libraries/statistics/src/fit/fit-engine';

// Options a statistic can depend on; everything else only repaints. errorModel is here because it
// weights the objective function and so changes the fitted parameters. mergeSeries is not: the
// renderer merges on a copy, and the statistics always read the original series.
const STATISTIC_AFFECTING_OPTIONS = ['logX', 'logY', 'allowXZeroes', 'fitFunction', 'errorModel'];

enum MANIPULATION_LEVEL {
  DATAFRAME = 'Dataframe',
  COLUMN = 'Column',
  CELL = 'Cell'
}

/** Chart properties with `showStatistics` offering the statistics this cell's fit functions produce,
 * rather than the fixed legacy list. */
export function chartPropertiesFor(chartData: IFitChartData): DG.Property[] {
  const names: string[] = [];
  for (const series of chartData.series ?? []) {
    for (const prop of getSeriesFitFunction(series).statisticsProperties) {
      if (!names.includes(prop.name))
        names.push(prop.name);
    }
  }
  if (names.length === 0)
    return fitChartDataProperties;
  return fitChartDataProperties.map((p) => p.name !== 'showStatistics' ? p :
    DG.Property.js('showStatistics', DG.TYPE.STRING_LIST, {description: p.description, choices: names,
      inputType: 'MultiChoice', friendlyName: 'Statistics'}));
}

/** Runs a curve statistic function so that its `join(table)` output action makes the platform add the
 * result as a calculated column - it writes the formula, marks it a vector func so the curve column is
 * bound as a column rather than a per-row value, subscribes it to changes, and (because the call is not
 * processed) appends the AddNewColumn line to the table's creation script. */
type CurveStatisticParams = {propName: string, seriesNumber: number} | {propName: string, aggrType: string};

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

/** `join(table)` appends the column the platform adds, so the position the legacy `addStatisticsColumn`
 * inserted at has to be restored afterwards. The grid only takes a column's position from the dataframe
 * when it first builds a GridColumn for it, and by now it has, so its order is set separately - off its
 * own sequence, which a manual reorder may have taken away from the dataframe's. */
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

/** Maps stored statistic names onto the names the choices use, so an option persisted under a legacy
 * name (`interceptX`) still ticks its current checkbox (`ic50`) instead of silently showing unchecked. */
export function normalizeStatisticNames(chartData: IFitChartData, names: string[]): string[] {
  const resolved: string[] = [];
  for (const name of names) {
    let mapped = name;
    for (const series of chartData.series ?? []) {
      const prop = getStatisticProperty(getSeriesFitFunction(series), name);
      if (prop) {
        mapped = prop.name;
        break;
      }
    }
    if (!resolved.includes(mapped))
      resolved.push(mapped);
  }
  return resolved;
}

/** Series properties with `fitFunction` able to name a custom JS function alongside the built-ins. */
function seriesPropertiesFor(customFitFunction: IFitFunctionDescription | null): DG.Property[] {
  if (!customFitFunction)
    return fitSeriesProperties;
  return fitSeriesProperties.map((p) => p.name !== 'fitFunction' ? p :
    DG.Property.js('fitFunction', DG.TYPE.STRING, {category: 'Fitting',
      choices: [customFitFunction.name, ...Object.keys(fitFunctions)], defaultValue: customFitFunction.name}));
}

type OptionsSection = 'chartOptions' | 'seriesOptions';

/** Records that the user set this option at this level, so it outranks the value the data declares.
 * Without it a column-level choice can only take effect by deleting the data's own value. */
function claim(chartData: IFitChartData, options: string, propertyName: string): void {
  const names = ((chartData.explicit ??= {})[options as OptionsSection] ??= []);
  if (!names.includes(propertyName))
    names.push(propertyName);
}

/** Drops the claim, leaving the value. Returns whether anything was claimed. */
function unclaim(chartData: IFitChartData, options: string, propertyName: string): boolean {
  const names = chartData.explicit?.[options as OptionsSection];
  if (!names?.includes(propertyName))
    return false;
  names.splice(names.indexOf(propertyName), 1);
  return true;
}

function changePlotOptions(chartData: IFitChartData, inputBase: DG.InputBase, options: string): void {
  const propertyName = inputBase.property.name as string;
  if (options === CHART_OPTIONS) {
    if (chartData.chartOptions === undefined) return;
    (chartData.chartOptions[propertyName as keyof IFitChartOptions] as any) = inputBase.value;
  } else if (options === SERIES_OPTIONS) {
    if (chartData.series === undefined) return;
    for (let i = 0; i < chartData.series.length; i++)
      (chartData.series[i][propertyName as keyof IFitSeries] as any) = inputBase.value;
  }
  claim(chartData, options, propertyName);
}

const DETECTED_AT_VERSION = 'fit-overrides-detected-at-version';

/** Records, per property, whether any cell of a curve column claims it - so a Column or Dataframe level
 * change knows whether it has to drop those claims. A value the data declares is not a claim: it loses
 * to an explicitly set level, so there is nothing to clear. Stamped with the column version it was
 * computed from: every cell write bumps that version, so the answer re-derives itself instead of going
 * stale the moment a cell gains a claim. */
function detectSettings(df: DG.DataFrame): void {
  const fitColumns = df.columns.bySemTypeAll(FitConstants.FIT_SEM_TYPE);
  for (let i = 0; i < fitColumns.length; i++) {
    fitColumns[i].temp[DETECTED_AT_VERSION] = fitColumns[i].version;
    fitChartDataProperties.map((prop) => {
      fitColumns[i].temp[`${CHART_OPTIONS}-custom-${prop.name}`] = false;
    });
    fitSeriesProperties.map((prop) => {
      fitColumns[i].temp[`${SERIES_OPTIONS}-custom-${prop.name}`] = false;
    });

    for (let j = 0; j < fitColumns[i].length; j++) {
      if (fitColumns[i].get(j) === '') continue;
      const chartData = parseCellValue(fitColumns[i].get(j), fitColumns[i]);
      for (const name of chartData.explicit?.chartOptions ?? [])
        fitColumns[i].temp[`${CHART_OPTIONS}-custom-${name}`] = true;
      for (const name of chartData.explicit?.seriesOptions ?? [])
        fitColumns[i].temp[`${SERIES_OPTIONS}-custom-${name}`] = true;
    }
  }
}

export function changeCurvesOptions(gridCell: DG.GridCell, inputBase: DG.InputBase, options: string, manipulationLevel: string): void {
  if (gridCell.cell.column.temp[DETECTED_AT_VERSION] !== gridCell.cell.column.version)
    detectSettings(gridCell.cell.dataFrame);
  const propertyName = inputBase.property.name as string;
  const chartOptions = manipulationLevel === MANIPULATION_LEVEL.DATAFRAME ?
    getDataFrameChartOptions(gridCell.cell.dataFrame) : getColumnChartOptions(gridCell.cell.column);
  if (options === CHART_OPTIONS)
    (chartOptions.chartOptions![propertyName as keyof IFitChartOptions] as any) = inputBase.value;
  else if (options === SERIES_OPTIONS)
    (chartOptions.seriesOptions![propertyName as keyof IFitSeriesOptions] as any) = inputBase.value;
  claim(chartOptions, options, propertyName);

  if (manipulationLevel === MANIPULATION_LEVEL.CELL) {
    // Cell-level editing only works for native JSON format
    if (!isNativeFormat(gridCell.cell.column))
      return;
    const value = gridCell.cell.value;
    if (value === '') return;
    const chartData: IFitChartData = JSON.parse(value ?? '{}') ?? {};
    changePlotOptions(chartData, inputBase, options);
    gridCell.cell.value = JSON.stringify(chartData);
  } else {
    let columns: DG.Column[];
    if (manipulationLevel === MANIPULATION_LEVEL.DATAFRAME) {
      gridCell.cell.dataFrame.tags[FitConstants.TAG_FIT] = JSON.stringify(chartOptions);
      columns = gridCell.cell.dataFrame.columns.bySemTypeAll(FitConstants.FIT_SEM_TYPE);
    } else {
      gridCell.cell.column.tags[FitConstants.TAG_FIT] = JSON.stringify(chartOptions);
      columns = [gridCell.cell.column];
    }

    for (let i = 0; i < columns.length; i++) {
      if (manipulationLevel === MANIPULATION_LEVEL.DATAFRAME) {
        const columnChartOptions = getColumnChartOptions(columns[i]);
        options === CHART_OPTIONS ? delete columnChartOptions.chartOptions![propertyName as keyof IFitChartOptions] :
          delete columnChartOptions.seriesOptions![propertyName as keyof IFitSeriesOptions];
        unclaim(columnChartOptions, options, propertyName);
        columns[i].tags[FitConstants.TAG_FIT] = JSON.stringify(columnChartOptions);
      }
      if (columns[i].temp[`${options}-custom-${propertyName}`] === false) continue;
      // Skip cell-level value mutation for non-native format columns
      if (!isNativeFormat(columns[i])) continue;

      // Only the claim is dropped - the value the data declares stays, and loses to this level from
      // now on. Rewritten without notifying: on its own it looks like a data change and would
      // recalculate every dependent statistic column, even for a colour. The notification below is
      // the deliberate one.
      const clearClaim = (j: number): string => {
        const value = columns[i].get(j);
        if (value === '') return value;
        const chartData = (JSON.parse(value ?? '{}') ?? {}) as IFitChartData;
        return unclaim(chartData, options, propertyName) ? JSON.stringify(chartData) : value;
      };
      for (let j = 0; j < columns[i].length; j++) {
        const updated = clearClaim(j);
        if (updated !== columns[i].get(j))
          columns[i].set(j, updated, false);
      }
      columns[i].temp[`${options}-custom-${propertyName}`] = false;
    }

    // the option lives in a tag rather than in the data, so nothing marks the curve column changed
    // and statistic columns extracted from it would keep stale numbers while the plot updates.
    // Only for options a statistic can depend on - a title edit would otherwise refit every column.
    const affectsStatistics = STATISTIC_AFFECTING_OPTIONS.includes(propertyName);
    if (affectsStatistics) {
      for (const column of columns)
        column.fireValuesChanged();
    }
    // the option and any claim the rewrite above dropped are both stored, but nothing has marked the
    // table modified - a tag write raises metadata, not data, and the rewrite is silent. This
    // notification carries no column, so the calculated columns' subscription, which filters on the
    // changed column, still ignores it.
    if (!affectsStatistics)
      gridCell.cell.dataFrame.fireValuesChanged();
  }
  gridCell.grid.invalidate();
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
      changeCurvesOptions(gridCell, inputBase, SERIES_OPTIONS, switchLevelInput.value)};
    const chartOptionsRefresh = {onValueChanged: (v: any, inputBase: DG.InputBase) =>
      changeCurvesOptions(gridCell, inputBase, CHART_OPTIONS, switchLevelInput.value)};

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
    const fitSeriesChildren = seriesPropertiesFor(customFitFunction).map((p) => {
      if (p.name !== 'fitFunction' || !customFitFunction)
        return ui.input.forProperty(p, seriesSource, seriesOptionsRefresh);
      // the choice input cannot hold the description object, so it is bound to its name - picking
      // anything else has to be written back to the series, or the input does nothing
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
      ui.input.forProperty(p, chartData.chartOptions, chartOptionsRefresh)));
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
