/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {
  getSeriesFit,
  getSeriesFitFunction,
  toFitStatistics,
  toDataSpace,
  LogOptions,
  getChartBounds,
} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {
  statisticsProperties,
  fitChartDataProperties,
  IFitChartData,
  IFitSeries,
  IFitChartOptions,
  IFitSeriesOptions,
  IFitFunctionDescription,
  FitStatistics,
  LEGACY_FIT_STATISTICS,
} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {
  getOrCreateParsedChartData,
  getColumnChartOptions,
  getDataFrameChartOptions,
  isColorValid,
  mergeProperties,
  substituteZeroes, getOrCreateCachedFitCurve, getOrCreateCachedCurvesDataPoints, FitChartCellRenderer
} from './fit-renderer';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {parseCellValue, isNativeFormat} from './curve-converter';
import {ColorType, SeriesColorType, getSeriesColor} from './render-utils';
import {Fit, FitFunction, fitFunctions, fitSeriesProperties, getStatistic, getStatisticProperty} from '@datagrok-libraries/statistics/src/fit/new-fit-API';


const CHART_OPTIONS = 'chartOptions';
const SERIES_OPTIONS = 'seriesOptions';
enum MANIPULATION_LEVEL {
  DATAFRAME = 'Dataframe',
  COLUMN = 'Column',
  CELL = 'Cell'
}
const AGGREGATION_TYPES: {[key: string]: string} = {
  'count': 'totalCount',
  'nulls': 'missingValueCount',
  'unique': 'uniqueCount',
  'values': 'valueCount',
  'min': 'min',
  'max': 'max',
  'sum': 'sum',
  'avg': 'avg',
  'stdev': 'stdev',
  'variance': 'variance',
  'skew': 'skew',
  'kurt': 'kurt',
  'med': 'med',
  'q1': 'q1',
  'q2': 'q2',
  'q3': 'q3'
};


/** Returns the typed fit for a series, with the inflection point reported in data space. */
export function calculateSeriesFit(series: IFitSeries, seriesIdx: number, chartLogOptions: LogOptions,
  tableCell?: DG.Cell, useCache: boolean = true, inDataSpace: boolean = true): Fit {
  const fitFunction = getSeriesFitFunction(series);
  if (series.parameters) {
    if (chartLogOptions.logX) {
      if (series.parameters[2] > 0)
        series.parameters[2] = Math.log10(series.parameters[2]);
    }
  } else {
    const params = getOrCreateCachedFitCurve(series, seriesIdx, fitFunction, chartLogOptions, tableCell, useCache).parameters;
    series.parameters = [...params];
  }

  const fit = getSeriesFit(series, fitFunction,
    getOrCreateCachedCurvesDataPoints(series, seriesIdx, chartLogOptions, false, tableCell), chartLogOptions);
  return inDataSpace ? toDataSpace(fit, chartLogOptions) : fit;
}

/** Returns series statistics in the legacy shape. Prefer {@link calculateSeriesFit}. */
export function calculateSeriesStats(series: IFitSeries, seriesIdx: number, chartLogOptions: LogOptions,
  tableCell?: DG.Cell, useCache: boolean = true): FitStatistics {
  return toFitStatistics(calculateSeriesFit(series, seriesIdx, chartLogOptions, tableCell, useCache));
}

/** Chart properties with `showStatistics` offering the statistics this cell's fit functions produce,
 * rather than the fixed legacy list. */
function chartPropertiesFor(chartData: IFitChartData): DG.Property[] {
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

async function addStatisticColumn(gridCell: DG.GridCell, funcName: string,
  params: CurveStatisticParams): Promise<void> {
  await DG.Func.find({name: funcName})[0]
    .prepare({table: gridCell.cell.dataFrame, curveColumn: gridCell.cell.column, ...params})
    .call(false, undefined, {processed: false});
}

/** Maps stored statistic names onto the names the choices use, so an option persisted under a legacy
 * name (`interceptX`) still ticks its current checkbox (`ic50`) instead of silently showing unchecked. */
function normalizeStatisticNames(chartData: IFitChartData, names: string[]): string[] {
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

export type AggregatedFitStatistics = FitStatistics & {[name: string]: number | undefined};

/** Statistics descriptors covering every fit function used in the cell, without duplicates. */
function aggregatedStatisticsProperties(chartData: IFitChartData): DG.Property[] {
  const props: DG.Property[] = [];
  for (const series of chartData.series ?? []) {
    for (const prop of getSeriesFitFunction(series).statisticsProperties) {
      if (!props.some((p) => p.name === prop.name))
        props.push(prop);
    }
  }
  return props.length ? props : statisticsProperties;
}

/** Aggregates across all series every statistic that the cell's fit functions produce - not just the
 * seven legacy ones. Values are aggregated in fit space and converted once at the end, so an averaged
 * IC50 is a geometric mean. Legacy names stay available so recorded transforms keep resolving. */
export function getChartDataAggrStats(chartData: IFitChartData, aggrType: string,
  tableCell?: DG.Cell): AggregatedFitStatistics {
  const chartLogOptions: LogOptions = {logX: chartData.chartOptions?.logX, logY: chartData.chartOptions?.logY};
  const values: Map<string, (number | undefined)[]> = new Map();
  const fitFunctionsUsed: FitFunction[] = [];

  for (let i = 0; i < chartData.series?.length!; i++) {
    const series = chartData.series![i];
    if (series.points.every((p) => p.outlier))
      continue;
    const fitFunction = getSeriesFitFunction(series);
    fitFunctionsUsed.push(fitFunction);
    const fit = calculateSeriesFit(series, i, chartLogOptions, tableCell, true, false);
    for (const prop of fitFunction.statisticsProperties) {
      if (!values.has(prop.name))
        values.set(prop.name, []);
      values.get(prop.name)!.push(getStatistic(fit, prop.name));
    }
  }

  const aggregated: AggregatedFitStatistics = {};
  for (const [name, seriesValues] of values) {
    aggregated[name] = seriesValues.some((v) => v === undefined || v === null || isNaN(v!)) ? undefined :
      DG.Stats.fromValues(seriesValues as number[])[AGGREGATION_TYPES[aggrType] as keyof DG.Stats] as number;
  }

  // convert once, after aggregating - the same boundary a single series goes through
  toDataSpace(aggregated as unknown as Fit, chartLogOptions);

  // legacy names resolve to whatever the series' fit functions call them today
  for (const legacyName of LEGACY_FIT_STATISTICS) {
    if (aggregated[legacyName] !== undefined)
      continue;
    for (const fitFunction of fitFunctionsUsed) {
      const prop = getStatisticProperty(fitFunction, legacyName);
      if (prop && aggregated[prop.name] !== undefined) {
        aggregated[legacyName] = aggregated[prop.name];
        break;
      }
    }
  }
  return aggregated;
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
}

function detectSettings(df: DG.DataFrame): void {
  const fitColumns = df.columns.bySemTypeAll(FitConstants.FIT_SEM_TYPE);
  for (let i = 0; i < fitColumns.length; i++) {
    fitChartDataProperties.map((prop) => {
      fitColumns[i].temp[`${CHART_OPTIONS}-custom-${prop.name}`] = false;
    });
    fitSeriesProperties.map((prop) => {
      fitColumns[i].temp[`${SERIES_OPTIONS}-custom-${prop.name}`] = false;
    });

    for (let j = 0; j < fitColumns[i].length; j++) {
      if (fitColumns[i].get(j) === '') continue;
      const chartData = parseCellValue(fitColumns[i].get(j), fitColumns[i]);

      fitChartDataProperties.map((prop) => {
        if (!chartData.chartOptions) return;
        if (chartData.chartOptions[prop.name as keyof IFitChartOptions] !== undefined)
          fitColumns[i].temp[`${CHART_OPTIONS}-custom-${prop.name}`] = true;
      });

      fitSeriesProperties.map((prop) => {
        if (!chartData.series) return;
        for (const series of chartData.series) {
          if (series[prop.name as keyof IFitSeriesOptions] !== undefined)
            fitColumns[i].temp[`${SERIES_OPTIONS}-custom-${prop.name}`] = true;
        }
      });
    }
  }
}

function changeCurvesOptions(gridCell: DG.GridCell, inputBase: DG.InputBase, options: string, manipulationLevel: string): void {
  if (gridCell.cell.column.temp[`${CHART_OPTIONS}-custom-title`] === undefined)
    detectSettings(gridCell.cell.dataFrame);
  const propertyName = inputBase.property.name as string;
  const chartOptions = manipulationLevel === MANIPULATION_LEVEL.DATAFRAME ?
    getDataFrameChartOptions(gridCell.cell.dataFrame) : getColumnChartOptions(gridCell.cell.column);
  if (options === CHART_OPTIONS)
    (chartOptions.chartOptions![propertyName as keyof IFitChartOptions] as any) = inputBase.value;
  else if (options === SERIES_OPTIONS)
    (chartOptions.seriesOptions![propertyName as keyof IFitSeriesOptions] as any) = inputBase.value;

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
        columns[i].tags[FitConstants.TAG_FIT] = JSON.stringify(columnChartOptions);
      }
      if (columns[i].temp[`${options}-custom-${propertyName}`] === false) continue;
      // Skip cell-level value mutation for non-native format columns
      if (!isNativeFormat(columns[i])) continue;

      columns[i].init((j) => {
        const value = columns[i].get(j);
        if (value === '') return value;
        const chartData = (JSON.parse(columns[i].get(j) ?? '{}') ?? {}) as IFitChartData;
        if (options === CHART_OPTIONS) {
          if (chartData.chartOptions === undefined) return value;
          if (chartData.chartOptions[propertyName as keyof IFitChartOptions] === undefined)
            return value;
          delete chartData.chartOptions[propertyName as keyof IFitChartOptions];
        } else {
          if (chartData.series === undefined) return value;
          let isSeriesChanged = false;
          for (const series of chartData.series) {
            if (series[propertyName as keyof IFitSeriesOptions] !== undefined) {
              delete series[propertyName as keyof IFitSeriesOptions];
              isSeriesChanged = true;
            }
          }
          if (chartData.seriesOptions)
            delete chartData.seriesOptions[propertyName as keyof IFitSeriesOptions];
          if (!isSeriesChanged) return value;
        }
        return JSON.stringify(chartData);
      });
      columns[i].temp[`${options}-custom-${propertyName}`] = false;
    }
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
    const fitSeriesChildren = seriesPropertiesFor(customFitFunction).map((p) =>
      ui.input.forProperty(p, p.name === 'fitFunction' && customFitFunction ?
        {...seriesSource, fitFunction: customFitFunction.name} : seriesSource, seriesOptionsRefresh));
    ui.forms.addGroup(form, 'Series options', fitSeriesChildren);
    ui.forms.addGroup(form, 'Chart options', chartPropertiesFor(chartData).map((p) =>
      ui.input.forProperty(p, p.name === 'showStatistics' ? {...chartData.chartOptions,
        showStatistics: normalizeStatisticNames(chartData, chartData.chartOptions?.showStatistics ?? [])} :
        chartData.chartOptions, chartOptionsRefresh)));
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
          ui.input.form(seriesStatistics, aggregatedStatisticsProperties(chartData), {
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
