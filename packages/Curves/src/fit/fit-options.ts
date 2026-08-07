/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';

import {getSeriesFitFunction} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {
  fitChartDataProperties,
  IFitChartData,
  IFitSeries,
  IFitFunctionDescription,
} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {fitFunctions, fitSeriesProperties, getStatisticProperty} from '@datagrok-libraries/statistics/src/fit/fit-engine';
import {getColumnChartOptions, getDataFrameChartOptions, CHART_OPTIONS, SERIES_OPTIONS} from './fit-chart-data';
import {isNativeFormat} from './curve-converter';

// Options a statistic can depend on; everything else only repaints. errorModel weights the objective
// function; mergeSeries does not, since the statistics read the original series.
const STATISTIC_AFFECTING_OPTIONS = ['logX', 'logY', 'allowXZeroes', 'fitFunction', 'errorModel'];

export enum MANIPULATION_LEVEL {
  DATAFRAME = 'Dataframe',
  COLUMN = 'Column',
  CELL = 'Cell'
}

export type ManipulationLevel = `${MANIPULATION_LEVEL}`;
export type OptionsSection = typeof CHART_OPTIONS | typeof SERIES_OPTIONS;

/** Chart properties with `showStatistics` offering what this cell's fit functions actually produce. */
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

/** Maps stored statistic names onto the current ones, so `interceptX` still ticks the `ic50` box. */
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
export function seriesPropertiesFor(customFitFunction: IFitFunctionDescription | null): DG.Property[] {
  if (!customFitFunction)
    return fitSeriesProperties;
  return fitSeriesProperties.map((p) => p.name !== 'fitFunction' ? p :
    DG.Property.js('fitFunction', DG.TYPE.STRING, {category: 'Fitting',
      choices: [customFitFunction.name, ...Object.keys(fitFunctions)], defaultValue: customFitFunction.name}));
}

/** Records that the user set this option here, so it outranks the value the data declares. */
function claim(chartData: IFitChartData, options: OptionsSection, propertyName: string): void {
  const names = ((chartData.explicit ??= {})[options] ??= []);
  if (!names.includes(propertyName))
    names.push(propertyName);
}

/** Drops the claim, leaving the value. Returns whether anything was claimed. */
function unclaim(chartData: IFitChartData, options: OptionsSection, propertyName: string): boolean {
  const names = chartData.explicit?.[options];
  if (!names?.includes(propertyName))
    return false;
  names.splice(names.indexOf(propertyName), 1);
  return true;
}

function changePlotOptions(chartData: IFitChartData, inputBase: DG.InputBase, options: OptionsSection): void {
  const propertyName = inputBase.property.name as string;
  if (options === CHART_OPTIONS) {
    if (chartData.chartOptions === undefined) return;
    (chartData.chartOptions as any)[propertyName] = inputBase.value;
  } else {
    if (chartData.series === undefined) return;
    for (const series of chartData.series)
      (series as IFitSeries as any)[propertyName] = inputBase.value;
  }
  claim(chartData, options, propertyName);
}

export function changeCurvesOptions(gridCell: DG.GridCell, inputBase: DG.InputBase, options: OptionsSection,
  manipulationLevel: ManipulationLevel): void {
  const propertyName = inputBase.property.name as string;
  const chartOptions = manipulationLevel === MANIPULATION_LEVEL.DATAFRAME ?
    getDataFrameChartOptions(gridCell.cell.dataFrame) : getColumnChartOptions(gridCell.cell.column);
  ((chartOptions[options] ??= {}) as any)[propertyName] = inputBase.value;
  claim(chartOptions, options, propertyName);

  if (manipulationLevel === MANIPULATION_LEVEL.CELL) {
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

    for (const column of columns) {
      if (manipulationLevel === MANIPULATION_LEVEL.DATAFRAME) {
        const columnChartOptions = getColumnChartOptions(column);
        const section = columnChartOptions[options];
        if (section)
          delete (section as any)[propertyName];
        unclaim(columnChartOptions, options, propertyName);
        column.tags[FitConstants.TAG_FIT] = JSON.stringify(columnChartOptions);
      }
      if (!isNativeFormat(column))
        continue;

      // only the claim is dropped, and silently - a rewrite reads as a data change and would
      // recalculate every dependent statistic column. The notification below is the deliberate one.
      for (let j = 0; j < column.length; j++) {
        const value = column.get(j);
        if (value === '')
          continue;
        const chartData = (JSON.parse(value) ?? {}) as IFitChartData;
        if (unclaim(chartData, options, propertyName))
          column.set(j, JSON.stringify(chartData), false);
      }
    }

    // the option lives in a tag, so nothing else marks the curve column changed
    const affectsStatistics = STATISTIC_AFFECTING_OPTIONS.includes(propertyName);
    if (affectsStatistics) {
      for (const column of columns)
        column.fireValuesChanged();
    }
    // nothing has marked the table modified yet: a tag write raises metadata, not data. Carrying no
    // column keeps the calculated columns' subscription, which filters on one, from firing.
    if (!affectsStatistics)
      gridCell.cell.dataFrame.fireValuesChanged();
  }
  gridCell.grid.invalidate();
}
