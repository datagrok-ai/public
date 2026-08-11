/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {getSeriesFitFunction} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {
  fitChartDataProperties,
  IFitChartData,
  IFitSeries,
  IFitFunctionDescription,
} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {fitFunctions, fitFunctionDescriptions, fitSeriesProperties, getStatisticProperty, DEFAULT_FIT_FUNCTION} from '@datagrok-libraries/statistics/src/fit/fit-engine';
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

/** Label names this cell carries, plot-level ones first. */
function labelNames(chartData: IFitChartData): string[] {
  const names = Object.keys(chartData.chartOptions?.labels ?? {});
  for (const series of chartData.series ?? []) {
    for (const name of Object.keys(series.labels ?? {})) {
      if (!names.includes(name))
        names.push(name);
    }
  }
  return names;
}

/** Chart properties with `showStatistics` and `showLabels` offering what this cell actually carries. */
export function chartPropertiesFor(chartData: IFitChartData): DG.Property[] {
  const names: string[] = [];
  for (const series of chartData.series ?? []) {
    for (const prop of getSeriesFitFunction(series).statisticsProperties) {
      if (!names.includes(prop.name))
        names.push(prop.name);
    }
  }
  const labels = labelNames(chartData);
  // an option with nothing to choose from is an empty row in the panel, and a single curve has
  // nothing to summarise
  const several = (chartData.series?.length ?? 0) > 1;
  return fitChartDataProperties.filter((p) => (p.name !== 'showLabels' || labels.length > 0) &&
    ((p.name !== 'statisticsMode' && p.name !== 'aggrType') || several)).map((p) => {
    if (p.name === 'showStatistics' && names.length > 0) {
      return DG.Property.js('showStatistics', DG.TYPE.STRING_LIST, {description: p.description, choices: names,
        inputType: 'MultiChoice', friendlyName: 'Statistics'});
    }
    if (p.name === 'showLabels') {
      return DG.Property.js('showLabels', DG.TYPE.STRING_LIST, {description: p.description, choices: labels,
        inputType: 'MultiChoice', friendlyName: 'Labels'});
    }
    return p;
  });
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

/** Series properties with `fitFunction` naming every registered function, custom ones included. */
export function seriesPropertiesFor(chartData: IFitChartData,
  customFitFunction: IFitFunctionDescription | null): DG.Property[] {
  // an ICxx name only resolves for a curve that levels off at both ends
  const properties = fitSeriesProperties.filter((p) => p.name !== 'droplines' ||
    (chartData.series ?? []).some((series) => getSeriesFitFunction(series).hasAsymptotes));
  // a custom function registers itself the first time it is used, so the registry is what there is to
  // pick from - read here rather than off fitSeriesProperties, which captured it before any registered
  const choices = [...new Set([...(customFitFunction ? [customFitFunction.name] : []),
    ...Object.keys(fitFunctions)])];
  return properties.map((p) => p.name !== 'fitFunction' ? p :
    DG.Property.js('fitFunction', DG.TYPE.STRING, {category: 'Fitting', choices: choices,
      nullable: false, defaultValue: customFitFunction?.name ?? DEFAULT_FIT_FUNCTION}));
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

/** A custom fit function is stored as the notation that describes it, not as its name: the name only
 * resolves while something else in the session carries the function, which a saved table may not. */
function storedValue(propertyName: string, value: any): any {
  return propertyName === 'fitFunction' && typeof value === 'string' && fitFunctionDescriptions[value] ?
    fitFunctionDescriptions[value] : value;
}

function changePlotOptions(chartData: IFitChartData, inputBase: DG.InputBase, options: OptionsSection): void {
  const propertyName = inputBase.property.name as string;
  const value = storedValue(propertyName, inputBase.value);
  if (options === CHART_OPTIONS) {
    if (chartData.chartOptions === undefined) return;
    (chartData.chartOptions as any)[propertyName] = value;
  } else {
    if (chartData.series === undefined) return;
    for (const series of chartData.series)
      (series as IFitSeries as any)[propertyName] = value;
  }
  claim(chartData, options, propertyName);
}

export function changeCurvesOptions(gridCell: DG.GridCell, inputBase: DG.InputBase, options: OptionsSection,
  manipulationLevel: ManipulationLevel): void {
  const propertyName = inputBase.property.name as string;
  const affectsStatistics = STATISTIC_AFFECTING_OPTIONS.includes(propertyName);

  if (manipulationLevel === MANIPULATION_LEVEL.CELL) {
    if (!isNativeFormat(gridCell.cell.column))
      return;
    const value = gridCell.cell.value;
    if (value === '') return;
    let chartData: IFitChartData;
    try {
      chartData = JSON.parse(value ?? '{}') ?? {};
    } catch (e) {
      grok.shell.error(`Curves: this cell holds malformed data, so the option was not applied: ${e}`);
      return;
    }
    changePlotOptions(chartData, inputBase, options);
    // notifying carries this row's index, so only its statistic recalculates - worth it when the
    // option can change one, wasted work when it cannot
    gridCell.cell.column.set(gridCell.cell.rowIndex, JSON.stringify(chartData), affectsStatistics);
    if (!affectsStatistics)
      gridCell.cell.dataFrame.fireValuesChanged();
    gridCell.grid.invalidate();
    return;
  }

  const chartOptions = manipulationLevel === MANIPULATION_LEVEL.DATAFRAME ?
    getDataFrameChartOptions(gridCell.cell.dataFrame) : getColumnChartOptions(gridCell.cell.column);
  ((chartOptions[options] ??= {}) as any)[propertyName] = storedValue(propertyName, inputBase.value);
  claim(chartOptions, options, propertyName);

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
      let chartData: IFitChartData;
      // a malformed row keeps whatever it has rather than aborting the sweep with the tag already written
      try {
        chartData = (JSON.parse(value) ?? {}) as IFitChartData;
      } catch (_) {
        continue;
      }
      if (unclaim(chartData, options, propertyName))
        column.set(j, JSON.stringify(chartData), false);
    }
  }

  // the option lives in a tag, so nothing else marks the curve column changed
  if (affectsStatistics) {
    for (const column of columns)
      column.fireValuesChanged();
  }
  // nothing has marked the table modified yet: a tag write raises metadata, not data. Carrying no
  // column keeps the calculated columns' subscription, which filters on one, from firing.
  if (!affectsStatistics)
    gridCell.cell.dataFrame.fireValuesChanged();
  gridCell.grid.invalidate();
}
