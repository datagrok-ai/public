import * as DG from 'datagrok-api/dg';

import {fitFunctions} from '@datagrok-libraries/statistics/src/fit/fit-engine';
import {
  IFitChartData,
  IFitChartOptions,
  IFitSeriesOptions,
  IFitSeries,
  IFitPoint,
  FIT_FUNCTION_SIGMOID,
  FitFunctionName,
} from '@datagrok-libraries/statistics/src/fit/fit-curve';


const AXES = {x: 'xAxis', y: 'yAxis'};
const EXTREMUMS = {min: 'min', max: 'max'};

/** XML attributes are strings, so `!!element.getAttribute(name)` is true for "False" as well. */
function xmlBool(element: Element | undefined, name: string): boolean {
  return element?.getAttribute(name)?.toLowerCase() === 'true';
}


/** Constructs {@link IFitChartOptions} from grid and settings xml tags.
 * @param {Element} grid XML grid tag
 * @param {Element} settings XML settings tag
 * @return {IFitChartOptions} IFitChartOptions for the fitted curve
*/
function getChartOptions(grid: Element, settings: Element): IFitChartOptions {
  return {
    minX: +grid.getElementsByTagName(AXES.x)[0].getAttribute(EXTREMUMS.min)!,
    minY: +grid.getElementsByTagName(AXES.y)[0].getAttribute(EXTREMUMS.min)!,
    maxX: +grid.getElementsByTagName(AXES.x)[0].getAttribute(EXTREMUMS.max)!,
    maxY: +grid.getElementsByTagName(AXES.y)[0].getAttribute(EXTREMUMS.max)!,

    xAxisName: settings.getAttribute('xLabel')!,
    yAxisName: settings.getAttribute('yLabel')!,
    logX: xmlBool(settings, 'logX'),
    // only when declared: an always-present key owns it, and a Column or Dataframe level option
    // can no longer cascade onto the cell
    ...(settings.getAttribute('logY') !== null ? {logY: xmlBool(settings, 'logY')} : {}),
  };
}

/** Constructs {@link IFitSeriesOptions} from the series xml tag.
 * @param {Element} series XML series tag
 * @return {IFitSeriesOptions} IFitSeriesOptions for the fitted curve
*/
function getSeriesOptions(series: Element): IFitSeriesOptions {
  const params = (series.getElementsByTagName('params')[0]?.childNodes[0].nodeValue)?.split(',')!.map(Number)!;
  const declaredType = series.getElementsByTagName('function')[0]?.getAttribute('type');
  // 'sigif' is this format's name for a sigmoid; anything we cannot fit is better rendered as one
  // than passed through as a name no fit function is registered under
  const funcType = !declaredType || declaredType === 'sigif' || !(declaredType in fitFunctions) ?
    FIT_FUNCTION_SIGMOID : declaredType;
  const markerColor = series.getElementsByTagName('settings')[0].getAttribute('markerColor')!;
  const lineColor = series.getElementsByTagName('settings')[0].getAttribute('color')!;
  const drawLine = xmlBool(series.getElementsByTagName('settings')[0], 'drawLine');
  const seriesName = series.getAttribute('name')!;

  const seriesOptions: IFitSeriesOptions = {
    name: seriesName,
    fitFunction: funcType as FitFunctionName,
    markerType: DG.MARKER_TYPE.CIRCLE,
    pointColor: markerColor,
    fitLineColor: lineColor,
    showFitLine: drawLine,
    showPoints: 'points',
    showCurveConfidenceInterval: false,
    clickToToggle: false
  };

  // params there are: [IC50, tan, max, min] - so we place them correctly: [max, tan, IC50, min]
  if (params)
    seriesOptions.parameters = [params[2], params[1], params[0], params[3]];

  return seriesOptions;
}

/** Constructs {@link IFitPoint} array from the grid series tag.
 * @param {Element} series XML series tag
 * @return {IFitPoint[]} IFitPoint array for current series
*/
function getPoints(series: Element): IFitPoint[] {
  const xCoords = (series.getElementsByTagName('x')[0].childNodes[0].nodeValue)?.split(',')!;
  const yCoords = (series.getElementsByTagName('y')[0].childNodes[0].nodeValue)?.split(',')!;
  const mask = (series.getElementsByTagName('mask')[0].childNodes[0].nodeValue)?.split('')!;

  const points: IFitPoint[] = [];
  for (let j = 0; j < xCoords.length; j++) {
    points[j] = {
      x: +xCoords[j],
      y: +yCoords[j],
      outlier: mask[j] === '0',
    };
  }

  return points;
}

/** Constructs {@link IFitSeries} from the grid series tag.
 * @param {Element} series XML series tag
 * @return {IFitSeries} IFitSeries for the current series
*/
function getSeries(series: Element): IFitSeries {
  const points: IFitPoint[] = getPoints(series);
  const currentFitSeriesOptions: IFitSeriesOptions = getSeriesOptions(series);

  const returnSeries: IFitSeries = {
    points: points,
  };
  Object.assign(returnSeries, currentFitSeriesOptions);

  return returnSeries;
}

/** Constructs {@link IFitSeries} array from the seriesCollection xml tag.
 * @param {Element} seriesCollection XML seriesCollection tag
 * @return {IFitSeries[]} IFitSeries array for the fitted curve
*/
function getSeriesArray(seriesCollection: Element): IFitSeries[] {
  const fitSeries: IFitSeries[] = [];

  const seriesCollectionLength = seriesCollection.getElementsByTagName('series').length;
  for (let i = 0; i < seriesCollectionLength; i++) {
    const currentSeries = seriesCollection.getElementsByTagName('series')[i];
    fitSeries[i] = getSeries(currentSeries);
  }

  return fitSeries;
}


/** Converts XML fitted curve chart document into {@link IFitChartData} interface.
 * @param {string} xmlText XML document
 * @return {IFitChartData} IFitChartData interface for the fitted curve
*/
export function convertXMLToIFitChartData(xmlText: string): IFitChartData {
  const parser = new DOMParser();
  const xmlDoc = parser.parseFromString(xmlText, 'text/xml');

  // get IFitChartOptions from grid
  const gridElement = xmlDoc.getElementsByTagName('grid')[0];
  const settingsElement = xmlDoc.getElementsByTagName('settings')[0];
  const fitChartOptions: IFitChartOptions = getChartOptions(gridElement, settingsElement);

  // get the whole series collection
  const seriesCollection = xmlDoc.getElementsByTagName('seriesCollection')[0];

  // get first series to make it the default options (IFitSeriesOptions)
  const firstSeries = seriesCollection.getElementsByTagName('series')[0];
  const fitSeriesOptions: IFitSeriesOptions = getSeriesOptions(firstSeries);

  // get IFitSeries[]
  const fitSeries = getSeriesArray(seriesCollection);

  return {
    chartOptions: fitChartOptions,
    seriesOptions: fitSeriesOptions,
    series: fitSeries
  };
}
