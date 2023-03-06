import {
  IFitChartData,
  IFitChartOptions,
  IFitSeriesOptions,
  IFitSeries,
  IFitPoint
} from './fit-data';

const AXES = {x: 'xAxis', y: 'yAxis'};
const EXTREMUMS = {min: 'min', max: 'max'};


/** Constructs {@link IFitChartOptions} from the grid xml tag.
 * @param {Element} grid XML grid tag
 * @return {IFitChartOptions} IFitChartOptions for the fitted curve
*/
function getChartOptions(grid: Element): IFitChartOptions {
  const fitChartOptions: IFitChartOptions = {
    minX: +grid.getElementsByTagName(AXES.x)[0].getAttribute(EXTREMUMS.min)!,
    minY: +grid.getElementsByTagName(AXES.y)[0].getAttribute(EXTREMUMS.min)!,
    maxX: +grid.getElementsByTagName(AXES.x)[0].getAttribute(EXTREMUMS.max)!,
    maxY: +grid.getElementsByTagName(AXES.y)[0].getAttribute(EXTREMUMS.max)!,
  };

  return fitChartOptions;
}

/** Constructs {@link IFitSeriesOptions} from the series xml tag.
 * @param {Element} series XML series tag
 * @return {IFitSeriesOptions} IFitSeriesOptions for the fitted curve
*/
function getSeriesOptions(series: Element): IFitSeriesOptions {
  const params = (series.getElementsByTagName('params')[0].childNodes[0].nodeValue)?.split(',')!.map(Number)!;
  let funcType = series.getElementsByTagName('function')[0].getAttribute('type')!;
  funcType = funcType === 'sigif' ? 'Sigmoid': funcType;
  const markerColor = series.getElementsByTagName('settings')[0].getAttribute('markerColor')!;
  const lineColor = series.getElementsByTagName('settings')[0].getAttribute('color')!;
  const drawLine = !!series.getElementsByTagName('settings')[0].getAttribute('drawLine')!;
  const seriesName = series.getAttribute('name')!;

  const fitSeriesOptions: IFitSeriesOptions = {
    parameters: params,
    fitFunction: funcType,
    pointColor: markerColor,
    fitLineColor: lineColor,
    showFitLine: drawLine,
    name: seriesName,
  };

  return fitSeriesOptions;
}

/** Constructs {@link IFitPoint} array from the grid series tag.
 * @param {Element} series XML series tag
 * @return {IFitPoint[]} IFitPoint array for current series
*/
function getPoints(series: Element): IFitPoint[] {
  const xCoords = (series.getElementsByTagName('x')[0].childNodes[0].nodeValue)?.split(',')!;
  const yCoords = (series.getElementsByTagName('y')[0].childNodes[0].nodeValue)?.split(',')!;
  // const mask = (series.getElementsByTagName('mask')[0].childNodes[0].nodeValue)?.split('')!;

  const points: IFitPoint[] = [];
  for (let j = 0; j < xCoords.length; j++) {
    points[j] = {
      x: +xCoords[j],
      y: +yCoords[j],
      // outlier: !!mask[j],
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
  Object.assign(series, currentFitSeriesOptions);

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
 * @return {IFitChartData} IFitChartData intefrace for the fitted curve
*/
export function convertXMLToIFitChartData(xmlText: string): IFitChartData {
  const parser = new DOMParser();
  const xmlDoc = parser.parseFromString(xmlText, 'text/xml');

  // get IFitChartOptions from grid
  const gridElement = xmlDoc.getElementsByTagName('grid')[0];
  const fitChartOptions: IFitChartOptions = getChartOptions(gridElement);

  // get the whole series collection
  const seriesCollection = xmlDoc.getElementsByTagName('seriesCollection')[0];

  // get first series to make it the default options (IFitSeriesOptions)
  const firstSeries = seriesCollection.getElementsByTagName('series')[0];
  const fitSeriesOptions: IFitSeriesOptions = getSeriesOptions(firstSeries);

  // get IFitSeries[]
  const fitSeries = getSeriesArray(seriesCollection);

  const fitChartData: IFitChartData = {
    chartOptions: fitChartOptions,
    seriesOptions: fitSeriesOptions,
    series: fitSeries
  };

  return fitChartData;
}
