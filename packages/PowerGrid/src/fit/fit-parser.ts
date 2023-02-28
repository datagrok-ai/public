import {IFitChartData, IFitChartOptions, IFitSeriesOptions, IFitSeries, IFitPoint} from './fit-data';

const AXES = {x: 'xAxis', y: 'yAxis'};
const EXTREMUMS = {min: 'min', max: 'max'};


function getChartOptions(grid: Element): IFitChartOptions {
  const fitChartOptions: IFitChartOptions = {
    minX: +grid.getElementsByTagName(AXES.x)[0].getAttribute(EXTREMUMS.min)!,
    minY: +grid.getElementsByTagName(AXES.y)[0].getAttribute(EXTREMUMS.min)!,
    maxX: +grid.getElementsByTagName(AXES.x)[0].getAttribute(EXTREMUMS.max)!,
    maxY: +grid.getElementsByTagName(AXES.y)[0].getAttribute(EXTREMUMS.max)!,
  };

  return fitChartOptions;
}

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

function getSeries(series: Element): IFitSeries {
  const points: IFitPoint[] = getPoints(series);
  const currentFitSeriesOptions: IFitSeriesOptions = getSeriesOptions(series);

  const returnSeries: IFitSeries = {
    points: points,
  };
  Object.assign(series, currentFitSeriesOptions);

  return returnSeries;
}


export function convertXMLToIFitChartData(xmlText: string): IFitChartData {
  const parser = new DOMParser();
  const xmlDoc = parser.parseFromString(xmlText, 'text/xml');

  const gridElement = xmlDoc.getElementsByTagName('grid')[0];
  const fitChartOptions: IFitChartOptions = getChartOptions(gridElement);

  // get first series to make it default options
  const firstSeries = xmlDoc.getElementsByTagName('series')[0];
  const fitSeriesOptions: IFitSeriesOptions = getSeriesOptions(firstSeries);

  const fitSeries: IFitSeries[] = [];
  for (let i = 0; i < xmlDoc.getElementsByTagName('series').length; i++) {
    const currentSeries = xmlDoc.getElementsByTagName('series')[i];
    fitSeries[i] = getSeries(currentSeries);
  }

  const fitChartData: IFitChartData = {
    chartOptions: fitChartOptions,
    seriesOptions: fitSeriesOptions,
    series: fitSeries
  };

  console.log(fitChartData);

  return fitChartData;
}
