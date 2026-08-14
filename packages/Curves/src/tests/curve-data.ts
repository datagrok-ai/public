import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {IFitChartData} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {FitChartCellRenderer} from '../fit/fit-renderer';
import {isTitleShown, layoutChart, shownAxesNames} from '../fit/fit-layout';
import {getOrCreateParsedChartData} from '../fit/fit-chart-data';

/** Renders a cell onto a canvas and returns every string the chart drew, with where it put it and
 * the plot it was laid out in - which is what the legend hover is asked about. */
export function renderedRows(json: string, name: string, width: number = 400, height: number = 300):
  {rows: {text: string, x: number, y: number}[], data: IFitChartData, dataBox: DG.Rect} {
  const col = DG.Column.fromStrings('curve', [json]);
  col.semType = FitConstants.FIT_SEM_TYPE;
  const df = DG.DataFrame.fromColumns([col]);
  df.name = name;

  const rows: {text: string, x: number, y: number}[] = [];
  const g = ui.canvas(width, height).getContext('2d')!;
  const fillText = g.fillText.bind(g);
  g.fillText = ((text: string, x: number, y: number) => {
    rows.push({text: text, x: x, y: y});
    fillText(text, x, y);
  }) as any;

  const rect = new DG.Rect(0, 0, width, height);
  const data = getOrCreateParsedChartData(df.cell(0, 'curve'));
  new FitChartCellRenderer().renderCurves(g, rect, data);
  const axesNames = shownAxesNames(rect, data);
  const dataBox = layoutChart(rect, axesNames.x, axesNames.y, isTitleShown(rect, data))[0];
  return {rows: rows, data: data, dataBox: dataBox};
}

/** Every string the chart drew. */
export function renderedTexts(json: string, name: string, width: number = 400, height: number = 300): string[] {
  return renderedRows(json, name, width, height).rows.map((row) => row.text);
}

export const CONCENTRATIONS = [1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4];

/** Sigmoid dose-response points with the inflection point at 10^`logIC50`. */
export function sigmoidPoints(logIC50: number): {x: number, y: number}[] {
  return CONCENTRATIONS.map((x) => ({x: x, y: 5 + 95 / (1 + Math.pow(10, Math.log10(x) - logIC50))}));
}

/** One sigmoid series in a cell, with logX baked into the cell's chart options. */
export function curveJson(logIC50: number): string {
  return JSON.stringify({
    chartOptions: {logX: true},
    series: [{fitFunction: 'sigmoid', name: 'series', points: sigmoidPoints(logIC50)}],
  } as IFitChartData);
}

/** One sigmoid series that declares a confidence interval, the way the demo data does. */
export function ciCurveJson(): string {
  return JSON.stringify({
    chartOptions: {logX: true},
    series: [{fitFunction: 'sigmoid', name: 'series', showCurveConfidenceInterval: true, points: sigmoidPoints(-6.5)}],
  } as IFitChartData);
}

/** One cell holding several sigmoid series, each with its own inflection point. */
export function multiSeriesCurveJson(logIC50s: number[]): string {
  return JSON.stringify({
    chartOptions: {logX: true},
    series: logIC50s.map((logIC50, i) => ({
      fitFunction: 'sigmoid', name: `series ${i}`, points: sigmoidPoints(logIC50),
    })),
  } as IFitChartData);
}

/** Two sigmoid series, one plot-level label and one per curve. */
export function labelledCurveJson(): string {
  return JSON.stringify({
    chartOptions: {logX: true, labels: {'Z prime': 0.72}, showLabels: ['Z prime', 'compound']},
    series: [
      {fitFunction: 'sigmoid', name: 's0', labels: {compound: 'GRK-1'}, points: sigmoidPoints(-6.5)},
      {fitFunction: 'sigmoid', name: 's1', labels: {compound: 'GRK-2'}, points: sigmoidPoints(-5.5)},
    ],
  } as IFitChartData);
}
