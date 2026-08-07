import {IFitChartData} from '@datagrok-libraries/statistics/src/fit/fit-curve';

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
