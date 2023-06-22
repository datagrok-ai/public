import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';
import {createSigmoidPoints, createDemoDataFrame} from '../fit/fit-demo';
import {category, test, delay} from '@datagrok-libraries/utils/src/test';
import {FitChartCellRenderer} from '../fit/fit-renderer';
import {IFitChartData} from '@datagrok-libraries/statistics/src/fit/fit-curve';


export function createIFitChartData(seriesLength: number): IFitChartData {
	const step = 0.5;
	const chartData: IFitChartData = {
		series: [],
		chartOptions: {showStatistics: ['auc']},
	};

	const points = createSigmoidPoints(seriesLength, step);
	const color = DG.Color.toHtml(DG.Color.scatterPlotMarker);
	chartData.series![0] = {
		fitLineColor: color,
		pointColor: color,
		showFitLine: true,
		showPoints: 'points',
		showCurveConfidenceInterval: true,
		fitFunction: 'sigmoid',
		points: wu.count().take(seriesLength)
			.map(function(i) { return {x: points.x[i], y: points.y[i]}; })
			.toArray()
	};

	return chartData;
}

export async function scrollTable(el: HTMLElement, delta: number, cycles: number, secDelay: number) {
  for (let i = 0; i < cycles; i++) {
    el.dispatchEvent(new WheelEvent('wheel', {deltaY: delta}));
    await delay(secDelay);
  }
}


category('creation', () => {
	test('fit curve creation', async () => {
		const df = createDemoDataFrame(DG.Test.isInBenchmark ? 1000 : 100, 5, 2);
		grok.shell.addTableView(df);
		await delay(50);
	});
});

category('rendering', () => {
	// TODO: return amount of rendering things and same for different benchmarks, same for fitTests also ask Ann
	test('rendering same fitChartData on canvas', async () => {
		const canvasWidth = 160;
		const canvasHeight = 100;

		const canvas = ui.canvas(canvasWidth, canvasHeight);
		grok.shell.newView('canvas', [canvas]);
		const g = canvas.getContext('2d')!;

		const fitChartCellRenderer = new FitChartCellRenderer();
		const fitChartData: IFitChartData = createIFitChartData(DG.Test.isInBenchmark ? 50 : 15);
		for (let i = 1; i <= 1000; i++) {
			fitChartCellRenderer.renderCurves(g, canvas.clientLeft, canvas.clientTop, canvasWidth, canvasHeight, fitChartData);
			g.clearRect(0, 0, canvasWidth, canvasHeight);
		}
	});

	test('rendering different fitChartData on canvas', async () => {
		const canvasWidth = 160;
		const canvasHeight = 100;

		const canvas = ui.canvas(canvasWidth, canvasHeight);
		grok.shell.newView('canvas', [canvas]);
		const g = canvas.getContext('2d')!;

		const fitChartCellRenderer = new FitChartCellRenderer();
		for (let i = 1; i <= 1000; i++) {
			const fitChartData = createIFitChartData(DG.Test.isInBenchmark ? 50 : 15);
			fitChartCellRenderer.renderCurves(g, canvas.clientLeft, canvas.clientTop, canvasWidth, canvasHeight, fitChartData);
			g.clearRect(0, 0, canvasWidth, canvasHeight);
		}
	});

	test('rendering in grid', async () => {
		const df = createDemoDataFrame(30, 5, 2);
		const tv = grok.shell.addTableView(df);
		const scrollCycles = 10;
		const scrollDelta = 300;
		await delay(100);
		const canvas = tv.grid.root.getElementsByTagName('canvas')[2];
		const start = new Date();
		await scrollTable(canvas, scrollDelta, scrollCycles, 10);
		const stop = new Date();
		console.log(`Time for curves rendering is ${stop.valueOf() - start.valueOf()} ms`);
	});
});
