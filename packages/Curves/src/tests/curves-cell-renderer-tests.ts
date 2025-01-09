import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';
import {createSigmoidPoints, createDemoDataFrame} from '../fit/fit-demo';
import {category, test, delay} from '@datagrok-libraries/utils/src/test';
import {FitChartCellRenderer} from '../fit/fit-renderer';
import {IFitChartData} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {FitConstants} from '../fit/const';


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
		markerType: DG.MARKER_TYPE.CIRCLE,
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

async function getTime(func: (chartData?: IFitChartData | null) => Promise<void>, funcChartData?: IFitChartData | null): Promise<number> {
	const start = Date.now();
	await func(funcChartData);
	const stop = Date.now();
	return stop.valueOf() - start.valueOf();
}


category('creation', () => {
	test('fit curve creation', async () => {
		const pointsAmount = DG.Test.isInBenchmark ? 1000 : 100;
		const func = async (chartData?: IFitChartData | null) => {
			const df = createDemoDataFrame(pointsAmount, 5, 2);
			grok.shell.addTableView(df);
		};
		const time = await getTime(func, null);
		console.log(`Creation took ${time} ms`)
		return `${pointsAmount} points, took ${time} ms`
	}, {benchmark: true});
});

category('rendering', () => {
	test('rendering same fitChartData on canvas', async () => {
		const canvasWidth = 160;
		const canvasHeight = 120;
		const renderingTimesAmount = 1000;

		const canvas = ui.canvas(canvasWidth, canvasHeight);
		grok.shell.newView('canvas', [canvas]);
		const g = canvas.getContext('2d')!;

		const fitChartCellRenderer = new FitChartCellRenderer();
		const fitChartData: IFitChartData = createIFitChartData(DG.Test.isInBenchmark ? 50 : 15);
		const times: number[] = [];
		const func = async (chartData?: IFitChartData | null) => {
			fitChartCellRenderer.renderCurves(g, FitChartCellRenderer.inflateScreenBounds(new DG.Rect(canvas.clientLeft, canvas.clientTop, canvasWidth, canvasHeight)), fitChartData);
			g.clearRect(0, 0, canvasWidth, canvasHeight);
		};
		for (let i = 1; i <= renderingTimesAmount; i++)
			times[times.length] = await getTime(func);
		return `rendering performed ${renderingTimesAmount} times, max time is ${Math.max(...times)} ms, min time is ${Math.min(...times)} ms, average time is ${times.reduce((a, b) => a + b, 0) / times.length} ms`;
	}, {benchmark: true});

	test('rendering different fitChartData on canvas', async () => {
		const canvasWidth = 160;
		const canvasHeight = 120;
		const renderingTimesAmount = 1000;

		const canvas = ui.canvas(canvasWidth, canvasHeight);
		grok.shell.newView('canvas', [canvas]);
		const g = canvas.getContext('2d')!;

		const fitChartCellRenderer = new FitChartCellRenderer();
		const times: number[] = [];
		const func = async (chartData?: IFitChartData | null) => {
			fitChartCellRenderer.renderCurves(g, FitChartCellRenderer.inflateScreenBounds(new DG.Rect(canvas.clientLeft, canvas.clientTop, canvasWidth, canvasHeight)), chartData!);
			g.clearRect(0, 0, canvasWidth, canvasHeight);
		}
		for (let i = 1; i <= renderingTimesAmount; i++) {
			const fitChartData = createIFitChartData(DG.Test.isInBenchmark ? 50 : 15);
			times[times.length] = await getTime(func, fitChartData);
		}
		return `rendering performed ${renderingTimesAmount} times, max time is ${Math.max(...times)} ms, min time is ${Math.min(...times)} ms, average time is ${times.reduce((a, b) => a + b, 0) / times.length} ms`;
	}, {benchmark: true});

	test('rendering in grid', async () => {
		const df = createDemoDataFrame(50, 5, 5);
		const tv = grok.shell.addTableView(df);
		const scrollCycles = 10;
		const scrollDeltaPlus = 300;
		const scrollDeltaMinus = -300;
		await delay(100);
		const canvas = tv.grid.root.getElementsByTagName('canvas')[2];
		const func = async (chartData?: IFitChartData | null) => {
			await scrollTable(canvas, scrollDeltaPlus, scrollCycles, 1);
			await scrollTable(canvas, scrollDeltaMinus, scrollCycles, 1);
			await scrollTable(canvas, scrollDeltaPlus, scrollCycles, 1);
		};
		const time = await getTime(func, null) - 30;
		console.log(`curves rendering took ${time} ms`);
		return time;
	}, {benchmark: true});
});
