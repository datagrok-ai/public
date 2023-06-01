import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, test, delay} from '@datagrok-libraries/utils/src/test';
import {createDemoDataFrame} from '../fit/fit-demo';

export async function scrollTable(el: HTMLElement, delta: number, cycles: number, secDelay: number) {
  for (let i = 0; i < cycles; i++) {
    el.dispatchEvent(new WheelEvent('wheel', {deltaY: delta}));
    await delay(secDelay);
  }
}


category('creation', () => {
	test('fit curve creation', async () => {
		const df = DG.Test.isInBenchmark ? createDemoDataFrame(1000, 5, 2) : createDemoDataFrame(100, 5, 2);
		grok.shell.addTableView(df);
	});
});

category('rendering', () => {
	test('fit curve cell rendering', async () => {
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
