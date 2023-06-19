import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {createDemoDataFrame} from '../fit/fit-demo';
import {category, test, delay} from '@datagrok-libraries/utils/src/test';
import {FitChartCellRenderer} from '../fit/fit-renderer';
import {IFitChartData} from '@datagrok-libraries/statistics/src/fit/fit-curve';


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
	test('rendering on canvas', async () => {
		const canvasWidth = 160;
		const canvasHeight = 100;

		const canvas = ui.canvas(canvasWidth, canvasHeight);
		grok.shell.newView('canvas', [canvas]);
		const g = canvas.getContext('2d')!;

		const fitChartCellRenderer = new FitChartCellRenderer();
		const fitChartData: IFitChartData = DG.Test.isInBenchmark ? {"series":[{"fitLineColor":"#2ca02c","pointColor":"#2ca02c","showCurveConfidenceInterval":true,"fitFunction": "sigmoid", "clickToToggle": false, "showFitLine": true, "showPoints": true, "showBoxPlot": true, "points":[{"x":0.10000000149011612,"y":1.8080450296401978},{"x":0.30000001192092896,"y":1.7073274850845337},{"x":0.5,"y":1.8684450387954712},{"x":0.699999988079071,"y":1.8083453178405762},{"x":0.8999999761581421,"y":1.871852993965149},{"x":1.100000023841858,"y":2.0314087867736816},{"x":1.2999999523162842,"y":1.8348426818847656},{"x":1.5,"y":1.716798186302185},{"x":1.7000000476837158,"y":1.9205678701400757},{"x":1.899999976158142,"y":2.040295124053955},{"x":2.0999999046325684,"y":1.915327787399292},{"x":2.299999952316284,"y":1.7997850179672241},{"x":2.5,"y":1.8200938701629639},{"x":2.700000047683716,"y":1.8472325801849365},{"x":2.9000000953674316,"y":1.892165184020996},{"x":3.0999999046325684,"y":1.6691789627075195},{"x":3.299999952316284,"y":1.7463487386703491},{"x":3.5,"y":1.856266975402832},{"x":3.700000047683716,"y":1.8447580337524414},{"x":3.9000000953674316,"y":2.044322967529297},{"x":4.099999904632568,"y":1.946313500404358},{"x":4.300000190734863,"y":1.7538775205612183},{"x":4.5,"y":1.7557778358459473},{"x":4.699999809265137,"y":1.743997573852539},{"x":4.900000095367432,"y":1.8651076555252075},{"x":5.099999904632568,"y":1.9675723314285278},{"x":5.300000190734863,"y":1.9454524517059326},{"x":5.5,"y":1.8253188133239746},{"x":5.699999809265137,"y":2.0103540420532227},{"x":5.900000095367432,"y":2.0355498790740967},{"x":6.099999904632568,"y":1.9354662895202637},{"x":6.300000190734863,"y":1.7890597581863403},{"x":6.5,"y":2.1066558361053467},{"x":6.699999809265137,"y":2.1674845218658447},{"x":6.900000095367432,"y":1.8644132614135742},{"x":7.099999904632568,"y":2.1932010650634766},{"x":7.300000190734863,"y":1.9899929761886597},{"x":7.5,"y":2.0531647205352783},{"x":7.699999809265137,"y":2.2773845195770264},{"x":7.900000095367432,"y":1.98076331615448},{"x":8.100000381469727,"y":2.2184560298919678},{"x":8.300000190734863,"y":2.323612689971924},{"x":8.5,"y":2.2961959838867188},{"x":8.699999809265137,"y":2.2926039695739746},{"x":8.899999618530273,"y":2.043086051940918},{"x":9.100000381469727,"y":2.1240429878234863},{"x":9.300000190734863,"y":1.9789459705352783},{"x":9.5,"y":2.192298173904419},{"x":9.699999809265137,"y":2.0946543216705322},{"x":0,"y":0.11970949918031693}]}],"chartOptions":{"showStatistics":["auc"]}} : {"series":[{"fitLineColor":"#2ca02c","pointColor":"#2ca02c","showCurveConfidenceInterval":true,"fitFunction":"sigmoid","parameters":[1.4638188806702541, -1.9110075693069182, 4.769806319123481, 0.10434338956506588],"points":[{"x":0.10000000149011612,"y":0.02634293958544731},{"x":0.6000000238418579,"y":0.08862859010696411},{"x":1.100000023841858,"y":0.12354502826929092},{"x":1.600000023841858,"y":0.21434032917022705},{"x":2.0999999046325684,"y":0.07175814360380173},{"x":2.5999999046325684,"y":0.11858399212360382},{"x":3.0999999046325684,"y":0.10206034034490585},{"x":3.5999999046325684,"y":0.1560754030942917},{"x":4.099999904632568,"y":0.09271059930324554},{"x":4.599999904632568,"y":0.5750495195388794},{"x":5.099999904632568,"y":1.175055742263794},{"x":5.599999904632568,"y":1.4829083681106567},{"x":6.099999904632568,"y":1.3193761110305786},{"x":6.599999904632568,"y":1.5237388610839844},{"x":7.099999904632568,"y":1.50946044921875}]}],"chartOptions":{"showStatistics":["auc"]}};
		fitChartCellRenderer.renderCurves(g, canvas.clientLeft, canvas.clientTop, canvasWidth, canvasHeight, fitChartData);
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
