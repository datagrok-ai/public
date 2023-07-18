import * as DG from 'datagrok-api/dg';

import {Viewport} from '@datagrok-libraries/utils/src/transform';
import {layoutChart} from '../fit/fit-renderer';

import {category, expect, test} from '@datagrok-libraries/utils/src/test';


category('viewport', () => {
  test('viewportMethods', async () => {
    const screenBounds = new DG.Rect(120, 20, 160, 100).inflate(-6, -6);
    const [dataBox] = layoutChart(screenBounds);
    const dataBounds = new DG.Rect(0.1, 0.06, 7, 5.94);
    const viewport = new Viewport(dataBounds, dataBox, false, false);

		expect(viewport.xToScreen(1.5), 179.6);
		expect(viewport.yToScreen(2.33), 73.01346801346801);
		expect(viewport.xToWorld(150), -0.2559322033898306);
		expect(viewport.yToWorld(67), 2.855294117647059);
  });
});
