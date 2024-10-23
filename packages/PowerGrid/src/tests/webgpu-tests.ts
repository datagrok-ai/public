import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {awaitCheck, category, expect, expectArray, test, testEvent} from '@datagrok-libraries/utils/src/test';
import { _scWebGPURender, _scWebGPUPointHitTest } from '../package';


category('webgpu', () => {
  test('gpuRender', async () => {
    const types = ['dot', 'circle'];
    for (let size = 10000000; size <= 10000000; size *= 10) {
    for (let type = 0; type < types.length; ++type) {
        const df: DG.DataFrame = grok.data.demo.randomWalk(size, 2);
        df.name = 'test';
        const view = grok.shell.addTableView(df);
        df.currentRowIdx = 0;
    
        const sp: DG.ScatterPlotViewer = view.scatterPlot();
        sp.props.markerType = types[type];
  
        for (let i = 1; i <= 3; ++i) {
          await DG.timeAsync(`${i} GPU Render ${size / 1000000.0}m ${types[type]}s`, async() => {
            await _scWebGPURender(sp, true);
          });
        }
      }
    }
  });
  
  test('gpuHitTest', async () => {
    const types = ['dot', 'circle'];
    for (let size = 10000000; size <= 10000000; size *= 10) {
    for (let type = 0; type < types.length; ++type) {
        const df: DG.DataFrame = grok.data.demo.randomWalk(size, 2);
        df.name = 'test';
        const view = grok.shell.addTableView(df);
        df.currentRowIdx = 0;
    
        const sp: DG.ScatterPlotViewer = view.scatterPlot();
        sp.props.markerType = types[type];
  
        for (let i = 1; i <= 3; ++i) {
          await DG.timeAsync(`${i} GPU Hit Test ${size / 1000000.0}m ${types[type]}s`, async() => {
            await _scWebGPUPointHitTest(sp, new DG.Point(0, 0));
          });
        }
      }
    }
  });
});
