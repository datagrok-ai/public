import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test} from '@datagrok-libraries/test/src/test';
import {getGPUDevice} from '@datagrok-libraries/math/src/webGPU/getGPUDevice';

const sizes = [10000, 1000000, 10000000];
const types = ['dot', 'circle'];
const consecutiveLaunches = 3;

// registers its tests synchronously: tests a category body registers after an await land in a later category
category('webgpu', () => {
  const render = DG.Func.byName('_scWebGPURender');
  const hitTest = DG.Func.byName('_scWebGPUPointHitTest');

  types.forEach((t) => {
    sizes.forEach((s) => {
      test(`GPURender: ${s / 1000000.0}m ${t}s`, async () => {
        if (!DG.Test.isInBenchmark || !(await getGPUDevice()))
          return;
        const view = grok.shell.addTableView(grok.data.demo.randomWalk(s, 2));
        const sp: DG.ScatterPlotViewer = view.scatterPlot();
        sp.props.markerType = t;
        for (let i = 1; i <= consecutiveLaunches; ++i) {
          await DG.timeAsync(`Render Call #${i}`, async () => {
            await render.apply([sp, true]);
          });
        }
      }, {benchmarkTimeout: 10000, benchmark: true, skipReason: 'Run manually in "Run Skipped" mode'});

      test(`GPURender Size: ${s / 1000000.0}m ${t}s`, async () => {
        if (!DG.Test.isInBenchmark || !(await getGPUDevice()))
          return;
        const view = grok.shell.addTableView(grok.data.demo.randomWalk(s, 2));
        const sp: DG.ScatterPlotViewer = view.scatterPlot();
        sp.props.markerType = t;
        sp.props.sizeColumnName = '#1';
        for (let i = 1; i <= consecutiveLaunches; ++i) {
          await DG.timeAsync(`Render Call #${i}`, async () => {
            await render.apply([sp, true]);
          });
        }
      }, {benchmarkTimeout: 10000, benchmark: true, skipReason: 'Run manually in "Run Skipped" mode'});

      test(`GPURender Color: ${s / 1000000.0}m ${t}s`, async () => {
        if (!DG.Test.isInBenchmark || !(await getGPUDevice()))
          return;
        const view = grok.shell.addTableView(grok.data.demo.randomWalk(s, 2));
        const sp: DG.ScatterPlotViewer = view.scatterPlot();
        sp.props.markerType = t;
        sp.props.colorColumnName = '#0';
        for (let i = 1; i <= consecutiveLaunches; ++i) {
          await DG.timeAsync(`Render Call #${i}`, async () => {
            await render.apply([sp, true]);
          });
        }
      }, {benchmarkTimeout: 10000, benchmark: true, skipReason: 'Run manually in "Run Skipped" mode'});

      test(`GPUHitTest: ${s / 1000000.0}m ${t}s`, async () => {
        if (!DG.Test.isInBenchmark || !(await getGPUDevice()))
          return;
        const view = grok.shell.addTableView(grok.data.demo.randomWalk(s, 2));
        const sp: DG.ScatterPlotViewer = view.scatterPlot();
        sp.props.markerType = t;
        for (let i = 1; i <= consecutiveLaunches; ++i) {
          await DG.timeAsync(`GPU HitTest call #${i}`, async () => {
            await hitTest.apply([sp, new DG.Point(0, 0)]);
          });
        }
      }, {benchmarkTimeout: 10000, benchmark: true, skipReason: 'Run manually in "Run Skipped" mode'});
    });
  });
}, {benchmarks: true});
