import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {getGPUDevice, isIntegratedGPU} from '@datagrok-libraries/math/src/webGPU/getGPUDevice';

category('webgpu: integrated GPU', () => {
  test('isIntegratedGPU', async () => {
    const d3d = (gpu: string) => `ANGLE (${gpu} (0x00001114) Direct3D11 vs_5_0 ps_5_0, D3D11)`;
    for (const gpu of ['AMD, AMD Radeon(TM) 860M Graphics', 'AMD, AMD Radeon(TM) Graphics',
      'Intel, Intel(R) UHD Graphics 770', 'Intel, Intel(R) Arc(TM) 140V GPU (16GB)'])
      expect(isIntegratedGPU(d3d(gpu)), true, gpu);
    for (const gpu of ['NVIDIA, NVIDIA GeForce RTX 5070 Laptop GPU', 'AMD, AMD Radeon RX 7600M XT',
      'AMD, AMD Radeon PRO W7900', 'Intel, Intel(R) Arc(TM) A770 Graphics',
      'Google, Vulkan 1.3.0 (SwiftShader Device)'])
      expect(isIntegratedGPU(d3d(gpu)), false, gpu);
  });
});

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
