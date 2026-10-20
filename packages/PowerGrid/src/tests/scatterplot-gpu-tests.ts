import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {awaitCheck, before, category, delay, expect, test, TestOptions} from '@datagrok-libraries/test/src/test';
import {getGPUDevice} from '@datagrok-libraries/math/src/webGPU/getGPUDevice';
import {isWebGLAvailable, scWebGLPointHitTest, scWebGLRender} from '../webgl/scatterplot';
import {scWebGPUPointHitTest, scWebGPURender} from '../webgpu/scatterplot';
import {isGpuRenderValid} from '../scatterplot-gpu';

const FLOAT_NONE = Math.fround(DG.FLOAT_NULL);

/** A GPU renderer as the tests drive it, directly: the test bundle has its own renderers, beside the package's. */
interface Backend {
  name: string;
  options: TestOptions;
  init(): Promise<void>;
  render(sp: DG.ScatterPlotViewer): HTMLCanvasElement;
  hitTest(sp: DG.ScatterPlotViewer, pt: DG.Point): Promise<number | null>;
  /** Errors the backend reported since the last call. */
  errors(): Promise<string[]>;
}

let device: GPUDevice | null = null;
const deviceErrors: string[] = [];
const backends: Backend[] = [{
  name: 'WebGL',
  options: {skipReason: isWebGLAvailable() ? undefined : 'No hardware WebGL 2 in this browser'},
  init: async () => {},
  render: (sp) => scWebGLRender(sp, true)!,
  hitTest: async (sp, pt) => scWebGLPointHitTest(sp, pt),
  errors: async () => [],
}, {
  name: 'WebGPU',
  options: {},
  init: async function(this: Backend) {
    device = await getGPUDevice();
    if (!device)
      this.options.skipReason = 'No WebGPU device in this browser';
    device?.addEventListener('uncapturederror', (e) => deviceErrors.push((e as GPUUncapturedErrorEvent).error.message));
  },
  render: (sp) => scWebGPURender(sp, true, device)!,
  hitTest: async (sp, pt) => scWebGPUPointHitTest(sp, pt),
  errors: async () => {
    await device!.queue.onSubmittedWorkDone();
    await delay(50);
    return deviceErrors.splice(0);
  },
}];

/** A scatter plot of a demog table in its own view; the test runner closes all views after each test. */
async function scatterPlot(rows: number, options: object = {}): Promise<DG.ScatterPlotViewer> {
  const view = grok.shell.addTableView(grok.data.demo.demog(rows));
  const sp = view.scatterPlot({x: 'height', y: 'weight', markerType: DG.MARKER_TYPE.DOT, ...options});
  await rendered(sp);
  return sp;
}

async function rendered(sp: DG.ScatterPlotViewer): Promise<void> {
  await delay(100);
  await awaitCheck(() => sp.viewBox.width > 100 && !sp.isRenderPending, 'the scatter plot did not render', 5000);
}

/** Renders [sp] and reads the canvas back in the same task, before the drawing is presented. */
function renderPixels(backend: Backend, sp: DG.ScatterPlotViewer):
  {pixels: Uint8ClampedArray, left: number, top: number, width: number} {
  const canvas = backend.render(sp);
  const copy = new OffscreenCanvas(canvas.width, canvas.height);
  const g = copy.getContext('2d')!;
  g.drawImage(canvas, 0, 0);
  const dpr = window.devicePixelRatio;
  return {pixels: g.getImageData(0, 0, canvas.width, canvas.height).data, width: canvas.width,
    left: Math.round(parseFloat(canvas.style.left) * dpr), top: Math.round(parseFloat(canvas.style.top) * dpr)};
}

/** Where the core's dot renderer puts each filtered row: one device pixel at the floored screen position. */
function dotRendererPixels(sp: DG.ScatterPlotViewer): Set<number> {
  const df = sp.dataFrame;
  const xs = df.col(sp.props.xColumnName)!.getRawData();
  const ys = df.col(sp.props.yColumnName)!.getRawData();
  const vb = sp.viewBox;
  const vp = sp.viewport;
  const dpr = window.devicePixelRatio;
  const pixels = new Set<number>();
  for (const i of sp.filter.getSelectedIndexes()) {
    if (xs[i] == FLOAT_NONE || ys[i] == FLOAT_NONE)
      continue;
    const x = vb.x + (xs[i] - vp.x) / vp.width * vb.width;
    const y = vb.y + (1 - (ys[i] - vp.y) / vp.height) * vb.height;
    if (x >= vb.x && x <= vb.x + vb.width && y >= vb.y && y <= vb.y + vb.height)
      pixels.add(Math.floor(y * dpr) * 65536 + Math.floor(x * dpr));
  }
  return pixels;
}

function expectDotRendererPixels(backend: Backend, sp: DG.ScatterPlotViewer): void {
  const {pixels, left, top, width} = renderPixels(backend, sp);
  const drawn = new Set<number>();
  for (let i = 0; i < pixels.length / 4; i++) {
    if (pixels[i * 4 + 3] > 0)
      drawn.add((top + Math.floor(i / width)) * 65536 + left + i % width);
  }
  const expected = dotRendererPixels(sp);
  let missing = 0;
  for (const p of expected) {
    if (!drawn.has(p))
      missing++;
  }
  // float32 positions may round a few rows on a pixel boundary the other way
  expect(missing <= expected.size * 0.001 && Math.abs(drawn.size - expected.size) <= expected.size * 0.001, true,
    `${drawn.size} pixels drawn, ${expected.size} expected, ${missing} of them missing`);
}

for (const backend of backends) {
  category(`Scatter plot: ${backend.name}`, () => {
    before(() => backend.init());
    // the package gives a backend up on its first error, even one that changes no pixel
    const gpuTest = (name: string, body: () => Promise<void>) => test(name, async () => {
      await body();
      expect((await backend.errors()).join('; '), '');
    }, backend.options);

    gpuTest('Dots land where the dot renderer puts them', async () => {
      expectDotRendererPixels(backend, await scatterPlot(200000));
    });

    gpuTest('Follows a change of the row source', async () => {
      // each row source has its own filter bitset, and bitset versions start over: only identity tells them apart
      const sp = await scatterPlot(200000);
      sp.dataFrame.filter.init((i) => i % 2 == 0);
      await rendered(sp);
      expectDotRendererPixels(backend, sp);
      sp.props.rowSource = 'All';
      await rendered(sp);
      expectDotRendererPixels(backend, sp);
    });

    gpuTest('Dots take the row colors', async () => {
      const sp = await scatterPlot(200000, {colorColumnName: 'race'});
      const rowColors = new Set(Array.from(sp.getMarkerColors(), (c) => c & 0xFFFFFF));
      const {pixels} = renderPixels(backend, sp);
      const drawn = new Set<number>();
      for (let i = 0; i < pixels.length; i += 4) {
        if (pixels[i + 3] > 0)
          drawn.add((pixels[i] << 16) | (pixels[i + 1] << 8) | pixels[i + 2]);
      }
      expect([...drawn].every((c) => rowColors.has(c)), true, 'a dot is drawn in a color no row has');
      expect(drawn.size, sp.dataFrame.col('race')!.categories.length);
    });

    gpuTest('Hit test matches the CPU', async () => {
      const sp = await scatterPlot(200000);
      for (const props of [{markerType: DG.MARKER_TYPE.DOT}, {markerType: DG.MARKER_TYPE.CIRCLE, sizeColumnName: 'age'},
        {markerType: DG.MARKER_TYPE.TRIANGLE_TOP, sizeColumnName: '', markerDefaultSize: 12}]) {
        Object.assign(sp.props, props);
        // the core's hit test takes its search radius from the layout of the next frame
        await rendered(sp);
        backend.render(sp);
        const rows = sp.filter.getSelectedIndexes();
        const vb = sp.viewBox;
        for (let k = 0; k < 100; k++) {
          const p = sp.pointToScreen(rows[(k * 7919) % rows.length]);
          // off the marker center by a fraction: a point exactly on a hit radius is a float32 vs float64 coin toss
          const near = [p.x + (k % 5) - 1.75, p.y + (k % 3) - 0.75];
          for (const [x, y] of [near, [vb.x + (k * 37) % vb.width, vb.y + (k * 53) % vb.height]]) {
            const cpu = sp.hitTest(x, y);
            expect((await backend.hitTest(sp, new DG.Point(x, y))) ?? cpu, cpu, `${props.markerType} at (${x}, ${y})`);
          }
        }
      }
    });

    gpuTest('Releases the canvas when the viewer closes', async () => {
      const sp = await scatterPlot(60000);
      const canvas = backend.render(sp);
      expect(canvas.isConnected, true);
      grok.shell.closeAll();
      await awaitCheck(() => !canvas.isConnected, 'the canvas outlived its viewer', 3000);
    });
  });
}

category('Scatter plot: GPU', () => {
  test('Leaves jittered, packed and category-ordered plots to the CPU', async () => {
    const sp = await scatterPlot(60000, {jitterSize: 5});
    expect(isGpuRenderValid(sp), false);
    Object.assign(sp.props, {jitterSize: 0, zoomAndFilter: 'pack and zoom by filter'});
    expect(isGpuRenderValid(sp), false);
    Object.assign(sp.props, {zoomAndFilter: 'filter by zoom', xColumnName: 'race',
      xCategoryOrder: '{"manual":["Other"]}'});
    expect(isGpuRenderValid(sp), false);
    // an order left on an axis turned numeric does not apply
    sp.props.xColumnName = 'height';
    expect(isGpuRenderValid(sp), true);
  });
});
