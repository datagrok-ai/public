import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {awaitCheck, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {isWebGLAvailable, scWebGLPointHitTest, scWebGLRender} from '../webgl/scatterplot';
import {isGpuRenderValid} from '../scatterplot-gpu';

const FLOAT_NONE = Math.fround(DG.FLOAT_NULL);
const webGL = {skipReason: isWebGLAvailable() ? undefined : 'No hardware WebGL 2 in this browser'};

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

/** Renders [sp] with WebGL and reads the canvas back in the same task, before the drawing buffer is presented. */
function renderPixels(sp: DG.ScatterPlotViewer): {pixels: Uint8ClampedArray, left: number, top: number, width: number} {
  // the test bundle has its own renderer, beside the package's: read the canvas it drew on
  const gl = scWebGLRender(sp, true)!;
  const copy = new OffscreenCanvas(gl.width, gl.height);
  const g = copy.getContext('2d')!;
  g.drawImage(gl, 0, 0);
  const dpr = window.devicePixelRatio;
  return {pixels: g.getImageData(0, 0, gl.width, gl.height).data, width: gl.width,
    left: Math.round(parseFloat(gl.style.left) * dpr), top: Math.round(parseFloat(gl.style.top) * dpr)};
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

function expectDotRendererPixels(sp: DG.ScatterPlotViewer): void {
  const {pixels, left, top, width} = renderPixels(sp);
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

category('Scatter plot: WebGL', () => {
  test('Dots land where the dot renderer puts them', async () => {
    expectDotRendererPixels(await scatterPlot(200000));
  }, webGL);

  test('Follows a change of the row source', async () => {
    // each row source has its own filter bitset, and bitset versions start over: only identity tells them apart
    const sp = await scatterPlot(200000);
    sp.dataFrame.filter.init((i) => i % 2 == 0);
    await rendered(sp);
    expectDotRendererPixels(sp);
    sp.props.rowSource = 'All';
    await rendered(sp);
    expectDotRendererPixels(sp);
  }, webGL);

  test('Dots take the row colors', async () => {
    const sp = await scatterPlot(200000, {colorColumnName: 'race'});
    const rowColors = new Set(Array.from(sp.getMarkerColors(), (c) => c & 0xFFFFFF));
    const {pixels} = renderPixels(sp);
    const drawn = new Set<number>();
    for (let i = 0; i < pixels.length; i += 4) {
      if (pixels[i + 3] > 0)
        drawn.add((pixels[i] << 16) | (pixels[i + 1] << 8) | pixels[i + 2]);
    }
    expect([...drawn].every((c) => rowColors.has(c)), true, 'a dot is drawn in a color no row has');
    expect(drawn.size, sp.dataFrame.col('race')!.categories.length);
  }, webGL);

  test('Hit test matches the CPU', async () => {
    const sp = await scatterPlot(200000);
    for (const props of [{markerType: DG.MARKER_TYPE.DOT}, {markerType: DG.MARKER_TYPE.CIRCLE, sizeColumnName: 'age'},
      {markerType: DG.MARKER_TYPE.TRIANGLE_TOP, sizeColumnName: '', markerDefaultSize: 12}]) {
      Object.assign(sp.props, props);
      // the core's hit test takes its search radius from the layout of the next frame
      await delay(100);
      await awaitCheck(() => !sp.isRenderPending, 'the scatter plot did not render', 5000);
      scWebGLRender(sp, true);
      const rows = sp.filter.getSelectedIndexes();
      const vb = sp.viewBox;
      for (let k = 0; k < 100; k++) {
        const p = sp.pointToScreen(rows[(k * 7919) % rows.length]);
        // off the marker center by a fraction: a point exactly on a hit radius is a float32 vs float64 coin toss
        const near = [p.x + (k % 5) - 1.75, p.y + (k % 3) - 0.75];
        for (const [x, y] of [near, [vb.x + (k * 37) % vb.width, vb.y + (k * 53) % vb.height]]) {
          const cpu = sp.hitTest(x, y);
          expect(scWebGLPointHitTest(sp, new DG.Point(x, y)) ?? cpu, cpu, `${props.markerType} at (${x}, ${y})`);
        }
      }
    }
  }, webGL);

  test('Leaves jittered and packed plots to the CPU', async () => {
    const sp = await scatterPlot(60000, {jitterSize: 5});
    expect(isGpuRenderValid(sp), false);
    Object.assign(sp.props, {jitterSize: 0, zoomAndFilter: 'pack and zoom by filter'});
    expect(isGpuRenderValid(sp), false);
  });

  test('Releases the canvas when the viewer closes', async () => {
    const sp = await scatterPlot(60000);
    const canvas = scWebGLRender(sp, true)!;
    expect(canvas.isConnected, true);
    grok.shell.closeAll();
    await awaitCheck(() => !canvas.isConnected, 'the WebGL canvas outlived its viewer', 3000);
  }, webGL);
});
