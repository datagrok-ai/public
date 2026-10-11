/* Chooses how the core's large scatter plots are drawn: WebGPU where the browser allows it, WebGL 2 where it does not,
 * and the core's own CPU renderers when neither works or the plot needs what they lack. The core reaches this through
 * the `isWebGPUAvailable`, `isWebGPURenderValid`, `_scWebGPURender` and `_scWebGPUPointHitTest` package functions,
 * named before the WebGL backend existed. */
import * as DG from 'datagrok-api/dg';
import {scWebGPUPointHitTest, scWebGPURender} from './webgpu/scatterplot';
import {disposeWebGLRenderers, isWebGLAvailable, isWebGLRenderValid, scWebGLPointHitTest,
  scWebGLRender} from './webgl/scatterplot';
import {isGpuRenderable, readLook} from './scatterplot-look';

type Backend = 'WebGPU' | 'WebGL';

let webGPUDevice = false;
const failed = new Set<Backend>();

export function setGPUDevice(device: GPUDevice | null): void {
  webGPUDevice = device != null;
  // the WebGPU renderer keeps its buffers on this device, and a lost device fails silently
  device?.lost.then(() => failed.add('WebGPU'));
}

function backend(): Backend | null {
  if (webGPUDevice && !failed.has('WebGPU'))
    return 'WebGPU';
  if (!failed.has('WebGL') && isWebGLAvailable())
    return 'WebGL';
  return null;
}

/** A backend that fails once is not used again; the core redraws the plot with the next one, or on the CPU. */
function fail(failedBackend: Backend, sc: DG.ScatterPlotViewer, error: unknown): void {
  failed.add(failedBackend);
  console.warn(`Scatter plot: ${failedBackend} rendering failed, falling back`, error);
  if (failedBackend == 'WebGL')
    disposeWebGLRenderers();
  sc.invalidateCanvas();
}

export function isGpuAvailable(): boolean {
  return backend() != null;
}

export function isGpuRenderValid(sc: DG.ScatterPlotViewer): boolean {
  const current = backend();
  const look = readLook(sc);
  if (!current || !isGpuRenderable(sc, look))
    return false;
  // per-row marker shapes are WebGL only
  return current == 'WebGPU' ? !look.markersColumnName : isWebGLRenderValid(sc, look);
}

export async function gpuRender(sc: DG.ScatterPlotViewer, show: boolean): Promise<void> {
  const current = show ? backend() : null;
  if (current != 'WebGL')
    scWebGLRender(sc, false);
  // hiding is synchronous: awaiting it would let a later call for this plot run in between
  if (current != 'WebGPU')
    void scWebGPURender(sc, false);
  if (!current)
    return;
  try {
    if (current == 'WebGL')
      scWebGLRender(sc, true);
    else
      await scWebGPURender(sc, true);
  }
  catch (error) {
    fail(current, sc, error);
  }
}

export async function gpuPointHitTest(sc: DG.ScatterPlotViewer, pt: DG.Point): Promise<number> {
  const current = backend();
  try {
    if (current == 'WebGPU')
      return await scWebGPUPointHitTest(sc, pt);
    if (current == 'WebGL')
      return scWebGLPointHitTest(sc, pt) ?? sc.hitTest(pt.x, pt.y);
  }
  catch (error) {
    fail(current!, sc, error);
  }
  return sc.hitTest(pt.x, pt.y);
}
