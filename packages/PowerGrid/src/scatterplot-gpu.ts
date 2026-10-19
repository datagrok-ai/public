/* Chooses how the core's large scatter plots are drawn: WebGPU where the browser allows it, WebGL 2 where it does not,
 * and the core's own CPU renderers when neither works or the plot needs what they lack. The core reaches this through
 * the `isWebGPUAvailable`, `isWebGPURenderValid`, `_scWebGPURender` and `_scWebGPUPointHitTest` package functions,
 * named before the WebGL backend existed. */
import * as DG from 'datagrok-api/dg';
import {disposeWebGPURenderers, scWebGPUPointHitTest, scWebGPURender} from './webgpu/scatterplot';
import {disposeWebGLRenderers, isWebGLAvailable, isWebGLRenderValid, scWebGLPointHitTest,
  scWebGLRender} from './webgl/scatterplot';
import {isGpuRenderable, readCheckLook} from './scatterplot-look';

type Backend = 'WebGPU' | 'WebGL';

let gpuDevice: GPUDevice | null = null;
const failed = new Set<Backend>();

export function setGPUDevice(device: GPUDevice | null): void {
  gpuDevice = device;
  // WebGPU reports errors as events, not exceptions
  device?.addEventListener('uncapturederror', (e) => fail('WebGPU', (e as GPUUncapturedErrorEvent).error));
  device?.lost.then((info) => fail('WebGPU', info.message));
}

function backend(): Backend | null {
  if (gpuDevice && !failed.has('WebGPU'))
    return 'WebGPU';
  if (!failed.has('WebGL') && isWebGLAvailable())
    return 'WebGL';
  return null;
}

/** A backend that fails once is not used again: the plots it drew are drawn again with the next one, or on the CPU. */
function fail(failedBackend: Backend, error: unknown): void {
  if (failed.has(failedBackend))
    return;
  failed.add(failedBackend);
  console.warn(`Scatter plot: ${failedBackend} rendering failed, falling back`, error);
  if (failedBackend == 'WebGPU')
    disposeWebGPURenderers();
  else
    disposeWebGLRenderers();
}

export function isGpuAvailable(): boolean {
  return backend() != null;
}

export function isGpuRenderValid(sc: DG.ScatterPlotViewer): boolean {
  const current = backend();
  const look = readCheckLook(sc);
  if (!current || !isGpuRenderable(sc, look))
    return false;
  return current == 'WebGPU' || isWebGLRenderValid(sc, look);
}

export function gpuRender(sc: DG.ScatterPlotViewer, show: boolean): void {
  const current = show ? backend() : null;
  if (current != 'WebGL')
    scWebGLRender(sc, false);
  if (current != 'WebGPU')
    scWebGPURender(sc, false, null);
  try {
    if (current == 'WebGL')
      scWebGLRender(sc, true);
    else if (current == 'WebGPU')
      scWebGPURender(sc, true, gpuDevice);
  }
  catch (error) {
    fail(current!, error);
    sc.invalidateCanvas();
  }
}

export async function gpuPointHitTest(sc: DG.ScatterPlotViewer, pt: DG.Point): Promise<number> {
  const current = backend();
  try {
    const row = current == 'WebGPU' ? await scWebGPUPointHitTest(sc, pt) :
      current == 'WebGL' ? scWebGLPointHitTest(sc, pt) : null;
    if (row != null)
      return row;
  }
  catch (error) {
    fail(current!, error);
  }
  return sc.hitTest(pt.x, pt.y);
}
