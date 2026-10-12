/* WebGPU scatter plot renderer. It draws what ../scatterplot-data.ts describes with the WebGL renderer's shapes,
 * ported to WGSL: dots as one-pixel points, other markers as quads, rows read from storage buffers by row index. A
 * frame is encoded and submitted synchronously, so the markers show in the same frame as the axes the core draws. */
import * as DG from 'datagrok-api/dg';
import {Frame, RowList} from '../scatterplot-data';
import {Renderers, ScatterRenderer} from '../scatterplot-renderer';
import {DRAW_SHADER, FLAGS, FRAME_BYTES, PICK_SHADER, PICK_WORKGROUP} from './shaders';

const STORAGE = ['x', 'y', 'colors', 'sizes', 'types', 'selection', 'filtered', 'filteredOut', 'selected'] as const;
type Storage = typeof STORAGE[number];
const PASSES: RowList[] = ['filteredOut', 'filtered', 'selected'];
const MAX_WORKGROUPS = 65535;

interface Pipelines {
  dot: GPURenderPipeline;
  shape: GPURenderPipeline;
  pick: GPUComputePipeline;
  drawLayout: GPUBindGroupLayout;
  pickLayout: GPUBindGroupLayout;
}

const pipelines = new WeakMap<GPUDevice, Pipelines>();

/** Compiled once per device: shaders and pipelines do not depend on the plot. */
function devicePipelines(device: GPUDevice): Pipelines {
  let p = pipelines.get(device);
  if (p)
    return p;
  const storage = (binding: number, visibility: number): GPUBindGroupLayoutEntry =>
    ({binding, visibility, buffer: {type: 'read-only-storage'}});
  const drawLayout = device.createBindGroupLayout({entries: [
    {binding: 0, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, buffer: {type: 'uniform'}},
    ...[1, 2, 3, 4, 5, 6, 7].map((b) => storage(b, GPUShaderStage.VERTEX))]});
  const pickLayout = device.createBindGroupLayout({entries: [
    {binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: {type: 'uniform'}},
    ...[1, 2, 3, 5, 6].map((b) => storage(b, GPUShaderStage.COMPUTE)),
    {binding: 8, visibility: GPUShaderStage.COMPUTE, buffer: {type: 'storage'}}]});
  const draw = device.createShaderModule({code: DRAW_SHADER});
  const format = navigator.gpu.getPreferredCanvasFormat();
  const layout = device.createPipelineLayout({bindGroupLayouts: [drawLayout]});
  const premultiplied: GPUBlendComponent = {srcFactor: 'one', dstFactor: 'one-minus-src-alpha', operation: 'add'};
  p = {
    drawLayout, pickLayout,
    dot: device.createRenderPipeline({layout, primitive: {topology: 'point-list'},
      vertex: {module: draw, entryPoint: 'dotVertex'},
      fragment: {module: draw, entryPoint: 'dotFragment', targets: [{format}]}}),
    shape: device.createRenderPipeline({layout, primitive: {topology: 'triangle-strip'},
      vertex: {module: draw, entryPoint: 'shapeVertex'},
      fragment: {module: draw, entryPoint: 'shapeFragment',
        targets: [{format, blend: {color: premultiplied, alpha: premultiplied}}]}}),
    pick: device.createComputePipeline({layout: device.createPipelineLayout({bindGroupLayouts: [pickLayout]}),
      compute: {module: device.createShaderModule({code: PICK_SHADER}), entryPoint: 'pick'}}),
  };
  pipelines.set(device, p);
  return p;
}

const renderers = new Renderers<ScatterPlotGPU>();

/** Draws [sc] with [device], or releases its renderer when not [show]; answers the canvas drawn on. */
export function scWebGPURender(sc: DG.ScatterPlotViewer, show: boolean,
  device: GPUDevice | null): HTMLCanvasElement | null {
  let renderer = renderers.get(sc);
  if (renderer && (!show || !device || renderer.device !== device || !sc.canvas.isConnected)) {
    renderer.dispose();
    renderer = undefined;
  }
  if (!show || !device || !sc.canvas.isConnected)
    return null;
  renderer ??= new ScatterPlotGPU(sc, device);
  renderer.render(sc);
  return renderer.canvas;
}

/** The row under [pt], or null when [sc] is not drawn with WebGPU and the CPU has to answer. */
export function scWebGPUPointHitTest(sc: DG.ScatterPlotViewer, pt: DG.Point): Promise<number> | null {
  return renderers.get(sc)?.hitTest(sc, pt) ?? null;
}

export function disposeWebGPURenderers(): void {
  renderers.disposeAll();
}


class ScatterPlotGPU extends ScatterRenderer {
  private readonly context: GPUCanvasContext;
  private readonly pipelines: Pipelines;
  private readonly storage = {} as Record<Storage, GPUBuffer>;
  private readonly uniforms: Record<RowList | 'pick', GPUBuffer>;
  private readonly result: GPUBuffer;
  private bindGroups: Record<RowList | 'pick', GPUBindGroup> | null = null;

  constructor(sc: DG.ScatterPlotViewer, readonly device: GPUDevice) {
    super(sc, 'webGPUCanvas', renderers);
    this.pipelines = devicePipelines(device);
    const context = this.canvas.getContext('webgpu') as GPUCanvasContext | null;
    if (!context)
      throw new Error('Scatter plot: no WebGPU canvas context');
    context.configure({device, format: navigator.gpu.getPreferredCanvasFormat(), alphaMode: 'premultiplied'});
    this.context = context;
    for (const name of STORAGE)
      this.storage[name] = this.buffer(16, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
    const uniform = () => this.buffer(FRAME_BYTES, GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST);
    this.uniforms = {filteredOut: uniform(), filtered: uniform(), selected: uniform(), pick: uniform()};
    this.result = this.buffer(4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC | GPUBufferUsage.COPY_DST);
    this.attach();
  }

  protected release(): void {
    for (const buffer of [...Object.values(this.storage), ...Object.values(this.uniforms), this.result])
      buffer.destroy();
    this.context.unconfigure();
  }

  setPositions(axis: 'x' | 'y', values: Float32Array): void {
    this.write(axis, values);
  }

  setRowValues(kind: 'colors' | 'sizes' | 'types', values: Uint32Array | Float32Array): void {
    this.write(kind, values);
  }

  setSelection(bits: Int32Array): void {
    this.write('selection', bits);
  }

  setRows(list: RowList, rows: Uint32Array): void {
    this.write(list, rows);
  }

  render(sc: DG.ScatterPlotViewer): void {
    const frame = this.data.update(sc);
    this.layout(frame);
    for (const list of PASSES)
      this.device.queue.writeBuffer(this.uniforms[list], 0, this.frameUniforms(frame, list));
    const groups = this.groups();
    const encoder = this.device.createCommandEncoder();
    const pass = encoder.beginRenderPass({colorAttachments: [{view: this.context.getCurrentTexture().createView(),
      clearValue: {r: 0, g: 0, b: 0, a: 0}, loadOp: 'clear', storeOp: 'store'}]});
    pass.setPipeline(frame.dots ? this.pipelines.dot : this.pipelines.shape);
    for (const list of PASSES) {
      const count = frame.rows[list];
      if (count == 0)
        continue;
      pass.setBindGroup(0, groups[list]);
      if (frame.dots)
        pass.draw(count);
      else
        pass.draw(4, count);
    }
    pass.end();
    this.device.queue.submit([encoder.finish()]);
  }

  async hitTest(sc: DG.ScatterPlotViewer, pt: DG.Point): Promise<number> {
    const frame = this.data.update(sc);
    const count = frame.rows.filtered;
    if (count == 0)
      return -1;
    this.device.queue.writeBuffer(this.uniforms.pick, 0, this.frameUniforms(frame, 'filtered', pt));
    this.device.queue.writeBuffer(this.result, 0, new Int32Array([-1]));
    const readback = this.buffer(4, GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST);
    const encoder = this.device.createCommandEncoder();
    const pass = encoder.beginComputePass();
    pass.setPipeline(this.pipelines.pick);
    pass.setBindGroup(0, this.groups().pick);
    const workgroups = Math.ceil(count / PICK_WORKGROUP);
    const x = Math.min(workgroups, MAX_WORKGROUPS);
    pass.dispatchWorkgroups(x, Math.ceil(workgroups / x));
    pass.end();
    encoder.copyBufferToBuffer(this.result, 0, readback, 0, 4);
    this.device.queue.submit([encoder.finish()]);
    try {
      await readback.mapAsync(GPUMapMode.READ);
      return new Int32Array(readback.getMappedRange())[0];
    }
    finally {
      readback.destroy();
    }
  }

  private buffer(size: number, usage: number): GPUBuffer {
    return this.device.createBuffer({size, usage});
  }

  private write(name: Storage, data: ArrayBufferView): void {
    const bytes = Math.max(16, Math.ceil(data.byteLength / 4) * 4);
    if (this.storage[name].size < bytes) {
      this.storage[name].destroy();
      this.storage[name] = this.buffer(bytes, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
      this.bindGroups = null;
    }
    if (data.byteLength > 0)
      this.device.queue.writeBuffer(this.storage[name], 0, data.buffer, data.byteOffset, data.byteLength);
  }

  private groups(): Record<RowList | 'pick', GPUBindGroup> {
    if (this.bindGroups)
      return this.bindGroups;
    const s = this.storage;
    const entry = (binding: number, buffer: GPUBuffer): GPUBindGroupEntry => ({binding, resource: {buffer}});
    const draw = (list: RowList) => this.device.createBindGroup({layout: this.pipelines.drawLayout, entries: [
      entry(0, this.uniforms[list]), entry(1, s[list]), entry(2, s.x), entry(3, s.y), entry(4, s.colors),
      entry(5, s.sizes), entry(6, s.types), entry(7, s.selection)]});
    const pick = this.device.createBindGroup({layout: this.pipelines.pickLayout, entries: [
      entry(0, this.uniforms.pick), entry(1, s.filtered), entry(2, s.x), entry(3, s.y), entry(5, s.sizes),
      entry(6, s.types), entry(8, this.result)]});
    return this.bindGroups = {filteredOut: draw('filteredOut'), filtered: draw('filtered'),
      selected: draw('selected'), pick};
  }

  /** The `Frame` struct of the shaders for [list]: filtered-out rows all take the filtered-out color, if any. */
  private frameUniforms(frame: Frame, list: RowList, mouse?: DG.Point): ArrayBuffer {
    const buffer = new ArrayBuffer(FRAME_BYTES);
    const f = new Float32Array(buffer);
    const u = new Uint32Array(buffer);
    const vb = frame.viewBox;
    const outColor = list == 'filteredOut' ? frame.filteredOutColor : null;
    f.set(frame.xform, 0);
    f.set([vb.x, vb.y, vb.x + vb.width, vb.y + vb.height], 4);
    f.set(frame.origin, 8);
    f.set(frame.size, 10);
    f.set([mouse?.x ?? 0, mouse?.y ?? 0], 12);
    f.set([frame.dpr, frame.opacity, frame.border, frame.markerSize], 14);
    u.set([outColor ?? frame.color, frame.markerType, frame.selectedColor], 18);
    u[21] = (frame.rowColors && outColor == null ? FLAGS.rowColors : 0) | (frame.rowSizes ? FLAGS.rowSizes : 0) |
      (frame.rowTypes ? FLAGS.rowTypes : 0) | (frame.selection && outColor == null ? FLAGS.selection : 0) |
      (frame.markerCache ? FLAGS.markerCache : 0);
    u[22] = frame.rows[list];
    return buffer;
  }
}
