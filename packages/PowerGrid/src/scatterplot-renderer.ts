/* What the WebGL and the WebGPU scatter plot renderers share: a canvas over the plot area, placed between the
 * viewer's canvas and its overlay, the plot's data, and a lifetime tied to the viewer. */
import * as DG from 'datagrok-api/dg';
import {Frame, RowList, ScatterData, ScatterSink} from './scatterplot-data';
import {scatterDebug} from './scatterplot-debug';
import {plotCanvases} from './scatterplot-look';

/** The renderers of one backend, one per viewer canvas. */
export class Renderers<R extends ScatterRenderer> {
  private readonly byCanvas = new WeakMap<HTMLCanvasElement, R>();
  private readonly live = new Set<R>();

  get count(): number {
    return this.live.size;
  }

  get(sc: DG.ScatterPlotViewer): R | undefined {
    return this.byCanvas.get(plotCanvases(sc).canvas);
  }

  add(sc: DG.ScatterPlotViewer, renderer: R): void {
    this.byCanvas.set(plotCanvases(sc).canvas, renderer);
    this.live.add(renderer);
  }

  delete(sc: DG.ScatterPlotViewer, renderer: R): void {
    const canvas = plotCanvases(sc).canvas;
    if (this.byCanvas.get(canvas) === renderer)
      this.byCanvas.delete(canvas);
    this.live.delete(renderer);
  }

  /** Drops every drawing and redraws its plot: once a backend is given up, the core asks it to hide nothing. */
  disposeAll(): void {
    for (const renderer of [...this.live])
      renderer.dispose(true);
  }
}

export abstract class ScatterRenderer implements ScatterSink {
  readonly canvas: HTMLCanvasElement;
  protected readonly data: ScatterData;
  /** The frame last drawn: what the hit test answers for, as it is what the mouse is over. */
  protected frame: Frame | null = null;
  /** The backend and the GPU it draws on, for the debug log. */
  protected abstract readonly engine: string;
  private readonly created = performance.now();
  private detachSub: {unsubscribe(): void} | null = null;
  private disposed = false;
  private layoutKey = '';

  constructor(protected readonly sc: DG.ScatterPlotViewer, name: string,
    private readonly renderers: Renderers<ScatterRenderer>) {
    const canvas = this.canvas = document.createElement('canvas');
    canvas.setAttribute('name', name);
    canvas.style.position = 'absolute';
    canvas.style.pointerEvents = 'none';
    this.data = new ScatterData(this);
  }

  /** Shows the canvas and ties the renderer to the viewer, once the backend has its resources. */
  protected attach(): void {
    const {canvas, overlay} = plotCanvases(this.sc);
    canvas.parentElement!.insertBefore(this.canvas, overlay);
    this.renderers.add(this.sc, this);
    this.detachSub = this.sc.onDetached.subscribe(() => this.dispose());
    scatterDebug.engineSetUp(this.sc, this.engine, this.created);
  }

  /** Releases the backend's resources and the canvas, and with [redraw] has the plot drawn again without them. */
  dispose(redraw: boolean = false): void {
    if (this.disposed)
      return;
    this.disposed = true;
    this.renderers.delete(this.sc, this);
    this.data.dispose();
    this.detachSub?.unsubscribe();
    scatterDebug.forget(plotCanvases(this.sc).canvas);
    this.release();
    this.canvas.remove();
    if (redraw)
      this.sc.invalidateCanvas();
  }

  render(sc: DG.ScatterPlotViewer): void {
    const debug = scatterDebug.on;
    const start = debug ? performance.now() : 0;
    const frame = this.data.update(sc);
    const updated = debug ? performance.now() : 0;
    this.layout(frame);
    this.drawFrame(frame);
    this.frame = frame;
    if (debug)
      scatterDebug.frame(sc, frame, this.engine, start, updated, this.gpuDone());
  }

  protected abstract release(): void;
  protected abstract drawFrame(frame: Frame): void;
  /** Settles once the GPU has finished the work submitted so far; the debug log times frames by it. */
  protected abstract gpuDone(): Promise<void>;
  abstract setPositions(axis: 'x' | 'y', values: Float32Array): void;
  abstract setRowValues(kind: 'colors' | 'sizes' | 'types', values: Uint32Array | Float32Array): void;
  abstract setSelection(bits: Int32Array): void;
  abstract setRows(list: RowList, rows: Uint32Array): void;

  /** Sizes the canvas to the frame's plot area, in device pixels, at the frame's place over the viewer canvas. */
  protected layout(frame: Frame): void {
    const key = [...frame.origin, ...frame.size, frame.dpr].join('|');
    if (key == this.layoutKey)
      return;
    this.layoutKey = key;
    const [[left, top], [width, height], dpr] = [frame.origin, frame.size, frame.dpr];
    this.canvas.width = width;
    this.canvas.height = height;
    Object.assign(this.canvas.style, {left: `${left / dpr}px`, top: `${top / dpr}px`,
      width: `${width / dpr}px`, height: `${height / dpr}px`});
  }
}
