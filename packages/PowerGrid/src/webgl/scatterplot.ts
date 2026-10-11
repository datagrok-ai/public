/* WebGL 2 scatter plot renderer: the GPU path where WebGPU is unavailable, for instance blocked by a policy.
 *
 * It reproduces the core's CPU renderers. DOT markers follow the pixel-buffer dot renderer (one opaque device pixel
 * per row, in row order); other markers follow the canvas renderer (shape, border, opacity, filtered-out rows first
 * and selected rows on top). Per-row data lives in vertex buffers addressed by row index and every pass draws a list
 * of row indexes as points, so filtering and selecting rebuild index lists only. */
import * as DG from 'datagrok-api/dg';
import {isDotRenderer, readLook, ScatterLook} from '../scatterplot-look';
import {DOT_FRAGMENT, DOT_VERTEX, PICK_FRAGMENT, PICK_VERTEX, SHAPE_FRAGMENT, SHAPE_VERTEX} from './shaders';

/** Marker names in the order of the core's `Marker.types`: `getMarkerTypes()` answers indexes into this list. */
const MARKER_TYPES = [
  'circle', 'circle border', 'square', 'square border', 'cross border', 'cross x border', 'outlier',
  'diamond', 'diamond border', 'triangle top', 'triangle right', 'triangle bottom', 'triangle left',
  'triangle top border', 'triangle right border', 'triangle bottom border', 'triangle left border',
  'star', 'dot', 'gradient'];
const CIRCLE_BORDER = MARKER_TYPES.indexOf('circle border');
/** The size the core draws a row with an empty Size value at, when *Show Markers With Empty Size* is on. */
const EMPTY_SIZE_MARKER_SIZE = 6;
/** Positions are float32 offsets from an origin near the viewport; once the viewport is this many of its widths
 * away from the origin, the buffer is rebuilt around it so that zooming deep in keeps sub-pixel precision. */
const REBASE_DISTANCE = 100;
const MIN_LOG = 1e-30;
const NO_ROWS = new Uint32Array(0);
/** Browsers keep about 16 WebGL contexts a page and drop the oldest past that; leave room for other viewers. */
const MAX_RENDERERS = 8;

const enum Attr { X, Y, Color, Size, Type, Selected }

/** Without a GPU, a software rasterizer draws markers slower than the core's canvas renderer does. */
const CONTEXT_ATTRIBUTES: WebGLContextAttributes = {alpha: true, premultipliedAlpha: true, antialias: false,
  depth: false, stencil: false, preserveDrawingBuffer: false, failIfMajorPerformanceCaveat: true};

let capabilities: {maxPointSize: number} | null | undefined;

function webGLCapabilities() {
  if (capabilities === undefined) {
    const gl = document.createElement('canvas').getContext('webgl2', CONTEXT_ATTRIBUTES);
    // the browser flags its own software fallback only; one picked by a flag or a policy shows in the renderer name
    const info = gl?.getExtension('WEBGL_debug_renderer_info');
    const software = info && /swiftshader|llvmpipe|softpipe|software|basic render/i
      .test(gl!.getParameter(info.UNMASKED_RENDERER_WEBGL));
    capabilities = gl && !software ? {maxPointSize: gl.getParameter(gl.ALIASED_POINT_SIZE_RANGE)[1]} : null;
    gl?.getExtension('WEBGL_lose_context')?.loseContext();
  }
  return capabilities;
}

export function isWebGLAvailable(): boolean {
  return webGLCapabilities() != null;
}

function largestMarker(sc: DG.ScatterPlotViewer, look: ScatterLook): number {
  return Math.max(look.sizeColumnName ? look.markerMaxSize : sc.getMarkerSize(0), EMPTY_SIZE_MARKER_SIZE);
}

/** The point sprite of a marker, as `SHAPE_VERTEX` sizes it for the widest stroke. */
function spriteSize(markerSize: number, border: number, dpr: number): number {
  return 2 * Math.ceil((markerSize / 2 + Math.max(border, 3)) * dpr + 1);
}

/** The largest marker fits in a point sprite, and a context is free for a plot that has none yet. */
export function isWebGLRenderValid(sc: DG.ScatterPlotViewer, look: ScatterLook): boolean {
  const sprite = spriteSize(largestMarker(sc, look), look.markerBorderWidth, window.devicePixelRatio);
  return sprite <= webGLCapabilities()!.maxPointSize && (renderers.has(sc.canvas) || live.size < MAX_RENDERERS);
}

const renderers = new WeakMap<HTMLCanvasElement, ScatterPlotGL>();
const live = new Set<ScatterPlotGL>();

/** Draws [sc], or releases its renderer when not [show]; answers the canvas drawn on. */
export function scWebGLRender(sc: DG.ScatterPlotViewer, show: boolean): HTMLCanvasElement | null {
  let renderer = renderers.get(sc.canvas);
  if (!show || !sc.canvas.isConnected) {
    renderer?.dispose();
    return null;
  }
  if (!renderer) {
    renderer = new ScatterPlotGL(sc);
    renderers.set(sc.canvas, renderer);
  }
  renderer.render(sc);
  return renderer.canvas;
}

/** The row under [pt], or null when [pt] is off the GPU canvas and the CPU has to answer. */
export function scWebGLPointHitTest(sc: DG.ScatterPlotViewer, pt: DG.Point): number | null {
  const renderer = renderers.get(sc.canvas);
  return renderer ? renderer.hitTest(sc, pt) : null;
}

/** Drops every WebGL drawing and redraws its plot: once WebGL is given up, the core asks it to hide nothing. */
export function disposeWebGLRenderers(): void {
  for (const renderer of [...live])
    renderer.dispose(true);
}

const objectIds = new WeakMap<object, number>();
let lastObjectId = 0;

/** Tells apart columns and bitsets that share a name and a version: versions start over on every new object. */
function objectId(x: {dart: object} | null | undefined): number {
  if (x == null)
    return 0;
  let id = objectIds.get(x.dart);
  if (id === undefined)
    objectIds.set(x.dart, id = ++lastObjectId);
  return id;
}


class Program {
  readonly program: WebGLProgram;
  private readonly locations = new Map<string, WebGLUniformLocation | null>();

  constructor(private readonly gl: WebGL2RenderingContext, vertex: string, fragment: string) {
    const program = gl.createProgram()!;
    for (const [type, source] of [[gl.VERTEX_SHADER, vertex], [gl.FRAGMENT_SHADER, fragment]] as [number, string][]) {
      const shader = gl.createShader(type)!;
      gl.shaderSource(shader, source);
      gl.compileShader(shader);
      if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS) && !gl.isContextLost())
        throw new Error(`Scatter plot WebGL shader: ${gl.getShaderInfoLog(shader)}`);
      gl.attachShader(program, shader);
    }
    gl.linkProgram(program);
    if (!gl.getProgramParameter(program, gl.LINK_STATUS) && !gl.isContextLost())
      throw new Error(`Scatter plot WebGL program: ${gl.getProgramInfoLog(program)}`);
    this.program = program;
  }

  uniform(name: string): WebGLUniformLocation | null {
    if (!this.locations.has(name))
      this.locations.set(name, this.gl.getUniformLocation(this.program, name));
    return this.locations.get(name)!;
  }
}


/** One axis: row positions as float32 offsets from [origin], in log space on a log axis. */
class Axis {
  key = '';
  origin = 0;
  log = false;
  /** 1 where the row has no position: empty, infinite, or not positive on a log axis. */
  invalid = new Uint8Array(0);
  version = 0;

  constructor(readonly buffer: WebGLBuffer) {}

  /** Rebuilds the buffer when the column changed, or when the viewport [min, max] drifted far from the origin. */
  update(gl: WebGL2RenderingContext, col: DG.Column, log: boolean, min: number, max: number): void {
    const categorical = col.type == DG.COLUMN_TYPE.STRING;
    const key = `${objectId(col)}|${col.version}|${col.length}|${log}`;
    this.log = log;
    const [lo, hi] = this.range(min, max);
    const center = (lo + hi) / 2;
    if (key == this.key && (categorical || !(Math.abs(center - this.origin) > REBASE_DISTANCE * (hi - lo))))
      return;

    this.key = key;
    this.origin = categorical || !Number.isFinite(center) ? 0 : center;
    this.version++;
    const n = col.length;
    const raw = col.getRawData();
    const values = new Float32Array(n);
    const invalid = this.invalid.length == n ? this.invalid : new Uint8Array(n);
    if (categorical) {
      for (let i = 0; i < n; i++)
        values[i] = raw[i];
      invalid.fill(0);
    }
    else {
      const none = noneValue(col, raw);
      for (let i = 0; i < n; i++) {
        const v = raw[i] as number;
        const bad = v === none || !Number.isFinite(v) || (log && v <= 0);
        invalid[i] = bad ? 1 : 0;
        values[i] = bad ? 0 : (log ? Math.log(Math.max(v, MIN_LOG)) : v) - this.origin;
      }
    }
    this.invalid = invalid;
    gl.bindBuffer(gl.ARRAY_BUFFER, this.buffer);
    gl.bufferData(gl.ARRAY_BUFFER, values, gl.STATIC_DRAW);
  }

  /** [start + scale * position] places a row on the screen as the core's `Transform.worldToScreen` does. */
  transform(min: number, max: number, start: number, length: number, inverse: boolean): [number, number] {
    const [lo, hi] = this.range(min, max);
    const scale = length / (hi - lo);
    const offset = (this.origin - lo) * scale;
    return inverse ? [start + length - offset, -scale] : [start + offset, scale];
  }

  private range(min: number, max: number): [number, number] {
    return this.log ? [Math.log(Math.max(min, MIN_LOG)), Math.log(max)] : [min, max];
  }
}

/** The raw value of an empty cell; `col.getRawData()` of a datetime column is float64, of a double column float32. */
function noneValue(col: DG.Column, raw: ArrayLike<number>): number {
  if (col.type == DG.COLUMN_TYPE.INT)
    return DG.INT_NULL;
  return raw instanceof Float32Array ? Math.fround(DG.FLOAT_NULL) : DG.FLOAT_NULL;
}


class RowList {
  key = '';
  rows: Uint32Array = NO_ROWS;
  version = 0;

  constructor(readonly buffer: WebGLBuffer) {}

  set(gl: WebGL2RenderingContext, key: string, rows: Uint32Array): void {
    this.key = key;
    this.rows = rows;
    this.version++;
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.buffer);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, rows, gl.STATIC_DRAW);
  }
}


/** A per-row attribute, or one value for all rows. */
class RowAttribute {
  key = '';
  perRow = false;
  constant: number[] = [0, 0, 0, 0];

  constructor(readonly buffer: WebGLBuffer, readonly index: Attr, readonly size: number, readonly type: number,
    readonly normalized: boolean) {}

  setRows(gl: WebGL2RenderingContext, key: string, data: ArrayBufferView): void {
    this.key = key;
    this.perRow = true;
    gl.bindBuffer(gl.ARRAY_BUFFER, this.buffer);
    gl.bufferData(gl.ARRAY_BUFFER, data, gl.STATIC_DRAW);
  }

  setConstant(key: string, ...values: number[]): void {
    this.key = key;
    this.perRow = false;
    this.constant = values;
  }

  bind(gl: WebGL2RenderingContext): void {
    if (this.perRow) {
      gl.enableVertexAttribArray(this.index);
      gl.bindBuffer(gl.ARRAY_BUFFER, this.buffer);
      gl.vertexAttribPointer(this.index, this.size, this.type, this.normalized, 0, 0);
    }
    else {
      const [a, b = 0, c = 0, d = 0] = this.constant;
      gl.disableVertexAttribArray(this.index);
      gl.vertexAttrib4f(this.index, a, b, c, d);
    }
  }
}

/** The core's ARGB color as the bytes a per-row color buffer holds: (b, g, r, a). */
function colorBytes(argb: number): [number, number, number, number] {
  return [argb & 0xFF, (argb >>> 8) & 0xFF, (argb >>> 16) & 0xFF, argb >>> 24].map((c) => c / 255) as
    [number, number, number, number];
}


interface Frame {
  look: ScatterLook;
  dots: boolean;
  dpr: number;
  viewBox: DG.Rect;
  xform: [number, number, number, number];
  /** The canvas over the plot area, in device pixels of the viewer canvas, snapped to them. */
  origin: [number, number];
  size: [number, number];
}


class ScatterPlotGL {
  readonly canvas: HTMLCanvasElement;
  private readonly gl: WebGL2RenderingContext;
  private readonly dotProgram: Program;
  private readonly shapeProgram: Program;
  private readonly pickProgram: Program;
  private readonly vao: WebGLVertexArrayObject;
  private readonly x: Axis;
  private readonly y: Axis;
  private readonly colors: RowAttribute;
  private readonly sizes: RowAttribute;
  private readonly types: RowAttribute;
  /** 1 for a selected row, when selected rows take the selection color. */
  private readonly selectedFlags: RowAttribute;
  private selectionColors = false;
  /** The selection when the row colors were taken: the core painted these rows in the selection color. */
  private colorsSelection: Int32Array | null = null;
  private readonly filtered: RowList;
  private readonly filteredOut: RowList;
  private readonly selected: RowList;
  private readonly pickTarget: WebGLFramebuffer;
  /** Size column state: 1 where the size is empty, 2 where it is not positive on a log scale. */
  private sizeFlags = new Uint8Array(0);
  private sizeFlagsKey = '';
  private sizeFlagsVersion = 0;
  private emptySizes = false;
  private dataFrame: DG.DataFrame | null = null;
  private metadataSub: {unsubscribe(): void} | null = null;
  private metadataVersion = 0;
  private readonly detachSub: {unsubscribe(): void};
  private disposed = false;
  private layoutKey = '';

  constructor(private readonly sc: DG.ScatterPlotViewer) {
    const canvas = this.canvas = document.createElement('canvas');
    canvas.setAttribute('name', 'webGLCanvas');
    canvas.style.position = 'absolute';
    canvas.style.pointerEvents = 'none';
    sc.canvas.parentElement!.insertBefore(canvas, sc.overlay);
    const gl = canvas.getContext('webgl2', CONTEXT_ATTRIBUTES);
    if (!gl) {
      canvas.remove();
      throw new Error('Scatter plot: no WebGL 2 context');
    }
    this.gl = gl;
    live.add(this);

    this.dotProgram = new Program(gl, DOT_VERTEX, DOT_FRAGMENT);
    this.shapeProgram = new Program(gl, SHAPE_VERTEX, SHAPE_FRAGMENT);
    this.pickProgram = new Program(gl, PICK_VERTEX, PICK_FRAGMENT);
    this.vao = gl.createVertexArray()!;
    const buffer = () => gl.createBuffer()!;
    this.x = new Axis(buffer());
    this.y = new Axis(buffer());
    this.colors = new RowAttribute(buffer(), Attr.Color, 4, gl.UNSIGNED_BYTE, true);
    this.sizes = new RowAttribute(buffer(), Attr.Size, 1, gl.FLOAT, false);
    this.types = new RowAttribute(buffer(), Attr.Type, 1, gl.UNSIGNED_INT, false);
    this.selectedFlags = new RowAttribute(buffer(), Attr.Selected, 1, gl.UNSIGNED_BYTE, false);
    this.filtered = new RowList(buffer());
    this.filteredOut = new RowList(buffer());
    this.selected = new RowList(buffer());

    this.pickTarget = gl.createFramebuffer()!;
    const pickPixel = gl.createRenderbuffer();
    gl.bindRenderbuffer(gl.RENDERBUFFER, pickPixel);
    gl.renderbufferStorage(gl.RENDERBUFFER, gl.RGBA8, 1, 1);
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.pickTarget);
    gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.RENDERBUFFER, pickPixel);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);

    // a lost context (driver reset, too many contexts) cannot draw again: the next frame starts on a new canvas
    canvas.addEventListener('webglcontextlost', () => this.dispose(true));
    this.detachSub = sc.onDetached.subscribe(() => this.dispose());
  }

  /** Releases the context and the canvas, and with [redraw] has the plot drawn again without them. */
  dispose(redraw: boolean = false): void {
    if (this.disposed)
      return;
    this.disposed = true;
    live.delete(this);
    if (renderers.get(this.sc.canvas) === this)
      renderers.delete(this.sc.canvas);
    this.metadataSub?.unsubscribe();
    this.detachSub.unsubscribe();
    this.gl.getExtension('WEBGL_lose_context')?.loseContext();
    this.canvas.remove();
    if (redraw)
      this.sc.invalidateCanvas();
  }

  render(sc: DG.ScatterPlotViewer): void {
    const gl = this.gl;
    const frame = this.prepare(sc);
    this.layout(frame);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.viewport(0, 0, frame.size[0], frame.size[1]);
    gl.clearColor(0, 0, 0, 0);
    gl.clear(gl.COLOR_BUFFER_BIT);
    gl.bindVertexArray(this.vao);
    this.bindPositions();
    this.selectedFlags.bind(gl);

    const look = frame.look;
    const outColor = look.rowSource == 'All' ? null : look.filteredOutRowsColor;
    if (frame.dots) {
      gl.disable(gl.BLEND);
      const program = this.use(this.dotProgram, frame);
      this.draw(program, this.filteredOut, outColor);
      this.draw(program, this.filtered, null);
    }
    else {
      gl.enable(gl.BLEND);
      gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
      const program = this.use(this.shapeProgram, frame);
      const viewer = sc.canvas;
      const opacity = Math.trunc(Math.min(Math.max(look.markerOpacity, 0), 100) * 255 / 100) / 255;
      gl.uniform1i(program.uniform('u_markerCache'), +(viewer.width > 250 && viewer.height > 250));
      gl.uniform1f(program.uniform('u_border'), look.markerDrawBorder ? look.markerBorderWidth : 0);
      gl.uniform1f(program.uniform('u_alpha'), opacity);
      this.sizes.bind(gl);
      this.types.bind(gl);
      this.draw(program, this.filteredOut, outColor);
      this.draw(program, this.filtered, null);
      this.draw(program, this.selected, null);
    }
    gl.bindVertexArray(null);
  }

  hitTest(sc: DG.ScatterPlotViewer, pt: DG.Point): number | null {
    const gl = this.gl;
    const frame = this.prepare(sc);
    // GL drops a point whose center is off the viewport, so the pick viewport reaches past the canvas by the largest
    // hit radius: a marker centered off the plot still takes the mouse on it, as on the CPU
    const margin = Math.ceil((largestMarker(sc, frame.look) / 2 + 2) * frame.dpr);
    frame.origin = [frame.origin[0] - margin, frame.origin[1] - margin];
    frame.size = [frame.size[0] + 2 * margin, frame.size[1] + 2 * margin];
    const px = Math.floor(pt.x * frame.dpr) - frame.origin[0];
    const py = Math.floor(pt.y * frame.dpr) - frame.origin[1];
    if (px < 0 || py < 0 || px >= frame.size[0] || py >= frame.size[1])
      return null;

    // a one-pixel target, with the viewport shifted so that the mouse pixel lands on it
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.pickTarget);
    gl.viewport(-px, py - frame.size[1] + 1, frame.size[0], frame.size[1]);
    gl.disable(gl.BLEND);
    gl.clearColor(0, 0, 0, 0);
    gl.clear(gl.COLOR_BUFFER_BIT);
    gl.bindVertexArray(this.vao);
    this.bindPositions();
    const program = this.use(this.pickProgram, frame);
    gl.uniform2f(program.uniform('u_mouse'), pt.x, pt.y);
    this.sizes.bind(gl);
    this.types.bind(gl);
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.filtered.buffer);
    gl.drawElements(gl.POINTS, this.filtered.rows.length, gl.UNSIGNED_INT, 0);
    const pixel = new Uint8Array(4);
    gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, pixel);
    gl.bindVertexArray(null);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    return (pixel[0] | (pixel[1] << 8) | (pixel[2] << 16) | (pixel[3] << 24)) - 1;
  }

  private bindPositions(): void {
    const gl = this.gl;
    for (const [axis, index] of [[this.x, Attr.X], [this.y, Attr.Y]] as [Axis, Attr][]) {
      gl.enableVertexAttribArray(index);
      gl.bindBuffer(gl.ARRAY_BUFFER, axis.buffer);
      gl.vertexAttribPointer(index, 1, gl.FLOAT, false, 0, 0);
    }
  }

  private use(program: Program, frame: Frame): Program {
    const gl = this.gl;
    const vb = frame.viewBox;
    gl.useProgram(program.program);
    gl.uniform4f(program.uniform('u_xform'), ...frame.xform);
    gl.uniform4f(program.uniform('u_viewBox'), vb.x, vb.y, vb.x + vb.width, vb.y + vb.height);
    gl.uniform2f(program.uniform('u_origin'), ...frame.origin);
    gl.uniform2f(program.uniform('u_size'), ...frame.size);
    gl.uniform1f(program.uniform('u_dpr'), frame.dpr);
    const [b, g, r] = colorBytes(frame.look.selectedRowsColor);
    gl.uniform3f(program.uniform('u_selectedColor'), r, g, b);
    return program;
  }

  /** Draws [list] in [color] when given, else in the row colors. */
  private draw(program: Program, list: RowList, color: number | null): void {
    if (list.rows.length == 0)
      return;
    const gl = this.gl;
    if (color == null)
      this.colors.bind(gl);
    else {
      gl.disableVertexAttribArray(Attr.Color);
      gl.vertexAttrib4f(Attr.Color, ...colorBytes(color));
    }
    gl.uniform1i(program.uniform('u_selection'), +(color == null && this.selectionColors));
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, list.buffer);
    gl.drawElements(gl.POINTS, list.rows.length, gl.UNSIGNED_INT, 0);
  }

  private prepare(sc: DG.ScatterPlotViewer): Frame {
    const gl = this.gl;
    const look = readLook(sc);
    const df = sc.dataFrame;
    const dpr = window.devicePixelRatio;
    const viewBox = sc.viewBox;
    const viewport = sc.viewport;
    this.watchMetadata(df);

    const xCol = df.col(look.xColumnName)!;
    const yCol = df.col(look.yColumnName)!;
    const xLog = xCol.isNumerical && look.xAxisType == DG.AxisType.logarithmic;
    const yLog = yCol.isNumerical && look.yAxisType == DG.AxisType.logarithmic;
    this.x.update(gl, xCol, xLog, viewport.x, viewport.x + viewport.width);
    this.y.update(gl, yCol, yLog, viewport.y, viewport.y + viewport.height);
    const xform = [
      ...this.x.transform(viewport.x, viewport.x + viewport.width, viewBox.x, viewBox.width, look.invertXAxis),
      // the core's Y axis grows upwards unless inverted
      ...this.y.transform(viewport.y, viewport.y + viewport.height, viewBox.y, viewBox.height, !look.invertYAxis),
    ] as [number, number, number, number];

    const dots = isDotRenderer(look);
    const sizeCol = look.sizeColumnName ? df.col(look.sizeColumnName) : null;
    this.updateSizeFlags(sizeCol, sizeCol != null && sizeCol.isNumerical &&
      look.markerSizeScaling == DG.AxisType.logarithmic);
    // the dot renderer ignores sizes and types, the hit test does not
    this.updateSizes(sc, look, sizeCol);
    this.updateTypes(sc, look, df);
    this.updateColors(sc, look, df);
    this.updateRows(sc, look, df, dots);

    const left = Math.floor(viewBox.x * dpr);
    const top = Math.floor(viewBox.y * dpr);
    const size: [number, number] = [Math.max(1, Math.ceil((viewBox.x + viewBox.width) * dpr) - left),
      Math.max(1, Math.ceil((viewBox.y + viewBox.height) * dpr) - top)];
    return {look, dots, dpr, viewBox, xform, origin: [left, top], size};
  }

  private layout(frame: Frame): void {
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

  /** Column color coding lives in tags, which change no column version. */
  private watchMetadata(df: DG.DataFrame): void {
    if (this.dataFrame?.dart === df.dart)
      return;
    this.metadataSub?.unsubscribe();
    this.dataFrame = df;
    this.metadataSub = df.onMetadataChanged.subscribe(() => this.metadataVersion++);
    this.metadataVersion++;
  }

  private updateSizeFlags(sizeCol: DG.Column | null, log: boolean): void {
    const key = sizeCol ? `${objectId(sizeCol)}|${sizeCol.version}|${sizeCol.length}|${log}` : '';
    if (key == this.sizeFlagsKey)
      return;
    this.sizeFlagsKey = key;
    this.sizeFlagsVersion++;
    if (!sizeCol) {
      this.sizeFlags = new Uint8Array(0);
      this.emptySizes = false;
      return;
    }
    const n = sizeCol.length;
    const flags = new Uint8Array(n);
    const raw = sizeCol.getRawData();
    const fast = sizeCol.type == DG.COLUMN_TYPE.INT || sizeCol.type == DG.COLUMN_TYPE.FLOAT;
    const none = noneValue(sizeCol, raw as ArrayLike<number>);
    let emptySizes = false;
    for (let i = 0; i < n; i++) {
      const empty = fast ? raw[i] === none : sizeCol.isNone(i);
      emptySizes ||= empty;
      flags[i] = empty ? 1 : log && (fast ? raw[i] as number : sizeCol.getNumber(i)) <= 0 ? 2 : 0;
    }
    this.sizeFlags = flags;
    this.emptySizes = emptySizes;
  }

  private updateSizes(sc: DG.ScatterPlotViewer, look: ScatterLook, sizeCol: DG.Column | null): void {
    if (!sizeCol) {
      const size = sc.getMarkerSize(0);
      this.sizes.setConstant(`${size}`, size);
      return;
    }
    const key = `${this.sizeFlagsKey}|${look.markerMinSize}|${look.markerMaxSize}`;
    if (key != this.sizes.key || !this.sizes.perRow)
      this.sizes.setRows(this.gl, key, sc.getMarkerSizes());
  }

  private updateTypes(sc: DG.ScatterPlotViewer, look: ScatterLook, df: DG.DataFrame): void {
    const markersCol = look.markersColumnName ? df.col(look.markersColumnName) : null;
    const emptySizes = look.showMarkersWithEmptySize && this.emptySizes;
    const type = Math.max(MARKER_TYPES.indexOf(look.markerType), 0);
    if (!markersCol && !emptySizes) {
      this.types.setConstant(`${type}`, type);
      return;
    }
    const key = [objectId(markersCol), markersCol?.version, look.markersMap, this.metadataVersion, type, emptySizes,
      this.sizeFlagsVersion].join('|');
    if (key == this.types.key && this.types.perRow)
      return;
    const types = markersCol ? sc.getMarkerTypes() : new Uint32Array(df.rowCount).fill(type);
    if (emptySizes) {
      for (let i = 0; i < types.length; i++) {
        if (this.sizeFlags[i] == 1)
          types[i] = CIRCLE_BORDER;
      }
    }
    this.types.setRows(this.gl, key, types);
  }

  /** The core paints selected rows in the selection color; the shaders do it from per-row flags instead, so a new
   * selection takes the core's colors again only when it deselects a row the core painted selected. */
  private updateColors(sc: DG.ScatterPlotViewer, look: ScatterLook, df: DG.DataFrame): void {
    const gl = this.gl;
    const selection = df.selection;
    this.selectionColors = look.showSelectedRows && look.rowSource != 'Selected' &&
      look.rowSource != 'FilteredSelected';
    const flagsKey = this.selectionColors && selection.trueCount > 0 ?
      `${objectId(selection)}|${selection.version}|${selection.length}` : '';
    if (flagsKey != this.selectedFlags.key) {
      if (flagsKey) {
        const bits = selection.getBuffer();
        const flags = new Uint8Array(selection.length);
        for (let i = 0; i < flags.length; i++)
          flags[i] = (bits[i >>> 5] >>> (i & 31)) & 1;
        this.selectedFlags.setRows(gl, flagsKey, flags);
      }
      else
        this.selectedFlags.setConstant(flagsKey, 0);
    }

    if (!look.colorColumnName) {
      this.colors.setConstant(`${look.filteredRowsColor}`, ...colorBytes(look.filteredRowsColor));
      return;
    }
    const colorCol = df.col(look.colorColumnName);
    const key = [objectId(colorCol), colorCol?.version, look.colorColumnName, look.colorMap, look.showSelectedRows,
      look.rowSource, look.selectedRowsColor, look.filteredRowsColor, look.missingValueColor, look.colorAxisType,
      look.colorMin, look.colorMax, look.invertColorScheme, look.linearColorScheme, look.categoricalColorScheme,
      this.metadataVersion, objectId(selection), df.rowCount].join('|');
    const bits = selection.getBuffer();
    const taken = this.colorsSelection;
    let deselected = false;
    for (let w = 0; taken && w < taken.length && !deselected; w++)
      deselected = (taken[w] & ~bits[w]) != 0;
    if (key == this.colors.key && this.colors.perRow && !deselected)
      return;
    this.colors.setRows(gl, key, sc.getMarkerColors());
    this.colorsSelection = this.selectionColors ? bits.slice(0, (selection.length + 31) >>> 5) : null;
  }

  private updateRows(sc: DG.ScatterPlotViewer, look: ScatterLook, df: DG.DataFrame, dots: boolean): void {
    const gl = this.gl;
    const filter = sc.filter;
    const skipEmptySize = !dots && !look.showMarkersWithEmptySize;
    const key = [objectId(filter), filter.version, filter.length, this.x.version, this.y.version,
      this.sizeFlagsVersion, look.showFilteredOutPoints, skipEmptySize].join('|');
    if (key != this.filtered.key) {
      const n = filter.length;
      const bits = filter.getBuffer();
      const xInvalid = this.x.invalid;
      const yInvalid = this.y.invalid;
      const sizeFlags = this.sizeFlags;
      const hasSize = sizeFlags.length == n;
      const filtered = new Uint32Array(n);
      const filteredOut = new Uint32Array(look.showFilteredOutPoints ? n : 0);
      let f = 0;
      let o = 0;
      for (let i = 0; i < n; i++) {
        if (xInvalid[i] || yInvalid[i] || (hasSize && (sizeFlags[i] == 2 || (skipEmptySize && sizeFlags[i] == 1))))
          continue;
        if ((bits[i >>> 5] >>> (i & 31)) & 1)
          filtered[f++] = i;
        else if (filteredOut.length)
          filteredOut[o++] = i;
      }
      this.filtered.set(gl, key, filtered.subarray(0, f));
      this.filteredOut.set(gl, key, filteredOut.subarray(0, o));
    }

    // the canvas renderer draws the selected rows once more, on top of the rest; the dot renderer does not
    const selection = df.selection;
    const onTop = !dots && this.selectionColors && selection.trueCount > 0;
    const selectedKey = onTop ? `${this.filtered.version}|${objectId(selection)}|${selection.version}` : '';
    if (selectedKey != this.selected.key) {
      const bits = selection.getBuffer();
      const rows = onTop ? this.filtered.rows.filter((i) => (bits[i >>> 5] >>> (i & 31)) & 1) : NO_ROWS;
      this.selected.set(gl, selectedKey, rows);
    }
  }
}
