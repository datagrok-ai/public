/* WebGL 2 scatter plot renderer: the GPU path where WebGPU is unavailable, for instance blocked by a policy.
 *
 * Per-row data lives in vertex buffers addressed by row index, and every pass draws a list of row indexes as point
 * sprites; the fragment shader shapes the markers as the core's `Marker.draw` does. */
import * as DG from 'datagrok-api/dg';
import {ScatterLook} from '../scatterplot-look';
import {Frame, largestMarker, RowList} from '../scatterplot-data';
import {Renderers, ScatterRenderer} from '../scatterplot-renderer';
import {ANY_KIND, DOT_FRAGMENT, DOT_VERTEX, PICK_FRAGMENT, PICK_VERTEX, SHAPE_VERTEX, shapeFragment} from './shaders';

/** Browsers keep about 16 WebGL contexts a page and drop the oldest past that; leave room for other viewers. */
const MAX_RENDERERS = 8;

const enum Attr { X, Y, Color, Size, Type, Selected }

/** Without a GPU, a software rasterizer draws markers slower than the core's canvas renderer does. */
const CONTEXT_ATTRIBUTES: WebGLContextAttributes = {alpha: true, premultipliedAlpha: true, antialias: false,
  depth: false, stencil: false, preserveDrawingBuffer: false, failIfMajorPerformanceCaveat: true};

let capabilities: {maxPointSize: number, renderer: string} | null | undefined;
/** Why the plots are not drawn with WebGL 2 here, once the probe found it unusable. */
let unavailable = '';

function webGLCapabilities() {
  if (capabilities === undefined) {
    const gl = document.createElement('canvas').getContext('webgl2', CONTEXT_ATTRIBUTES);
    // the browser flags its own software fallback only; one picked by a flag or a policy shows in the renderer name
    const info = gl?.getExtension('WEBGL_debug_renderer_info');
    const renderer: string = info ? gl!.getParameter(info.UNMASKED_RENDERER_WEBGL) : 'an unnamed GPU';
    const software = info && /swiftshader|llvmpipe|softpipe|software|basic render/i.test(renderer);
    capabilities = gl && !software ? {maxPointSize: gl.getParameter(gl.ALIASED_POINT_SIZE_RANGE)[1], renderer} : null;
    unavailable = !gl ? 'no hardware WebGL 2' : software ? `software WebGL 2 (${renderer})` : '';
    gl?.getExtension('WEBGL_lose_context')?.loseContext();
  }
  return capabilities;
}

export function isWebGLAvailable(): boolean {
  return webGLCapabilities() != null;
}

export function webGLUnavailable(): string {
  webGLCapabilities();
  return unavailable;
}

/** The point sprite of a marker, as `SHAPE_VERTEX` sizes it for the widest stroke. */
function spriteSize(markerSize: number, border: number, dpr: number): number {
  return 2 * Math.ceil((markerSize / 2 + Math.max(border, 3)) * dpr + 1);
}

const renderers = new Renderers<ScatterPlotGL>();

/** Why WebGL cannot draw [sc], or null: the largest marker has to fit in a point sprite, and a context be free for a
 * plot that has none yet. */
export function whyNotWebGL(sc: DG.ScatterPlotViewer, look: ScatterLook): string | null {
  const sprite = spriteSize(largestMarker(sc, look), look.markerBorderWidth, window.devicePixelRatio);
  const maxPointSize = webGLCapabilities()!.maxPointSize;
  if (sprite > maxPointSize)
    return `markers of ${sprite} px exceed the largest WebGL point sprite, ${maxPointSize} px`;
  return renderers.get(sc) != null || renderers.count < MAX_RENDERERS ? null :
    `${MAX_RENDERERS} other plots already draw with WebGL`;
}

/** Draws [sc], or releases its renderer when not [show]; answers the canvas drawn on. */
export function scWebGLRender(sc: DG.ScatterPlotViewer, show: boolean): HTMLCanvasElement | null {
  let renderer = renderers.get(sc);
  if (!show || !sc.canvas.isConnected) {
    renderer?.dispose();
    return null;
  }
  renderer ??= new ScatterPlotGL(sc);
  renderer.render(sc);
  return renderer.canvas;
}

/** The row under [pt], or null when [pt] is off the GPU canvas and the CPU has to answer. */
export function scWebGLPointHitTest(sc: DG.ScatterPlotViewer, pt: DG.Point): number | null {
  return renderers.get(sc)?.hitTest(sc, pt) ?? null;
}

export function disposeWebGLRenderers(): void {
  renderers.disposeAll();
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

/** The core's ARGB color as the bytes a per-row color buffer holds: (b, g, r, a). */
function colorBytes(argb: number): [number, number, number, number] {
  return [argb & 0xFF, (argb >>> 8) & 0xFF, (argb >>> 16) & 0xFF, argb >>> 24].map((c) => c / 255) as
    [number, number, number, number];
}


class ScatterPlotGL extends ScatterRenderer {
  protected readonly engine = `WebGL 2 on ${webGLCapabilities()!.renderer}`;
  private readonly gl: WebGL2RenderingContext;
  private readonly dotProgram: Program;
  /** By the marker kind they fix, or `ANY_KIND`; compiled on first use. */
  private readonly shapePrograms = new Map<number, Program>();
  private readonly pickProgram: Program;
  private readonly vao: WebGLVertexArrayObject;
  private readonly buffers: Record<'x' | 'y' | 'colors' | 'sizes' | 'types' | 'selection', WebGLBuffer>;
  private readonly rowBuffers: Record<RowList, WebGLBuffer>;
  private readonly pickTarget: WebGLFramebuffer;

  constructor(sc: DG.ScatterPlotViewer) {
    super(sc, 'webGLCanvas', renderers);
    const gl = this.canvas.getContext('webgl2', CONTEXT_ATTRIBUTES);
    if (!gl)
      throw new Error('Scatter plot: no WebGL 2 context');
    this.gl = gl;
    this.dotProgram = new Program(gl, DOT_VERTEX, DOT_FRAGMENT);
    this.pickProgram = new Program(gl, PICK_VERTEX, PICK_FRAGMENT);
    this.vao = gl.createVertexArray()!;
    const buffer = () => gl.createBuffer()!;
    this.buffers = {x: buffer(), y: buffer(), colors: buffer(), sizes: buffer(), types: buffer(), selection: buffer()};
    this.rowBuffers = {filtered: buffer(), filteredOut: buffer(), selected: buffer()};

    this.pickTarget = gl.createFramebuffer()!;
    const pickPixel = gl.createRenderbuffer();
    gl.bindRenderbuffer(gl.RENDERBUFFER, pickPixel);
    gl.renderbufferStorage(gl.RENDERBUFFER, gl.RGBA8, 1, 1);
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.pickTarget);
    gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.RENDERBUFFER, pickPixel);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);

    // a lost context (driver reset, too many contexts) cannot draw again: the next frame starts on a new canvas
    this.canvas.addEventListener('webglcontextlost', () => this.dispose(true));
    this.attach();
  }

  protected release(): void {
    this.gl.getExtension('WEBGL_lose_context')?.loseContext();
  }

  setPositions(axis: 'x' | 'y', values: Float32Array): void {
    this.upload(this.buffers[axis], values);
  }

  setRowValues(kind: 'colors' | 'sizes' | 'types', values: Uint32Array | Float32Array): void {
    this.upload(this.buffers[kind], values);
  }

  setSelection(bits: Int32Array): void {
    const flags = new Uint8Array(bits.length * 32);
    for (let i = 0; i < flags.length; i++)
      flags[i] = (bits[i >>> 5] >>> (i & 31)) & 1;
    this.upload(this.buffers.selection, flags);
  }

  setRows(list: RowList, rows: Uint32Array): void {
    this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, this.rowBuffers[list]);
    this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, rows, this.gl.STATIC_DRAW);
  }

  private upload(buffer: WebGLBuffer, data: ArrayBufferView): void {
    this.gl.bindBuffer(this.gl.ARRAY_BUFFER, buffer);
    this.gl.bufferData(this.gl.ARRAY_BUFFER, data, this.gl.STATIC_DRAW);
  }

  protected drawFrame(frame: Frame): void {
    const gl = this.gl;
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.viewport(0, 0, frame.size[0], frame.size[1]);
    gl.clearColor(0, 0, 0, 0);
    gl.clear(gl.COLOR_BUFFER_BIT);
    gl.bindVertexArray(this.vao);
    this.bindRows(frame);
    if (frame.dots) {
      gl.disable(gl.BLEND);
      const program = this.use(this.dotProgram, frame);
      this.draw(program, frame, 'filteredOut');
      this.draw(program, frame, 'filtered');
    }
    else {
      gl.enable(gl.BLEND);
      gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
      const program = this.use(this.shapeProgram(frame.rowTypes ? ANY_KIND : frame.markerType), frame);
      gl.uniform1i(program.uniform('u_markerCache'), +frame.markerCache);
      gl.uniform1f(program.uniform('u_border'), frame.border);
      gl.uniform1f(program.uniform('u_alpha'), frame.opacity);
      this.draw(program, frame, 'filteredOut');
      this.draw(program, frame, 'filtered');
      this.draw(program, frame, 'selected');
    }
    gl.bindVertexArray(null);
  }

  /** Polls a fence: WebGL has no promise for the GPU finishing, and a sync wait would stall the frame it times. */
  protected gpuDone(): Promise<void> {
    const gl = this.gl;
    const fence = gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0);
    if (!fence)
      return Promise.resolve();
    gl.flush();
    return new Promise((resolve) => {
      const poll = () => {
        if (gl.isContextLost() || gl.getSyncParameter(fence, gl.SYNC_STATUS) == gl.SIGNALED) {
          gl.deleteSync(fence);
          resolve();
        }
        else
          setTimeout(poll, 1);
      };
      poll();
    });
  }

  hitTest(sc: DG.ScatterPlotViewer, pt: DG.Point): number | null {
    const gl = this.gl;
    const frame = this.data.update(sc);
    // GL drops a point whose center is off the viewport, so the pick viewport reaches past the canvas by the largest
    // hit radius: a marker centered off the plot still takes the mouse on it, as on the CPU
    const margin = Math.ceil((frame.largestMarker / 2 + 2) * frame.dpr);
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
    this.bindRows(frame);
    const program = this.use(this.pickProgram, frame);
    gl.uniform2f(program.uniform('u_mouse'), pt.x, pt.y);
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.rowBuffers.filtered);
    gl.drawElements(gl.POINTS, frame.rows.filtered, gl.UNSIGNED_INT, 0);
    const pixel = new Uint8Array(4);
    gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, pixel);
    gl.bindVertexArray(null);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    return (pixel[0] | (pixel[1] << 8) | (pixel[2] << 16) | (pixel[3] << 24)) - 1;
  }

  private shapeProgram(kind: number): Program {
    let program = this.shapePrograms.get(kind);
    if (!program)
      this.shapePrograms.set(kind, program = new Program(this.gl, SHAPE_VERTEX, shapeFragment(kind)));
    return program;
  }

  /** Points the attributes at the per-row buffers, or sets the one value every row shares. */
  private bindRows(frame: Frame): void {
    const gl = this.gl;
    const attribute = (index: Attr, buffer: WebGLBuffer, perRow: boolean, size: number, type: number,
      normalized: boolean, constant: [number, number, number, number]) => {
      if (perRow) {
        gl.enableVertexAttribArray(index);
        gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
        gl.vertexAttribPointer(index, size, type, normalized, 0, 0);
      }
      else {
        gl.disableVertexAttribArray(index);
        gl.vertexAttrib4f(index, ...constant);
      }
    };
    attribute(Attr.X, this.buffers.x, true, 1, gl.FLOAT, false, [0, 0, 0, 0]);
    attribute(Attr.Y, this.buffers.y, true, 1, gl.FLOAT, false, [0, 0, 0, 0]);
    attribute(Attr.Color, this.buffers.colors, frame.rowColors, 4, gl.UNSIGNED_BYTE, true, colorBytes(frame.color));
    attribute(Attr.Size, this.buffers.sizes, frame.rowSizes, 1, gl.FLOAT, false, [frame.markerSize, 0, 0, 0]);
    attribute(Attr.Type, this.buffers.types, frame.rowTypes, 1, gl.UNSIGNED_INT, false, [frame.markerType, 0, 0, 0]);
    attribute(Attr.Selected, this.buffers.selection, frame.selection, 1, gl.UNSIGNED_BYTE, false, [0, 0, 0, 0]);
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
    const [b, g, r] = colorBytes(frame.selectedColor);
    gl.uniform3f(program.uniform('u_selectedColor'), r, g, b);
    return program;
  }

  /** Draws the rows of [list]; filtered-out rows all take the filtered-out color when there is one. */
  private draw(program: Program, frame: Frame, list: RowList): void {
    const count = frame.rows[list];
    if (count == 0)
      return;
    const gl = this.gl;
    const color = list == 'filteredOut' ? frame.filteredOutColor : null;
    if (color != null) {
      gl.disableVertexAttribArray(Attr.Color);
      gl.vertexAttrib4f(Attr.Color, ...colorBytes(color));
    }
    gl.uniform1i(program.uniform('u_selection'), +(color == null && frame.selection));
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.rowBuffers[list]);
    gl.drawElements(gl.POINTS, count, gl.UNSIGNED_INT, 0);
    if (color != null)
      this.bindRows(frame);
  }
}
