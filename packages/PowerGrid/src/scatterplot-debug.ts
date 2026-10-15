/* An opt-in timing log of the GPU scatter plot, for finding what is slow on a given machine. In the browser console:
 *   localStorage.setItem('PowerGrid.scatterPlotDebug', '1')    turns it on, no reload needed
 *   localStorage.removeItem('PowerGrid.scatterPlotDebug')      turns it off
 * Every second a plot is drawn, it prints a summary: the engine and the GPU, PowerGrid's JavaScript time per frame by
 * step, how long the GPU took to finish the frames, hit tests, how long the plot took to follow wheel and drag input,
 * and the main thread's long frames with the scripts that filled them, so the core's work shows next to the plugin's.
 * Engine set-ups, rebuilds of the per-row data, slow frames and the reasons a plot stays on the CPU print at once. */
import * as DG from 'datagrok-api/dg';
import type {Frame} from './scatterplot-data';
import {MARKER_TYPES, plotCanvases} from './scatterplot-look';

const KEY = 'PowerGrid.scatterPlotDebug';
const PREFIX = '[Scatter plot GPU]';
const SUMMARY_MS = 1000;
/** A frame that keeps PowerGrid's JavaScript busy longer than this is reported on its own. */
const SLOW_FRAME_MS = 16;
/** Input that has not moved the view for this long is not waiting for a frame any more. */
const MAX_INPUT_LAG_MS = 3000;
/** The core's work after a gesture (filter by zoom, the frame after it) lands within this of a plot's last report. */
const LONG_FRAME_GRACE_MS = 2000;

/** An entry of the Long Animation Frames API, which the TypeScript DOM library does not have yet. */
interface LongFrame extends PerformanceEntry {
  blockingDuration: number;
  scripts: {duration: number, invoker: string, invokerType: string, sourceURL: string, sourceFunctionName: string}[];
}

const fmt = (ms: number) => ms < 10 ? ms.toFixed(1) : Math.round(ms).toLocaleString();
const count = (n: number, noun: string) => `${n.toLocaleString()} ${noun}${n == 1 ? '' : 's'}`;

/** Times by median and maximum: one frame stuck behind the core's work would make an average say little. */
class Stat {
  private readonly values: number[] = [];

  get count(): number {
    return this.values.length;
  }

  add(ms: number): void {
    this.values.push(ms);
  }

  toString(): string {
    const sorted = [...this.values].sort((a, b) => a - b);
    return `${fmt(sorted[sorted.length >> 1])} median / ${fmt(sorted[sorted.length - 1])} max ms`;
  }
}

/** A script's file, or the package a published package file belongs to. */
function source(url: string): string {
  const pkg = /\/files\/([^/]+)\/.*\/([^/?]+)(\?|$)/.exec(url);
  return pkg ? `${pkg[1]}/${pkg[2]}` : url.split('/').pop()!.split('?')[0];
}

class Summary {
  readonly start = performance.now();
  frames = 0;
  readonly plots = new Set<HTMLCanvasElement>();
  engine = '';
  plot = '';
  readonly check = new Stat();
  readonly data = new Stat();
  readonly draw = new Stat();
  readonly gpu = new Stat();
  readonly hitTests = new Stat();
  readonly hitEngines = new Set<string>();
  inputs = 0;
  viewChanges = 0;
  readonly lag = new Stat();
  readonly longFrames: LongFrame[] = [];

  print(): void {
    const lines = [`${fmt((performance.now() - this.start) / 1000)} s` + (this.frames ?
      `: ${count(this.frames, 'frame')} of ${count(this.plots.size, 'plot')} with ${this.engine}; ${this.plot}` :
      ', no GPU frames')];
    if (this.frames)
      lines.push(`  PowerGrid JavaScript per frame: data ${this.data}, drawing ${this.draw}`);
    if (this.check.count) {
      lines.push(`  the core asked ${count(this.check.count, 'time')} whether the GPU draws a plot: ` +
        `${this.check} each`);
    }
    if (this.gpu.count) {
      lines.push(`  the GPU was done with the frames ${this.gpu} after they were submitted (as the main thread ` +
        'saw it, so a busy main thread adds to it)');
    }
    if (this.hitTests.count) {
      lines.push(`  hit tests: ${this.hitTests.count} on ${[...this.hitEngines].join(' and ')}, ` +
        `${this.hitTests} each`);
    }
    if (this.inputs) {
      lines.push(`  input: ${count(this.inputs, 'wheel or drag event')}; the view moved ` +
        `${count(this.viewChanges, 'time')}` + (this.lag.count ? `, ${this.lag} after the input that moved it` : ''));
    }
    if (this.longFrames.length)
      lines.push(`  main thread: ${longFramesText(this.longFrames)}`);
    console.log(`${PREFIX} ${lines.join('\n')}`);
  }
}

function longFramesText(frames: LongFrame[]): string {
  const scripts = new Map<string, number>();
  let total = 0;
  let blocking = 0;
  let inScripts = 0;
  for (const frame of frames) {
    total += frame.duration;
    blocking += frame.blockingDuration;
    for (const s of frame.scripts) {
      const file = source(s.sourceURL);
      const invoker = /^https?:/.test(s.invoker) ? source(s.invoker) : s.invoker;
      const key = `${s.invokerType} ${invoker}${s.sourceFunctionName ? ` ${s.sourceFunctionName}` : ''}` +
        (file && file != invoker ? ` (${file})` : '');
      scripts.set(key, (scripts.get(key) ?? 0) + s.duration);
      inScripts += s.duration;
    }
  }
  const top = [...scripts].sort((a, b) => b[1] - a[1]).slice(0, 4).map(([key, ms]) => `${fmt(ms)} ms ${key}`);
  return `${count(frames.length, 'long frame')}, ${fmt(total)} ms (${fmt(blocking)} ms blocking input), ` +
    `${fmt(Math.max(total - inScripts, 0))} ms of it outside scripts (style, layout, paint, waits for the GPU)` +
    (top.length ? `; scripts: ${top.join('; ')}` : '');
}

let current: Summary | null = null;
let longFrameObserver: PerformanceObserver | null = null;
let lastActivity = -Infinity;

/** The summary of this second, started on the first thing a plot reports. */
function summary(): Summary {
  lastActivity = performance.now();
  return openSummary();
}

function openSummary(): Summary {
  if (!current) {
    const s = current = new Summary();
    setTimeout(() => {
      current = null;
      s.print();
      if (!enabled()) {
        longFrameObserver?.disconnect();
        longFrameObserver = null;
      }
    }, SUMMARY_MS);
  }
  if (!longFrameObserver && PerformanceObserver.supportedEntryTypes?.includes('long-animation-frame')) {
    // long frames well after the plots went quiet belong to something else on the page
    longFrameObserver = new PerformanceObserver((list) => {
      if (enabled() && performance.now() - lastActivity < LONG_FRAME_GRACE_MS)
        openSummary().longFrames.push(...list.getEntries() as LongFrame[]);
    });
    longFrameObserver.observe({type: 'long-animation-frame'});
  }
  return current;
}

function enabled(): boolean {
  try {
    return localStorage.getItem(KEY) != null;
  }
  catch {
    return false;
  }
}

/** When the last wheel or drag input on a plot arrived that has not moved its view yet, and the view it moved. */
interface InputWatch {
  pending: number | null;
  view: string;
  remove(): void;
}

const inputWatches = new WeakMap<HTMLCanvasElement, InputWatch>();

function watchInput(canvas: HTMLCanvasElement): InputWatch {
  let watch = inputWatches.get(canvas);
  if (!watch) {
    const onInput = (e: Event) => {
      if (e.type == 'pointermove' && !(e as PointerEvent).buttons)
        return;
      w.pending ??= performance.now();
      if (enabled())
        summary().inputs++;
    };
    const w: InputWatch = watch = {pending: null, view: '', remove: () => {
      canvas.removeEventListener('wheel', onInput);
      canvas.removeEventListener('pointermove', onInput);
    }};
    canvas.addEventListener('wheel', onInput, {passive: true});
    canvas.addEventListener('pointermove', onInput, {passive: true});
    inputWatches.set(canvas, w);
  }
  return watch;
}

function markers(frame: Frame): string {
  return frame.dots ? 'dots' : frame.rowTypes ? 'markers by a column' : MARKER_TYPES[frame.markerType];
}

const lastReasons = new WeakMap<object, string | null>();

export const scatterDebug = {
  get on(): boolean {
    return enabled();
  },

  /** The start of a step worth timing, or 0 when the log is off. */
  now(): number {
    return enabled() ? performance.now() : 0;
  },

  log(message: string): void {
    if (enabled())
      console.log(`${PREFIX} ${message}`);
  },

  /** A rebuild of per-row data that ran from [start] until now. */
  rebuilt(what: string, start: number, details: string): void {
    if (start)
      this.log(`${what} rebuilt in ${fmt(performance.now() - start)} ms: ${details}`);
  },

  engineSetUp(sc: DG.ScatterPlotViewer, engine: string, start: number): void {
    this.log(`${engine} set up in ${fmt(performance.now() - start)} ms for a plot of ` +
      `${sc.dataFrame.rowCount.toLocaleString()} rows`);
  },

  /** The core's check whether the plot is drawn on the GPU; [reason] says why not, and is printed when it changes. */
  checked(sc: DG.ScatterPlotViewer, start: number, reason: string | null): void {
    if (!start)
      return;
    summary().check.add(performance.now() - start);
    if (reason !== (lastReasons.get(sc.dart) ?? null)) {
      this.log(reason ? `a plot of ${sc.dataFrame.rowCount.toLocaleString()} rows stays on the CPU: ${reason}` :
        'a plot moves back to the GPU');
    }
    lastReasons.set(sc.dart, reason);
  },

  /** A frame drawn in [start, updated] (the data) and [updated, now] (the drawing); [gpuDone] settles when the GPU is
   * through with it. */
  frame(sc: DG.ScatterPlotViewer, frame: Frame, engine: string, start: number, updated: number,
    gpuDone: Promise<void>): void {
    const drawn = performance.now();
    const s = summary();
    s.frames++;
    const canvas = plotCanvases(sc).canvas;
    s.plots.add(canvas);
    s.engine = engine;
    s.plot = `last: ${sc.dataFrame.rowCount.toLocaleString()} rows, ` +
      `${(frame.rows.filtered + frame.rows.filteredOut).toLocaleString()} drawn, ${markers(frame)}, ` +
      `${frame.size[0]}×${frame.size[1]} px at dpr ${frame.dpr}`;
    s.data.add(updated - start);
    s.draw.add(drawn - updated);
    if (drawn - start > SLOW_FRAME_MS) {
      this.log(`slow frame: ${fmt(drawn - start)} ms of PowerGrid JavaScript ` +
        `(data ${fmt(updated - start)}, drawing ${fmt(drawn - updated)})`);
    }
    gpuDone.then(() => s.gpu.add(performance.now() - drawn));

    const watch = watchInput(canvas);
    const view = frame.xform.join();
    if (view != watch.view) {
      if (watch.view)
        s.viewChanges++;
      if (watch.pending != null && drawn - watch.pending < MAX_INPUT_LAG_MS)
        s.lag.add(drawn - watch.pending);
      watch.pending = null;
      watch.view = view;
    }
  },

  hitTest(engine: string, start: number): void {
    if (!start)
      return;
    const s = summary();
    s.hitTests.add(performance.now() - start);
    s.hitEngines.add(engine);
  },

  /** Stops watching the input of a plot whose renderer is gone. */
  forget(canvas: HTMLCanvasElement): void {
    inputWatches.get(canvas)?.remove();
    inputWatches.delete(canvas);
  },
};
