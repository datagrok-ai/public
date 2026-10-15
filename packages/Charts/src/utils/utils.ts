/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Observable, Subject, Subscription} from 'rxjs';
import {debounceTime, finalize, tap} from 'rxjs/operators';

export const ERROR_CLASS = 'd4-viewer-error';

/** A hit area, in CSS px of the part it is reported against. */
export type Box = {x: number, y: number, width: number, height: number};

export type Readings = {[name: string]: number | string | boolean};

export function unsubscribeAll(subs: Subscription[]): void {
  subs.forEach((sub) => sub.unsubscribe());
  subs.length = 0;
}

/** Runs [done] once [owed] stops reporting a frame outstanding, re-checking on up to [frames] frames.
 * `setOption` returns before echarts lays a series out, and zrender paints it a frame later, so
 * neither the promise nor the canvas says the picture is there — what the layout left does. */
export class LayoutSettler {
  private timer: any = null;

  constructor(private readonly frames = 20) {}

  settle(owed: () => boolean, done: () => void, attempt = 0): void {
    this.timer = null;
    if (owed() && attempt < this.frames) {
      this.timer = setTimeout(() => requestAnimationFrame(() => this.settle(owed, done, attempt + 1)));
      return;
    }
    done();
  }

  cancel(): void {
    if (this.timer !== null)
      clearTimeout(this.timer);
    this.timer = null;
  }
}

/** What automation settles a viewer on: [pending] from the moment a change arrives until the frame
 * that shows it, and [rendered] after every render. [frameOwed] says whether the picture still owes
 * a frame once the render code has run. */
export class RenderSignals {
  readonly rendered = new Subject<void>();
  private owed = 0;
  private settling = false;
  private readonly settler: LayoutSettler;

  constructor(private readonly frameOwed: () => boolean = () => false, frames = 20) {
    this.settler = new LayoutSettler(frames);
  }

  get pending(): boolean {return this.owed > 0 || this.settling;}

  hold(): () => void {
    this.owed++;
    let held = true;
    return () => {
      if (held)
        this.owed--;
      held = false;
    };
  }

  /** `DG.debounce` that holds [pending] from the first change of a burst, not from the render that
   * comes [ms] later — and lets go if the subscription ends first. */
  debounce<T>(stream: Observable<T>, ms: number): Observable<T> {
    let release: (() => void) | null = null;
    const pay = () => {
      release?.();
      release = null;
    };
    return stream.pipe(tap(() => release ??= this.hold()), debounceTime(ms), tap(pay), finalize(pay));
  }

  render(draw: () => void): void {
    const release = this.hold();
    try {
      draw();
    } finally {
      release();
      this.settle();
    }
  }

  /** Holds [pending] until an async render ends either way, then settles. The failure is logged, not
   * rethrown: a render queue chained on a rejected promise would skip every later render. */
  track(work: Promise<unknown>): Promise<void> {
    const release = this.hold();
    const settle = () => {
      release();
      this.settle();
    };
    return work.then(settle, (e) => {
      settle();
      console.error(e);
    });
  }

  /** Fires [rendered] once [frameOwed] stops reporting a frame outstanding; [pending] until then. */
  settle(): void {
    this.settling = true;
    this.settler.cancel();
    this.settler.settle(this.frameOwed, () => {
      this.settling = false;
      this.rendered.next();
    });
  }

  reset(): void {
    this.settler.cancel();
    this.settling = false;
  }
}

/** Whether echarts still owes a frame: a lazy `setOption` it has not applied, or a change zrender has
 * not painted (a dispatched action or a resize marks the elements and paints on the next frame). */
export function echartsFramePending(chart: any): boolean {
  return chart?.__pendingUpdate != null || chart?.getZr?.()?._needsRefresh === true;
}

/** An empty category value as an automation name — `''` and the `' '` a tree gives a null would
 * otherwise make `segment A | ` read as `segment A`. */
export function categoryName(name: unknown): string {
  const text = String(name ?? '');
  return text.trim() === '' ? '(empty)' : text;
}

/** The names of an echarts tree node and its ancestors, from the top down; the virtual root echarts
 * wraps the data in (the node with no parent) is not one of them. */
export function treePath(node: any): string[] {
  const names: string[] = [];
  for (; node != null && node.parentNode != null; node = node.parentNode)
    names.unshift(categoryName(node.name));
  return names;
}

export function drawnItems(data: any): number[] {
  return Array.from({length: data?.count?.() ?? 0}, (_, i) => i).filter((i) => data.getItemGraphicEl(i) != null);
}

/** A point of [el]'s user space in px of its outermost svg, which is what an svg viewer's hit
 * areas are relative to. */
export function svgPoint(el: SVGGraphicsElement, x: number, y: number): {x: number, y: number} {
  const m = el.getCTM();
  return m === null ? {x, y} : {x: m.a * x + m.c * y + m.e, y: m.b * x + m.d * y + m.f};
}

/** Whether a click at [p] (px of [svg]) lands on something [accepts] takes for its target. A point
 * the document cannot hit-test there — the viewer scrolled away or covered — is taken on its
 * geometry alone. */
export function svgHits(svg: SVGSVGElement, p: {x: number, y: number}, accepts: (hit: Element) => boolean): boolean {
  const box = svg.getBoundingClientRect();
  const hit = document.elementFromPoint(box.left + p.x, box.top + p.y);
  return hit === null || !svg.contains(hit) || accepts(hit);
}

export namespace ts {
  /** A type guard function.
   * See https://stackoverflow.com/questions/64616994/typescript-type-narrowing-not-working-for-in-when-key-is-stored-in-a-variable*/
  export function hasProp<T extends object>(obj: T, key: PropertyKey): key is keyof T {
    return key in obj;
  }
}

export namespace data {
  export function mapToRange(x: number, min1: number, max1: number, min2: number, max2: number): number {
    const range1 = max1 - min1;
    const range2 = max2 - min2;
    if (range1 === 0)
      return Math.min(Math.max(x, min2), max2);
    return (((x - min1) * range2) / range1) + min2;
  }

  export function aggToStat(dataframe: DG.DataFrame, columnName: string,
    aggregation: DG.AggregationType): number | null {
    const colStatsCall = 'dataframe.getCol(columnName).stats.';
    const stats = {
      'avg': colStatsCall + 'avg',
      'count': colStatsCall + 'totalCount',
      'kurt': colStatsCall + 'kurt',
      'max': colStatsCall + 'max',
      'med': colStatsCall + 'med',
      'min': colStatsCall + 'min',
      'nulls': colStatsCall + 'missingValueCount',
      'q1': colStatsCall + 'q1',
      'q2': colStatsCall + 'q2',
      'q3': colStatsCall + 'q3',
      'skew': colStatsCall + 'skew',
      'stdev': colStatsCall + 'stdev',
      'sum': colStatsCall + 'sum',
      'unique': colStatsCall + 'uniqueCount',
      'values': colStatsCall + 'valueCount',
      'variance': colStatsCall + 'variance',
      '#selected': 'dataframe.selection.trueCount',
      'first': 'dataframe.getCol(columnName).get(0)',
    };
    return ts.hasProp(stats, aggregation) ? eval(stats[aggregation]) : null;
  }
}

export namespace MessageHandler {
  export function _showMessage(root: HTMLElement, msg: string, className: string) {
    _removeMessage(root, className);
    root.appendChild(ui.divText(msg, className));
  }

  export function _removeMessage(root: HTMLElement, className: string) {
    const divTextElement = root.getElementsByClassName(className)[0];
    if (divTextElement)
      root.removeChild(divTextElement);
  }

  export function _getMessage(root: HTMLElement, className: string = ERROR_CLASS): string | null {
    return root.querySelector(`.${className}`)?.textContent ?? null;
  }
}
