/** Data model + condition/target helpers for the in-app guide system
 *  (interactive tutorials and how-to answers).
 *
 *  A `Guide` is an ordered list of `GuideStep`s. Each step highlights a UI
 *  element (via `ui.hints`) and waits for the user to actually perform an
 *  action — click something, type a value, add/connect a node — before
 *  advancing. Targets are addressed by the `data-testid` / `data-*` attributes
 *  stamped across the UI (see utils/test-ids.ts), which is what makes the
 *  highlight-and-wait approach reliable.
 *
 *  This file is DOM-light and side-effect-free except for the poll/listen
 *  helpers, so the conditions can be unit-tested. */

import {tid} from '../utils/test-ids';
import type {FlowEditor} from '../rete/flow-editor';

/** What the guide runner needs from its host view. */
export interface GuideHost {
  /** The live editor (may be undefined very early during view init). */
  getFlow(): FlowEditor | undefined;
  /** Ensure the function browser (toolbox) is visible. */
  showFunctionBrowser(): void;
  /** An always-present element intro/outro popups can anchor to (the help button). */
  readonly anchorEl: HTMLElement;
}

export interface GuideContext {
  host: GuideHost;
  /** Aborted when the user exits the guide — every wait must respect it. */
  signal: AbortSignal;
}

export type GuideKind = 'tutorial' | 'question';

export interface GuideStep {
  /** Short heading shown in the instruction popup. */
  title: string;
  /** Instruction body (plain text). */
  text: string;
  /** The element the popup anchors to (and, by default, the one highlighted).
   *  Re-resolved continuously; return null to center the popup. */
  target?: (ctx: GuideContext) => HTMLElement | null;
  /** Elements to highlight (pulse + tint). Defaults to `[target]`. Use this to
   *  highlight more than one thing — e.g. both pins the user must connect. */
  highlights?: (ctx: GuideContext) => Array<HTMLElement | null>;
  /** Where the popup appears relative to the target. Default 'right'. */
  position?: 'top' | 'bottom' | 'left' | 'right';
  /** Optional setup run before the step shows (e.g. open the toolbox). */
  setup?: (ctx: GuideContext) => void | Promise<void>;
  /** Resolves when the step's action is done. Omit for a manual "Next" step. */
  until?: (ctx: GuideContext) => Promise<void>;
}

export interface Guide {
  id: string;
  kind: GuideKind;
  title: string;
  /** One-line summary shown in the launcher menu. */
  summary: string;
  steps: GuideStep[];
}

// ---------- target resolvers ----------

/** Find an element by raw CSS selector (document-wide). */
export function el(selector: string): HTMLElement | null {
  return document.querySelector(selector);
}

/** Target by data-testid value (already ff-namespaced — pass the parts). */
export function byTid(...parts: Array<string | number>): (ctx: GuideContext) => HTMLElement | null {
  const sel = `[data-testid="${tid(...parts)}"]`;
  return () => el(sel);
}

/** Target by an arbitrary selector. */
export function bySel(selector: string): (ctx: GuideContext) => HTMLElement | null {
  return () => el(selector);
}

/** Target a canvas node by the DG function it runs (e.g. 'OpenFile'). Matches on
 *  the node's `data-func` attribute. */
export function byNodeFunc(funcName: string): (ctx: GuideContext) => HTMLElement | null {
  const lc = funcName.toLowerCase();
  return () => (Array.from(document.querySelectorAll('.ff-node')) as HTMLElement[])
    .find((n) => (n.dataset.func ?? '').toLowerCase() === lc) ?? null;
}

/** Target a canvas node by its registered type name (`data-node-type-name`,
 *  e.g. 'Outputs/Table Output' — for built-ins that have no `data-func`). */
export function byNodeType(typeName: string): (ctx: GuideContext) => HTMLElement | null {
  return () => (Array.from(document.querySelectorAll('.ff-node')) as HTMLElement[])
    .find((n) => n.dataset.nodeTypeName === typeName) ?? null;
}

/** Target a function-browser item by the DG function it adds (`data-func`). */
export function byBrowserFunc(funcName: string): (ctx: GuideContext) => HTMLElement | null {
  const lc = funcName.toLowerCase();
  return () => (Array.from(document.querySelectorAll('[data-testid^="ff-browser-item"]')) as HTMLElement[])
    .find((it) => (it.dataset.func ?? '').toLowerCase() === lc) ?? null;
}

/** Target a specific socket dot inside the node found by `nodeResolver`.
 *  `key` is the raw socket key (e.g. 'result', 'table', 'table__pt'). */
export function socketOf(
  nodeResolver: (ctx: GuideContext) => HTMLElement | null, side: 'input' | 'output', key: string,
): (ctx: GuideContext) => HTMLElement | null {
  const part = side === 'input' ? 'socket-input' : 'socket-output';
  const sel = `[data-testid="${tid(part, key)}"]`;
  return (ctx) => {
    const node = nodeResolver(ctx);
    return (node?.querySelector(sel) as HTMLElement | null) ?? null;
  };
}

/** Target a specific property-panel input row by its raw parameter name
 *  (`data-param`, e.g. 'fullPath', 'expression', 'name'). */
export function byParam(paramName: string): (ctx: GuideContext) => HTMLElement | null {
  return () => el(`[data-param="${cssEscape(paramName)}"]`);
}

/** The editable field (textarea/input) selector for a named property input. */
export function paramFieldSelector(paramName: string): string {
  return `[data-param="${cssEscape(paramName)}"] textarea, [data-param="${cssEscape(paramName)}"] input`;
}

// ---------- low-level waits ----------

class AbortedError extends Error {
  constructor() {
    super('guide-aborted');
    this.name = 'AbortedError';
  }
}
export function isAborted(e: unknown): boolean {
  return e instanceof AbortedError;
}

/** Resolve when `pred()` becomes true; reject (AbortedError) if the signal fires.
 *  Polls on a timer — robust for state that has no dedicated event. */
export function poll(pred: () => boolean, signal: AbortSignal, intervalMs = 300): Promise<void> {
  return new Promise<void>((resolve, reject) => {
    if (signal.aborted) return reject(new AbortedError());
    if (pred()) return resolve();
    const timer = setInterval(() => {
      if (signal.aborted) {
        cleanup();
        reject(new AbortedError());
      } else if (pred()) {
        cleanup();
        resolve();
      }
    }, intervalMs);
    const onAbort = (): void => {
      cleanup();
      reject(new AbortedError());
    };
    const cleanup = (): void => {
      clearInterval(timer);
      signal.removeEventListener('abort', onAbort);
    };
    signal.addEventListener('abort', onAbort);
  });
}

/** Resolve on the next click of the element returned by `getEl` (waiting for it
 *  to appear first if necessary). */
export function waitForClick(
  getEl: (ctx: GuideContext) => HTMLElement | null, ctx: GuideContext,
): Promise<void> {
  return new Promise<void>((resolve, reject) => {
    const {signal} = ctx;
    let attached: HTMLElement | null = null;
    const onClick = (): void => {
      cleanup();
      resolve();
    };
    const onAbort = (): void => {
      cleanup();
      reject(new AbortedError());
    };
    const timer = setInterval(tryAttach, 200);
    function tryAttach(): void {
      if (attached) return;
      const node = getEl(ctx);
      if (node) {
        attached = node;
        node.addEventListener('click', onClick, {once: true});
      }
    }
    function cleanup(): void {
      clearInterval(timer);
      attached?.removeEventListener('click', onClick);
      signal.removeEventListener('abort', onAbort);
    }
    if (signal.aborted) return reject(new AbortedError());
    signal.addEventListener('abort', onAbort);
    tryAttach();
  });
}

// ---------- high-level conditions (compose into `until`) ----------

/** Wait for the user to click the step's highlighted target. */
export function untilClick(getEl: (ctx: GuideContext) => HTMLElement | null) {
  return (ctx: GuideContext): Promise<void> => waitForClick(getEl, ctx);
}

/** Wait until a node of the given registered type name (e.g. 'Inputs/Table Input')
 *  exists on the canvas. */
export function untilNodeType(typeName: string) {
  return (ctx: GuideContext): Promise<void> =>
    poll(() => (ctx.host.getFlow()?.getNodes() ?? []).some((n) => n.dgTypeName === typeName), ctx.signal);
}

/** Wait until a node whose underlying DG function name matches exists. */
export function untilFuncNode(funcName: string) {
  const lc = funcName.toLowerCase();
  return (ctx: GuideContext): Promise<void> =>
    poll(() => (ctx.host.getFlow()?.getNodes() ?? [])
      .some((n) => (n.dgFuncName ?? '').toLowerCase().includes(lc)), ctx.signal);
}

/** Wait until the node count grows beyond its value when the step started. */
export function untilMoreNodes() {
  return (ctx: GuideContext): Promise<void> => {
    const base = ctx.host.getFlow()?.getNodeCount() ?? 0;
    return poll(() => (ctx.host.getFlow()?.getNodeCount() ?? 0) > base, ctx.signal);
  };
}

/** Wait until at least one more connection exists than when the step started. */
export function untilMoreConnections() {
  return (ctx: GuideContext): Promise<void> => {
    const base = ctx.host.getFlow()?.getConnectionCount() ?? 0;
    return poll(() => (ctx.host.getFlow()?.getConnectionCount() ?? 0) > base, ctx.signal);
  };
}

/** Wait until the canvas has at least `n` nodes (absolute threshold). */
export function untilNodeCountAtLeast(n: number) {
  return (ctx: GuideContext): Promise<void> =>
    poll(() => (ctx.host.getFlow()?.getNodeCount() ?? 0) >= n, ctx.signal);
}

/** Wait until fewer nodes than at step start (a delete happened). */
export function untilFewerNodes() {
  return (ctx: GuideContext): Promise<void> => {
    const base = ctx.host.getFlow()?.getNodeCount() ?? 0;
    return poll(() => (ctx.host.getFlow()?.getNodeCount() ?? 0) < base, ctx.signal);
  };
}

/** Wait until a browser section (by its `data-section`) is expanded. */
export function untilSectionExpanded(title: string) {
  return (ctx: GuideContext): Promise<void> =>
    poll(() => {
      const header = el(`.funcflow-section-header[data-section="${cssEscape(title)}"]`);
      return !!header && !header.classList.contains('collapsed');
    }, ctx.signal);
}

/** Wait until an element matching the selector appears in the DOM. */
export function untilExists(selector: string) {
  return (ctx: GuideContext): Promise<void> => poll(() => !!el(selector), ctx.signal);
}

/** Wait until the input/select matching `selector` contains `substr` (case-insensitive). */
export function untilValueContains(selector: string, substr: string) {
  const needle = substr.toLowerCase();
  return (ctx: GuideContext): Promise<void> =>
    poll(() => {
      const node = el(selector) as HTMLInputElement | null;
      return !!node && (node.value ?? '').toLowerCase().includes(needle);
    }, ctx.signal);
}

/** Wait until a node is selected (the property panel shows a node's title row). */
export function untilNodeSelected() {
  return (ctx: GuideContext): Promise<void> =>
    untilExists(`[data-testid="${tid('property-title-row')}"]`)(ctx);
}

/** Wait until a node on the canvas is collapsed. */
export function untilNodeCollapsed() {
  return (ctx: GuideContext): Promise<void> =>
    poll(() => !!el('.ff-node.ff-node-collapsed'), ctx.signal);
}

/** Wait until a node running the given DG function is selected (its settings
 *  panel opens). Keys off the node's `data-selected`/`data-func` attributes. */
export function untilNodeSelectedOfFunc(funcName: string) {
  const lc = funcName.toLowerCase();
  return (ctx: GuideContext): Promise<void> =>
    poll(() => (Array.from(document.querySelectorAll('.ff-node[data-selected="true"]')) as HTMLElement[])
      .some((n) => (n.dataset.func ?? '').toLowerCase().includes(lc)), ctx.signal);
}

/** Wait until the input matching `selector` contains `term`, ignoring case AND
 *  whitespace — so "Open File", "open file", and "openfile" all satisfy
 *  `untilValueMatches(sel, 'openfile')`. */
export function untilValueMatches(selector: string, term: string) {
  const needle = term.toLowerCase().replace(/\s+/g, '');
  return (ctx: GuideContext): Promise<void> =>
    poll(() => {
      const node = el(selector) as HTMLInputElement | null;
      return !!node && (node.value ?? '').toLowerCase().replace(/\s+/g, '').includes(needle);
    }, ctx.signal);
}

/** Wait until the node found by `rightResolver` sits at least `minDx` screen-px
 *  to the right of the node found by `leftResolver` (user dragged it clear). */
export function untilNodeRightOf(
  rightResolver: (ctx: GuideContext) => HTMLElement | null,
  leftResolver: (ctx: GuideContext) => HTMLElement | null,
  minDx = 200,
) {
  return (ctx: GuideContext): Promise<void> =>
    poll(() => {
      const a = rightResolver(ctx);
      const b = leftResolver(ctx);
      if (!a || !b) return false;
      return a.getBoundingClientRect().left - b.getBoundingClientRect().left >= minDx;
    }, ctx.signal);
}

/** Wait until the input/textarea matching `selector` has a non-empty value. */
export function untilValueNonEmpty(selector: string) {
  return (ctx: GuideContext): Promise<void> =>
    poll(() => {
      const node = el(selector) as HTMLInputElement | null;
      return !!node && (node.value ?? '').trim().length > 0;
    }, ctx.signal);
}

/** Best-effort copy to the clipboard so a step can say "paste it". Returns
 *  whether it succeeded — callers should still show the literal text as a
 *  fallback for when clipboard access is denied. */
export async function copyToClipboard(text: string): Promise<boolean> {
  try {
    await navigator.clipboard.writeText(text);
    return true;
  } catch {
    return false;
  }
}

/** Prefill the function-browser search box (and fire its filter) so a single
 *  named function is shown, ready to double-click. DOM-only — no host API. */
export function prefillSearch(text: string): void {
  const input = el('[data-testid="ff-browser-search"]') as HTMLInputElement | null;
  if (!input) return;
  input.value = text;
  input.dispatchEvent(new Event('input', {bubbles: true}));
}

/** Minimal CSS attribute-value escaper (titles here are simple ASCII). */
function cssEscape(s: string): string {
  return s.replace(/["\\]/g, '\\$&');
}

// ---------- popup placement (pure, unit-tested) ----------

export type Side = 'top' | 'bottom' | 'left' | 'right';
export interface PlaceRect {left: number; top: number; right: number; bottom: number; width: number; height: number}
export interface Placement {side: Side | 'center'; x: number; y: number}

const OPPOSITE: Record<Side, Side> = {right: 'left', left: 'right', top: 'bottom', bottom: 'top'};

/** Choose where a popup of size pw×ph goes next to `target` within a vw×vh
 *  viewport: honor `preferred` if it fits, else its opposite, else
 *  right→left→bottom→top; then clamp fully on screen. With no target, center
 *  horizontally / upper-third vertically. Pure so it can be unit-tested. */
export function computePlacement(
  target: PlaceRect | null, pw: number, ph: number, vw: number, vh: number,
  preferred?: Side, gap = 14, margin = 10,
): Placement {
  if (!target)
    return {side: 'center', x: Math.round((vw - pw) / 2), y: Math.round((vh - ph) / 3)};

  const cx = target.left + target.width / 2;
  const cy = target.top + target.height / 2;
  const at: Record<Side, {x: number; y: number}> = {
    right: {x: target.right + gap, y: cy - ph / 2},
    left: {x: target.left - gap - pw, y: cy - ph / 2},
    bottom: {x: cx - pw / 2, y: target.bottom + gap},
    top: {x: cx - pw / 2, y: target.top - gap - ph},
  };
  const fits = (s: Side): boolean => {
    if (s === 'right') return at.right.x + pw <= vw - margin;
    if (s === 'left') return at.left.x >= margin;
    if (s === 'bottom') return at.bottom.y + ph <= vh - margin;
    return at.top.y >= margin;
  };

  const order: Side[] = [];
  const add = (s: Side): void => {
    if (!order.includes(s)) order.push(s);
  };
  if (preferred) {
    add(preferred);
    add(OPPOSITE[preferred]);
  }
  (['right', 'left', 'bottom', 'top'] as Side[]).forEach(add);

  let side: Side = order[0];
  for (const s of order) {
    if (fits(s)) {
      side = s;
      break;
    }
  }
  const x = Math.max(margin, Math.min(at[side].x, vw - pw - margin));
  const y = Math.max(margin, Math.min(at[side].y, vh - ph - margin));
  return {side, x: Math.round(x), y: Math.round(y)};
}
