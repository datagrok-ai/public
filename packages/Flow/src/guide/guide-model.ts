/** Data model + condition/target helpers for the in-app guide system.
 *  DOM-light and side-effect-free except the poll/listen helpers, so conditions are unit-testable. */

import * as DG from 'datagrok-api/dg';
import {tid} from '../utils/test-ids';
import type {FlowEditor} from '../rete/flow-editor';

/** What the guide runner needs from its host view. */
export interface GuideHost {
  /** May be undefined very early during view init. */
  getFlow(): FlowEditor | undefined;
  showFunctionBrowser(): void;
  showToolboxTab(name: 'Files' | 'Spaces' | 'Queries' | 'Workflows' | 'Favorites'): void;
  /** Optional so bare test hosts don't need it. */
  hideStartPanel?(): void;
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
  title: string;
  text: string;
  /** Re-resolved continuously; return null to center the popup. */
  target?: (ctx: GuideContext) => HTMLElement | null;
  /** Defaults to `[target]` — set to highlight several things (e.g. both pins to connect). */
  highlights?: (ctx: GuideContext) => Array<HTMLElement | null>;
  /** Extra elements the card must not cover (highlights are avoided automatically). */
  avoid?: (ctx: GuideContext) => Array<HTMLElement | null>;
  /** Default 'right'. */
  position?: 'top' | 'bottom' | 'left' | 'right';
  setup?: (ctx: GuideContext) => void | Promise<void>;
  /** Skipped entirely at step start when true — for already-satisfied prerequisites. */
  skipIf?: (ctx: GuideContext) => boolean;
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

export function el(selector: string): HTMLElement | null {
  return document.querySelector(selector);
}

/** Target by data-testid value — pass the parts, `tid` namespaces them. */
export function byTid(...parts: Array<string | number>): (ctx: GuideContext) => HTMLElement | null {
  const sel = `[data-testid="${tid(...parts)}"]`;
  return () => el(sel);
}

export function bySel(selector: string): (ctx: GuideContext) => HTMLElement | null {
  return () => el(selector);
}

/** Target a canvas node by its `data-func` (e.g. 'OpenFile'). */
export function byNodeFunc(funcName: string): (ctx: GuideContext) => HTMLElement | null {
  const lc = funcName.toLowerCase();
  return () => (Array.from(document.querySelectorAll('.ff-node')) as HTMLElement[])
    .find((n) => (n.dataset.func ?? '').toLowerCase() === lc) ?? null;
}

/** Target a canvas node by its registered type name — for built-ins with no `data-func`. */
export function byNodeType(typeName: string): (ctx: GuideContext) => HTMLElement | null {
  return () => (Array.from(document.querySelectorAll('.ff-node')) as HTMLElement[])
    .find((n) => n.dataset.nodeTypeName === typeName) ?? null;
}

/** The `index`-th node running a function (DOM order ≈ creation order) — disambiguates identical nodes. */
export function byNodeFuncNth(funcName: string, index: number): (ctx: GuideContext) => HTMLElement | null {
  const lc = funcName.toLowerCase();
  return () => (Array.from(document.querySelectorAll('.ff-node')) as HTMLElement[])
    .filter((n) => (n.dataset.func ?? '').toLowerCase() === lc)[index] ?? null;
}

/** Top-most open dialog — the card sits above dialogs, so a step opening one must re-anchor to it. */
export function openDialogEl(): HTMLElement | null {
  const dialogs = DG.Dialog.getOpenDialogs();
  return dialogs.length ? dialogs[dialogs.length - 1].root : null;
}

/** Anchors to an open dialog when there is one, else to `resolver`. */
export function preferDialog(
  resolver: (ctx: GuideContext) => HTMLElement | null,
): (ctx: GuideContext) => HTMLElement | null {
  return (ctx) => openDialogEl() ?? resolver(ctx);
}

/** Target a function-browser item by its `data-func`. */
export function byBrowserFunc(funcName: string): (ctx: GuideContext) => HTMLElement | null {
  const lc = funcName.toLowerCase();
  return () => (Array.from(document.querySelectorAll('[data-testid^="ff-browser-item"]')) as HTMLElement[])
    .find((it) => (it.dataset.func ?? '').toLowerCase() === lc) ?? null;
}

/** A socket dot inside the resolved node; `key` is the raw socket key (e.g. 'table__pt'). */
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

/** A property-panel input row by its raw parameter name (`data-param`). */
export function byParam(paramName: string): (ctx: GuideContext) => HTMLElement | null {
  return () => el(`[data-param="${cssEscape(paramName)}"]`);
}

/** The editable field (textarea/input) selector for a named property input. */
export function paramFieldSelector(paramName: string): string {
  return `[data-param="${cssEscape(paramName)}"] textarea, [data-param="${cssEscape(paramName)}"] input`;
}

export function byFileTreeConn(connName: string): (ctx: GuideContext) => HTMLElement | null {
  return byTid('files-conn', connName);
}

export function byFileTreeFile(fileName: string): (ctx: GuideContext) => HTMLElement | null {
  return byTid('files-file', fileName);
}

class AbortedError extends Error {
  constructor() {
    super('guide-aborted');
    this.name = 'AbortedError';
  }
}
export function isAborted(e: unknown): boolean {
  return e instanceof AbortedError;
}

/** Polls on a timer — robust for state that has no dedicated event. */
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

/** Resolves on the next click of the element, waiting for it to appear first if necessary. */
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

export function untilClick(getEl: (ctx: GuideContext) => HTMLElement | null) {
  return (ctx: GuideContext): Promise<void> => waitForClick(getEl, ctx);
}

export function untilNodeType(typeName: string) {
  return (ctx: GuideContext): Promise<void> =>
    poll(() => (ctx.host.getFlow()?.getNodes() ?? []).some((n) => n.dgTypeName === typeName), ctx.signal);
}

export function untilFuncNode(funcName: string) {
  const lc = funcName.toLowerCase();
  return (ctx: GuideContext): Promise<void> =>
    poll(() => (ctx.host.getFlow()?.getNodes() ?? [])
      .some((n) => (n.dgFuncName ?? '').toLowerCase().includes(lc)), ctx.signal);
}

export function untilMoreNodes() {
  return (ctx: GuideContext): Promise<void> => {
    const base = ctx.host.getFlow()?.getNodeCount() ?? 0;
    return poll(() => (ctx.host.getFlow()?.getNodeCount() ?? 0) > base, ctx.signal);
  };
}

export function untilMoreConnections() {
  return (ctx: GuideContext): Promise<void> => {
    const base = ctx.host.getFlow()?.getConnectionCount() ?? 0;
    return poll(() => (ctx.host.getFlow()?.getConnectionCount() ?? 0) > base, ctx.signal);
  };
}

export function untilFewerNodes() {
  return (ctx: GuideContext): Promise<void> => {
    const base = ctx.host.getFlow()?.getNodeCount() ?? 0;
    return poll(() => (ctx.host.getFlow()?.getNodeCount() ?? 0) < base, ctx.signal);
  };
}

/** An expanded accordion pane's header carries the `expanded` class (set by DG.AccordionPane). */
export function untilSectionExpanded(title: string) {
  return (ctx: GuideContext): Promise<void> =>
    poll(() => {
      const header = el(`.d4-accordion-pane-header[data-section="${cssEscape(title)}"]`);
      return !!header && header.classList.contains('expanded');
    }, ctx.signal);
}

export function untilExists(selector: string) {
  return (ctx: GuideContext): Promise<void> => poll(() => !!el(selector), ctx.signal);
}

export function untilValueContains(selector: string, substr: string) {
  const needle = substr.toLowerCase();
  return (ctx: GuideContext): Promise<void> =>
    poll(() => {
      const node = el(selector) as HTMLInputElement | null;
      return !!node && (node.value ?? '').toLowerCase().includes(needle);
    }, ctx.signal);
}

/** NOTE: satisfied by a stale panel — prefer untilNodeOfTypeSelected / untilNodeSelectedOfFunc. */
export function untilNodeSelected() {
  return (ctx: GuideContext): Promise<void> =>
    untilExists(`[data-testid="${tid('property-title-row')}"]`)(ctx);
}

/** Keys off the canvas node's `data-selected` — immune to a stale property panel. */
export function untilNodeOfTypeSelected(typeName: string) {
  return (ctx: GuideContext): Promise<void> =>
    poll(() => (Array.from(document.querySelectorAll('.ff-node[data-selected="true"]')) as HTMLElement[])
      .some((n) => n.dataset.nodeTypeName === typeName), ctx.signal);
}

/** Detects the user's collapse action, not a node already collapsed at step start. */
export function untilMoreCollapsed() {
  return (ctx: GuideContext): Promise<void> => {
    const base = document.querySelectorAll('.ff-node.ff-node-collapsed').length;
    return poll(() => document.querySelectorAll('.ff-node.ff-node-collapsed').length > base, ctx.signal);
  };
}

export function nodeCount(ctx: GuideContext): number {
  return ctx.host.getFlow()?.getNodeCount() ?? 0;
}

export function hasFuncNode(funcName: string): (ctx: GuideContext) => boolean {
  const lc = funcName.toLowerCase();
  return (ctx) => (ctx.host.getFlow()?.getNodes() ?? [])
    .some((n) => (n.dgFuncName ?? '').toLowerCase().includes(lc));
}

export function hasNodeType(typeName: string): (ctx: GuideContext) => boolean {
  return (ctx) => (ctx.host.getFlow()?.getNodes() ?? []).some((n) => n.dgTypeName === typeName);
}

export function untilNodeSelectedOfFunc(funcName: string) {
  const lc = funcName.toLowerCase();
  return (ctx: GuideContext): Promise<void> =>
    poll(() => (Array.from(document.querySelectorAll('.ff-node[data-selected="true"]')) as HTMLElement[])
      .some((n) => (n.dataset.func ?? '').toLowerCase().includes(lc)), ctx.signal);
}

/** Ignores case AND whitespace — "Open File" satisfies `untilValueMatches(sel, 'openfile')`. */
export function untilValueMatches(selector: string, term: string) {
  const needle = term.toLowerCase().replace(/\s+/g, '');
  return (ctx: GuideContext): Promise<void> =>
    poll(() => {
      const node = el(selector) as HTMLInputElement | null;
      return !!node && (node.value ?? '').toLowerCase().replace(/\s+/g, '').includes(needle);
    }, ctx.signal);
}

/** Doubles as a `skipIf` so "drag it clear" steps skip when the node already sits clear. */
export function nodeIsRightOf(
  rightResolver: (ctx: GuideContext) => HTMLElement | null,
  leftResolver: (ctx: GuideContext) => HTMLElement | null,
  minDx = 200,
): (ctx: GuideContext) => boolean {
  return (ctx) => {
    const a = rightResolver(ctx);
    const b = leftResolver(ctx);
    if (!a || !b) return false;
    return a.getBoundingClientRect().left - b.getBoundingClientRect().left >= minDx;
  };
}

export function untilNodeRightOf(
  rightResolver: (ctx: GuideContext) => HTMLElement | null,
  leftResolver: (ctx: GuideContext) => HTMLElement | null,
  minDx = 200,
) {
  return (ctx: GuideContext): Promise<void> =>
    poll(() => nodeIsRightOf(rightResolver, leftResolver, minDx)(ctx), ctx.signal);
}

/** Freshly added nodes can land half-overlapping, hiding the dots a "connect" step highlights. */
export function nodeIsApart(
  resolver: (ctx: GuideContext) => HTMLElement | null, margin = 8,
): (ctx: GuideContext) => boolean {
  return (ctx) => {
    const el = resolver(ctx);
    if (!el) return false;
    const r = el.getBoundingClientRect();
    if (r.width === 0 || r.height === 0) return false;
    for (const other of Array.from(document.querySelectorAll('.ff-node')) as HTMLElement[]) {
      if (other === el || el.contains(other) || other.contains(el)) continue;
      const o = other.getBoundingClientRect();
      if (o.width === 0) continue;
      if (r.left < o.right + margin && r.right > o.left - margin &&
          r.top < o.bottom + margin && r.bottom > o.top - margin)
        return false;
    }
    return true;
  };
}

export function untilNodeApart(
  resolver: (ctx: GuideContext) => HTMLElement | null, margin = 8,
) {
  return (ctx: GuideContext): Promise<void> =>
    poll(() => nodeIsApart(resolver, margin)(ctx), ctx.signal);
}

export function untilNodeMovedBy(
  resolver: (ctx: GuideContext) => HTMLElement | null, minPx = 60,
) {
  return (ctx: GuideContext): Promise<void> => {
    const r0 = resolver(ctx)?.getBoundingClientRect();
    return poll(() => {
      const el = resolver(ctx);
      if (!el || !r0) return false;
      const r = el.getBoundingClientRect();
      return Math.hypot(r.left - r0.left, r.top - r0.top) >= minPx;
    }, ctx.signal);
  };
}

/** `untilExists` alone is fooled by present-but-hidden hosts. */
export function untilVisible(selector: string) {
  return (ctx: GuideContext): Promise<void> =>
    poll(() => Array.from(document.querySelectorAll(selector)).some((node) => {
      const r = node.getBoundingClientRect();
      return r.width > 0 && r.height > 0;
    }), ctx.signal);
}

export function untilFileTreeConnExpanded(connName: string) {
  const sel = `[data-testid="${tid('files-conn', connName)}"] .d4-tree-view-tri-expanded`;
  return (ctx: GuideContext): Promise<void> => poll(() => !!el(sel), ctx.signal);
}

export function isScrolledIntoView(node: HTMLElement): boolean {
  const r = node.getBoundingClientRect();
  if (r.height === 0 || r.bottom <= 0 || r.top >= window.innerHeight) return false;
  let p = node.parentElement;
  while (p) {
    const oy = getComputedStyle(p).overflowY;
    if (/(auto|scroll|hidden)/.test(oy) && p.scrollHeight > p.clientHeight + 1) {
      const pr = p.getBoundingClientRect();
      return r.top >= pr.top - 2 && r.bottom <= pr.bottom + 2;
    }
    p = p.parentElement;
  }
  return true; // no scroll container — viewport visibility already checked
}

export function untilScrolledIntoView(selector: string) {
  return (ctx: GuideContext): Promise<void> =>
    poll(() => {
      const node = el(selector);
      return !!node && isScrolledIntoView(node);
    }, ctx.signal);
}

export function untilFuncNodeWithInput(funcName: string, inputKey: string, substr: string) {
  const lcFn = funcName.toLowerCase();
  const lcSub = substr.toLowerCase();
  return (ctx: GuideContext): Promise<void> =>
    poll(() => (ctx.host.getFlow()?.getNodes() ?? []).some((n) =>
      (n.dgFuncName ?? '').toLowerCase().includes(lcFn) &&
      String((n.inputValues ?? {})[inputKey] ?? '').toLowerCase().includes(lcSub)), ctx.signal);
}

/** The user picked at least `n` columns (a comma-separated `column_list` field). */
export function untilColumnCountAtLeast(selector: string, n: number) {
  return (ctx: GuideContext): Promise<void> =>
    poll(() => {
      const node = el(selector) as HTMLInputElement | null;
      if (!node) return false;
      const tokens = (node.value ?? '').split(',').map((s) => s.trim()).filter(Boolean);
      return tokens.length >= n;
    }, ctx.signal);
}

export function untilValueNonEmpty(selector: string) {
  return (ctx: GuideContext): Promise<void> =>
    poll(() => {
      const node = el(selector) as HTMLInputElement | null;
      return !!node && (node.value ?? '').trim().length > 0;
    }, ctx.signal);
}

/** Returns whether it succeeded — callers should show the literal text as a fallback. */
export async function copyToClipboard(text: string): Promise<boolean> {
  try {
    await navigator.clipboard.writeText(text);
    return true;
  } catch {
    return false;
  }
}

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

export type Side = 'top' | 'bottom' | 'left' | 'right';
export interface PlaceRect {left: number; top: number; right: number; bottom: number; width: number; height: number}
export interface Placement {side: Side | 'center'; x: number; y: number}

const OPPOSITE: Record<Side, Side> = {right: 'left', left: 'right', top: 'bottom', bottom: 'top'};

/** Picks the best side (preferred → opposite → right/left/bottom/top), clamps on-screen,
 *  and prefers sides whose popup covers no `avoid` rect. Pure — unit-tested. */
export function computePlacement(
  target: PlaceRect | null, pw: number, ph: number, vw: number, vh: number,
  preferred?: Side, gap = 14, margin = 10, avoid: PlaceRect[] = [],
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
  const clamped = (s: Side): {x: number; y: number} => ({
    x: Math.max(margin, Math.min(at[s].x, vw - pw - margin)),
    y: Math.max(margin, Math.min(at[s].y, vh - ph - margin)),
  });
  const clear = (s: Side): boolean => {
    const c = clamped(s);
    return !avoid.some((a) =>
      c.x < a.right && c.x + pw > a.left && c.y < a.bottom && c.y + ph > a.top);
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

  const side: Side = order.find((s) => fits(s) && clear(s)) ??
    order.find(fits) ?? order[0];
  const {x, y} = clamped(side);
  return {side, x: Math.round(x), y: Math.round(y)};
}
