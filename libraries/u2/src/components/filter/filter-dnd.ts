/* Where a dragged condition or group may land: the pure half — the rendered nodes come in as
   rects and the tree decides legality, so the shim can drive it without layout — and the gesture
   layer over it (a press on a row handle, four pixels start the drag, an indicator box shows the
   target, Escape cancels, the drop is one `Filters.move`). The designer's dnd.ts/drag.ts pair,
   without registry metadata: every group takes children, a condition splits its parent. */
import {div} from '../../core/elements.js';
import {Filters} from '../../core/filter/index.js';
import type {FilterGroup} from '../../core/filter/index.js';

const LINE = 2;
const DRAG_THRESHOLD = 4;

export interface DropRect {
  x: number;
  y: number;
  width: number;
  height: number;
}

/** A rendered node: its id and where it is, viewport coordinates. */
export interface FilterDropHit {
  id: string;
  rect: DropRect;
}

export interface FilterDropTarget {
  parentId: string;
  /** Already in `Filters.move` semantics — a same-parent move counts the node as gone. */
  index: number;
  kind: 'into' | 'line';
  /** Indicator geometry, viewport coordinates. */
  rect: DropRect;
}

interface Placed {
  parent: FilterGroup | null;
  index: number;
  depth: number;
}

function contains(r: DropRect, x: number, y: number): boolean {
  return x >= r.x && x <= r.x + r.width && y >= r.y && y <= r.y + r.height;
}

/** What a pointer at `x, y` means for the node `movingId`: the deepest hit under it is the
 * target — into a group (appended), before or after a condition in its parent (split at the
 * midpoint; advanced mode is always vertical). Refused: the node itself, its own descendants,
 * and a drop that would leave the tree as it is. */
export function resolveFilterDrop(root: FilterGroup, hits: FilterDropHit[], x: number, y: number,
  movingId: string): FilterDropTarget | null {
  const placed = new Map<string, Placed>();
  Filters.walk(root, (n, parent, index) =>
    placed.set(n.id, {parent, index, depth: parent ? placed.get(parent.id)!.depth + 1 : 0}));
  const moving = placed.get(movingId);
  if (!moving)
    return null;
  let hit: FilterDropHit | undefined;
  for (const h of hits) {
    if (placed.has(h.id) && contains(h.rect, x, y) && (!hit || placed.get(h.id)!.depth >= placed.get(hit.id)!.depth))
      hit = h;
  }
  if (!hit || hit.id === movingId)
    return null;
  const movingNode = Filters.find(root, movingId)!;
  if (Filters.isGroup(movingNode) && Filters.find(movingNode, hit.id))
    return null;
  const node = Filters.find(root, hit.id)!;
  if (Filters.isGroup(node))
    return target(moving, node.id, node.nodes.length, 'into', hit.rect);
  const {parent, index} = placed.get(node.id)!;
  const after = y > hit.rect.y + hit.rect.height / 2;
  const line = {x: hit.rect.x, y: after ? hit.rect.y + hit.rect.height : hit.rect.y,
    width: hit.rect.width, height: LINE};
  return target(moving, parent!.id, index + (after ? 1 : 0), 'line', line);
}

/** The index a move needs: the node is spliced out first, so an insertion point past it shifts
 * down by one — and landing back where it already is means nothing to do. */
function target(moving: Placed, parentId: string, index: number, kind: 'into' | 'line',
  rect: DropRect): FilterDropTarget | null {
  if (moving.parent?.id !== parentId)
    return {parentId, index, kind, rect};
  const to = index > moving.index ? index - 1 : index;
  return to === moving.index ? null : {parentId, index: to, kind, rect};
}

export interface FilterDragHost {
  /** Positioned; hosts the indicator and hears the handle presses of everything under it. */
  root: HTMLElement;
  tree(): FilterGroup;
  hits(): FilterDropHit[];
  /** Whether `parentId` takes the node `movingId` right now — a locked group takes nothing, a
   * locked node stays, a fixed shape (no adds) keeps it under its parent. */
  accepts(movingId: string, parentId: string): boolean;
  drop(id: string, target: FilterDropTarget): void;
}

interface Drag {
  id: string;
  row: HTMLElement;
  x: number;
  y: number;
  active: boolean;
  /** Where the nodes were when the drag went active: the layout holds still under a drag. */
  hits: FilterDropHit[];
  target: FilterDropTarget | null;
}

/** The pointer gesture on `[data-u2-part="handle"]`: the handle captures the pointer, so the rows
 * underneath never see the drag; the move/up listeners live on the document, so a row disposed
 * mid-drag (another builder on the same signal) still ends the gesture. */
export class FilterDragLayer {
  readonly indicator: HTMLElement;
  private _drag: Drag | undefined;
  private readonly _onDown = (e: Event) => this._down(e as PointerEvent);
  private readonly _onMove = (e: Event) => this._move(e as PointerEvent);
  private readonly _onUp = () => this._end(true);
  private readonly _onCancel = () => this._end(false);
  private readonly _onKey = (e: Event) => {
    if ((e as KeyboardEvent).key === 'Escape')
      this._end(false);
  };

  constructor(private readonly _host: FilterDragHost) {
    this.indicator = div([], 'u2-fb-drop');
    this.indicator.hidden = true;
    _host.root.append(this.indicator);
    _host.root.addEventListener('pointerdown', this._onDown);
  }

  dispose(): void {
    this._end(false);
    this._host.root.removeEventListener('pointerdown', this._onDown);
    this.indicator.remove();
  }

  private _down(e: PointerEvent): void {
    if (e.button !== 0 || this._drag)
      return;
    const handle = (e.target as Element).closest('[data-u2-part="handle"]') as HTMLElement | null;
    const row = handle?.closest('[data-u2-node]') as HTMLElement | null;
    if (!handle || !row)
      return;
    // without this the press starts a text selection that fights the drag
    e.preventDefault();
    this._drag = {id: row.dataset.u2Node!, row, x: e.clientX, y: e.clientY, active: false, hits: [], target: null};
    handle.setPointerCapture(e.pointerId);
    document.addEventListener('pointermove', this._onMove);
    document.addEventListener('pointerup', this._onUp);
    document.addEventListener('pointercancel', this._onCancel);
    document.addEventListener('keydown', this._onKey);
  }

  private _move(e: PointerEvent): void {
    const drag = this._drag;
    if (!drag)
      return;
    if (!drag.active) {
      if (Math.abs(e.clientX - drag.x) < DRAG_THRESHOLD && Math.abs(e.clientY - drag.y) < DRAG_THRESHOLD)
        return;
      drag.active = true;
      drag.hits = this._host.hits();
      drag.row.classList.add('u2-fb-dragging');
    }
    const target = resolveFilterDrop(this._host.tree(), drag.hits, e.clientX, e.clientY, drag.id);
    drag.target = target && this._host.accepts(drag.id, target.parentId) ? target : null;
    this._show(drag.target);
  }

  private _end(commit: boolean): void {
    const drag = this._drag;
    if (!drag)
      return;
    this._drag = undefined;
    document.removeEventListener('pointermove', this._onMove);
    document.removeEventListener('pointerup', this._onUp);
    document.removeEventListener('pointercancel', this._onCancel);
    document.removeEventListener('keydown', this._onKey);
    if (!drag.active)
      return;
    drag.row.classList.remove('u2-fb-dragging');
    this._show(null);
    if (commit && drag.target)
      this._host.drop(drag.id, drag.target);
  }

  /** Viewport geometry laid over the host, which scrolls under it. */
  private _show(target: FilterDropTarget | null): void {
    const box = this.indicator;
    box.hidden = target === null;
    if (!target)
      return;
    const root = this._host.root;
    const base = root.getBoundingClientRect();
    box.classList.toggle('u2-fb-drop-into', target.kind === 'into');
    box.classList.toggle('u2-fb-drop-line', target.kind === 'line');
    box.style.left = `${target.rect.x - base.left + root.scrollLeft}px`;
    box.style.top = `${target.rect.y - base.top + root.scrollTop}px`;
    box.style.width = `${target.rect.width}px`;
    box.style.height = `${target.rect.height}px`;
  }
}
