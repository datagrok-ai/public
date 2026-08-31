/* Where a drag may land (D-7): the pure half of drag-and-drop — what a node accepts, what a
   pointer over a rendered node resolves to, and where a box over the canvas is placed for a rect.
   Legality comes from registry child metadata alone, so no component carries designer code;
   geometry comes from the rendered elements. Platform-free, like node-ref.ts — the view is the
   gesture layer over this. */
import {SpecInstance, isHtmlTag} from '../../spec/spec.js';
import type {SpecNode} from '../../spec/spec.js';
import type {Registry} from '../../spec/registry.js';
import {SpecNodeRef} from './node-ref.js';

const LINE = 2;

export interface DropRect {
  x: number;
  y: number;
  width: number;
  height: number;
}

export interface DropTarget {
  parent: SpecNode;
  /** Already in `add`/`move` patch semantics — a same-parent move counts the node as gone. */
  index: number;
  kind: 'into' | 'line';
  /** Indicator geometry, viewport coordinates. */
  rect: DropRect;
}

/** Viewport geometry laid over the canvas, which scrolls under it. */
export function place(root: HTMLElement, box: HTMLElement, rect: DropRect | DOMRect | null): void {
  if (!rect) {
    box.style.display = 'none';
    return;
  }
  const base = root.getBoundingClientRect();
  box.style.display = 'block';
  box.style.left = `${rect.x - base.left + root.scrollLeft}px`;
  box.style.top = `${rect.y - base.top + root.scrollTop}px`;
  box.style.width = `${rect.width}px`;
  box.style.height = `${rect.height}px`;
}

/** A container in the eyes of the spec renderer: it appends, adopts or consumes its children. */
export function accepts(reg: Registry, node: SpecNode): boolean {
  const meta = reg.get(node.tag);
  if (!meta)
    return isHtmlTag(node.tag);
  return meta.acceptsChildren === true || meta.adopt !== undefined || meta.createWithChildren !== undefined;
}

/** What a pointer at `x, y` over the rendered `hit` means: dropping into it when it takes children,
 * before or after it in its parent otherwise. A null `hit` — the pointer is on the canvas but over
 * no rendered node (designer padding, the area past a small form) — appends into the root when the
 * root takes children. `moving` — the node being dragged — additionally refuses itself, its own
 * descendants and drops that would change nothing. */
export function resolveDrop(instance: SpecInstance, hit: SpecNode | null, x: number, y: number,
  moving?: SpecNode): DropTarget | null {
  const reg = instance.registry;
  if (!hit) {
    const root = instance.spec.root;
    if (!accepts(reg, root))
      return null;
    const rect = rectOf(instance, root);
    return rect ? target(instance, root, (root.children ?? []).length, 'into', rect, moving) : null;
  }
  if (moving && (moving === hit || SpecInstance.contains(moving, hit)))
    return null;
  if (accepts(reg, hit)) {
    const rect = rectOf(instance, hit);
    return rect ? target(instance, hit, (hit.children ?? []).length, 'into', rect, moving) : null;
  }
  const parent = instance.parentOf(hit);
  if (!parent || !accepts(reg, parent))
    return null;
  const rect = rectOf(instance, hit);
  const parentEl = new SpecNodeRef(instance, parent).element();
  if (!rect || !parentEl)
    return null;
  const horizontal = (getComputedStyle(parentEl).flexDirection ?? '').startsWith('row');
  const after = horizontal ? x > rect.x + rect.width / 2 : y > rect.y + rect.height / 2;
  const index = (parent.children ?? []).indexOf(hit) + (after ? 1 : 0);
  const line: DropRect = horizontal ?
    {x: after ? rect.x + rect.width : rect.x, y: rect.y, width: LINE, height: rect.height} :
    {x: rect.x, y: after ? rect.y + rect.height : rect.y, width: rect.width, height: LINE};
  return target(instance, parent, index, 'line', line, moving);
}

/** The patch index a move needs: the node is spliced out first, so an insertion point past it
 * shifts down by one — and landing back where it already is means nothing to do. */
function target(instance: SpecInstance, parent: SpecNode, index: number, kind: 'into' | 'line',
  rect: DropRect, moving?: SpecNode): DropTarget | null {
  if (!moving || instance.parentOf(moving) !== parent)
    return {parent, index, kind, rect};
  const from = (parent.children ?? []).indexOf(moving);
  const to = index > from ? index - 1 : index;
  return to === from ? null : {parent, index: to, kind, rect};
}

function rectOf(instance: SpecInstance, node: SpecNode): DropRect | null {
  const el = new SpecNodeRef(instance, node).element();
  if (!el)
    return null;
  const rect = el.getBoundingClientRect();
  return {x: rect.left, y: rect.top, width: rect.width, height: rect.height};
}
