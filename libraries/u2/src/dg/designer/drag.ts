/* The mouse-drag layer over the glass pane — a press anchors, four pixels start a drag, Escape
   cancels, and the drop outlines show where `resolveDrop` would land it — plus where every drop
   onto the tray lands: the palette drag, the `+` wizard and the platform's own channel. HTML5 dnd
   was rejected: it fights `elementsFromPoint` hit-testing, is untestable headless, and its ghost
   and thresholds are uncontrollable. */
import * as grok from 'datagrok-api/grok';
import {div} from '../../core/elements.js';
import type {ReadonlySignal} from '../../core/signals.js';
import {nameForTag, seedNode} from '../../spec/editor.js';
import type {SpecEditor, SpecPatch} from '../../spec/editor.js';
import type {SpecInstance, SpecNode} from '../../spec/spec.js';
import {place, resolveDrop} from './dnd.js';
import type {DropTarget} from './dnd.js';
import {dropNode} from './drop.js';
import type {DropReading} from './drop.js';
import {nodeLabel} from './node-ref.js';
import type {Tray} from './tray.js';

/** A drag in flight: an insert payload from the palette (`tag`) or a move payload from the canvas
 * (`node`), plus what the pointer currently resolves to. */
export interface Drag {
  tag?: string;
  node?: SpecNode;
  label: string;
  x: number;
  y: number;
  active: boolean;
  target: DropTarget | null;
  /** A non-visual tag: the tray is its only target, wherever over the designer it is dropped. */
  tray?: boolean;
  /** Whether the pointer is over the designer at all — read at mouseup, which carries no point. */
  over?: boolean;
}

export interface DragHost {
  readonly root: HTMLElement;
  readonly pane: HTMLElement;
  readonly tray: Tray;
  readonly editor: ReadonlySignal<SpecEditor | undefined>;
  instance(): SpecInstance | undefined;
  hit(e: MouseEvent): SpecNode | null;
  refocus(): void;
  pending(node: SpecNode | null): void;
  apply(patch: SpecPatch): boolean;
  refresh(node: SpecNode): void;
}

/** Pixels the pointer must travel before a press on the pane becomes a drag instead of a click. */
const DRAG_THRESHOLD = 4;
const GHOST_OFFSET = 12;

export class DragLayer {
  readonly into = div([], 'u2-designer-drop-into');
  readonly line = div([], 'u2-designer-drop-line');
  private readonly _ghost = div([], 'u2-designer-drag-ghost');
  drag: Drag | undefined;
  /** A platform drag (a file or a query out of Browse) is over the designer — the other half of
   * {@link dragging}, since that gesture never presses a mouse button inside our DOM. */
  platformDrag = false;
  dragged = false;

  constructor(private readonly _host: DragHost) {}

  /** A drag of either kind is in flight — ours or the platform's. */
  get dragging(): boolean {
    return this.drag !== undefined || this.platformDrag;
  }

  onPaneDown(e: MouseEvent): void {
    if (e.button !== 0)
      return;
    this.dragged = false;
    this._host.refocus();
    const node = this._host.hit(e);
    const instance = this._host.instance();
    if (node && instance && node !== instance.spec.root)
      this.drag = {node, label: nodeLabel(node), x: e.clientX, y: e.clientY, active: false, target: null};
  }

  onPaletteDown(e: MouseEvent): void {
    if (e.button !== 0)
      return;
    const item = (e.target as Element).closest('.u2-palette-item') as HTMLElement | null;
    if (!item)
      return;
    // without this the press starts a text selection that fights the drag
    e.preventDefault();
    const tag = item.dataset.u2Tag;
    this.drag = {tag, label: item.textContent ?? '', x: e.clientX, y: e.clientY,
      active: false, target: null,
      tray: tag !== undefined && this._host.instance()?.registry.get(tag)?.visual === false};
  }

  onMove(e: MouseEvent): void {
    const drag = this.drag;
    const instance = this._host.instance();
    if (!drag || !instance)
      return;
    const root = this._host.root;
    if (!drag.active) {
      if (Math.abs(e.clientX - drag.x) < DRAG_THRESHOLD && Math.abs(e.clientY - drag.y) < DRAG_THRESHOLD)
        return;
      drag.active = true;
      this._ghost.textContent = drag.label;
      root.append(this._ghost);
    }
    const base = root.getBoundingClientRect();
    this._ghost.style.left = `${e.clientX - base.left + root.scrollLeft + GHOST_OFFSET}px`;
    this._ghost.style.top = `${e.clientY - base.top + root.scrollTop + GHOST_OFFSET}px`;
    // the Delphi gesture: a data source dropped anywhere over the designer lands on the tray,
    // which highlights itself for the whole drag
    if (drag.tray) {
      const box = root.getBoundingClientRect();
      drag.over = e.clientX >= box.left && e.clientX <= box.right &&
        e.clientY >= box.top && e.clientY <= box.bottom;
      this._host.tray.setDropTarget(drag.over);
      return;
    }
    // a pointer on the canvas but over no rendered node still means the root (resolveDrop's
    // null-hit fallback) — "drop somewhere in the empty middle" must not silently do nothing
    const hit = this._host.hit(e);
    drag.target = hit !== null || this._overPane(e) ?
      resolveDrop(instance, hit, e.clientX, e.clientY, drag.node) : null;
    this.showDrop(drag.target);
  }

  end(commit: boolean): void {
    const drag = this.drag;
    this.drag = undefined;
    if (!drag || !drag.active)
      return;
    this._ghost.remove();
    this.showDrop(null);
    this._host.tray.setDropTarget(false);
    // the click that follows a real drag is not a selection gesture
    this.dragged = true;
    // a palette drag starts with focus in the filter box; finished, Ctrl+Z/Delete are the canvas's
    this._host.refocus();
    if (!commit)
      return;
    if (drag.tray)
      this._dropComponent(drag);
    else if (drag.target)
      this.drop(drag, drag.target);
  }

  /** Every tray drop appends: chip order is cosmetic, and there is no move-component. */
  private _dropComponent(drag: Drag): void {
    if (drag.over !== true)
      return;
    const tag = drag.tag!;
    this.addComponent(nameForTag(tag), (name, unique) =>
      seedNode(this._host.instance()?.registry.get(tag), tag, name, unique));
  }

  /** What both ways onto the tray end in — the palette drag and the `+` wizard: one seeded,
   * uniquely named component, one patch, one undo step. */
  addComponent(base: string,
    seed: (name: string, unique: (base: string) => string) => SpecNode): void {
    const editor = this._host.editor.peek();
    const instance = this._host.instance();
    if (!editor || !instance)
      return;
    const unique = editor.uniqueNames();
    const node = seed(unique(base), unique);
    this._host.pending(node);
    this._host.apply({op: 'add-component', index: (instance.spec.components ?? []).length, node});
  }

  /** Every dropped item as one seeded source, one compound patch, one undo step. A file is run once
   * on arrival; a function or query follows its own kind policy, untouched. */
  dropSources(reading: DropReading): void {
    if (reading.refused.length > 0)
      grok.shell.warning(reading.refused.join('; '));
    const editor = this._host.editor.peek();
    const instance = this._host.instance();
    if (!editor || !instance || reading.items.length === 0)
      return;
    const unique = editor.uniqueNames();
    const nodes = reading.items.map((item) => dropNode(item, instance.registry, unique));
    const base = (instance.spec.components ?? []).length;
    const patches = nodes.map((node, i): SpecPatch => ({op: 'add-component', index: base + i, node}));
    // every member appends, so only the first can be judged against the tray on screen — the rest
    // carry indices the members before them create. `applyAll` validates each one in the state it
    // actually lands in, and they are the same tag as the one checked here
    const refusal = editor.canApply(patches[0]);
    if (refusal) {
      grok.shell.warning(refusal);
      return;
    }
    this._host.pending(nodes[0]);
    editor.applyAll(patches);
    for (let i = 0; i < nodes.length; i++) {
      if (reading.items[i].runs)
        this._host.refresh(nodes[i]);
    }
  }

  drop(drag: Drag, target: DropTarget): void {
    const editor = this._host.editor.peek();
    if (!editor)
      return;
    if (drag.node) {
      this._host.pending(drag.node);
      this._host.apply({op: 'move', node: drag.node, parent: target.parent, index: target.index});
      return;
    }
    const tag = drag.tag!;
    const unique = editor.uniqueNames();
    const node = seedNode(this._host.instance()?.registry.get(tag), tag, unique(nameForTag(tag)), unique);
    this._host.pending(node);
    this._host.apply({op: 'add', parent: target.parent, index: target.index, node});
  }

  showDrop(target: DropTarget | null): void {
    place(this._host.root, this.into, target && target.kind === 'into' ? target.rect : null);
    place(this._host.root, this.line, target && target.kind === 'line' ? target.rect : null);
  }

  /** Whether the pointer is over the canvas at all: the glass pane is the topmost element there
   * (the adorners are pointer-events: none), so its presence under the point is the answer. */
  private _overPane(e: MouseEvent): boolean {
    return this._host.pane.isConnected &&
      document.elementsFromPoint(e.clientX, e.clientY).includes(this._host.pane);
  }
}
