/* The selection model: the lead, the multi-selection and the hovered node, the adorners that draw
   them over the canvas, and the one `grok.shell.o` assignment every gesture ends in — plus the
   canvas hit-test (`hit`) every gesture starts from, and the tab/accordion headers a click on the
   glass pane is forwarded to. */
import * as grok from 'datagrok-api/grok';
import {div} from '../../core/elements.js';
import type {ReadonlySignal} from '../../core/signals.js';
import type {Action} from '../../components/actions/actions.js';
import {Accordion} from '../../components/containers/accordion.js';
import {TabStrip} from '../../components/containers/tabs.js';
import type {VirtualTree} from '../../components/collections/tree.js';
import type {SpecEditor} from '../../spec/editor.js';
import type {SpecInstance, SpecNode} from '../../spec/spec.js';
import {place} from './dnd.js';
import {refreshPanel} from './handler.js';
import {SpecNodeRef, SpecNodesRef, idPath} from './node-ref.js';
import type {SpecTree} from './node-ref.js';
import type {Tray} from './tray.js';

export interface SelectionHost {
  readonly root: HTMLElement;
  readonly surface: HTMLElement;
  readonly pane: HTMLElement;
  readonly tree: VirtualTree<SpecNode>;
  readonly tray: Tray;
  readonly editor: ReadonlySignal<SpecEditor | undefined>;
  instance(): SpecInstance | undefined;
  model(): SpecTree | undefined;
  dragged(): boolean;
  actions(): Action[];
  report(): void;
  touch(): void;
  own(dispose: () => void): void;
}

export class Selection {
  readonly box = div([], 'u2-designer-selected');
  /** One adorner per selected node; grows on demand, the spares stay hidden. */
  private readonly _adorners = [this.box];
  readonly hover = div([], 'u2-designer-hover');
  selected: SpecNode | null = null;
  /** The whole selection, lead included, in selection order — more than one only after a
   * Ctrl/Shift gesture; every single-select path collapses it back to `[selected]`. */
  multi: SpecNode[] = [];
  hovered: SpecNode | null = null;
  /** Whether the panel was last told the lead is broken — what a re-render may flip. */
  private _shownBroken = false;
  /** The lead's element follows its own size: a splitter panel dragged, a viewer that grew. */
  private readonly _resize = new ResizeObserver(() => this.reposition());
  private _observed: Element | null = null;

  constructor(private readonly _host: SelectionHost) {
    _host.own(() => this._resize.disconnect());
  }

  /** The topmost element under the pointer that the spec owns — the pane itself sits above the
   * canvas, so it is skipped, and the adorners never take pointer events at all. */
  hit(e: MouseEvent): SpecNode | null {
    const instance = this._host.instance();
    if (!instance)
      return null;
    for (const el of document.elementsFromPoint(e.clientX, e.clientY)) {
      if (!this._host.surface.contains(el))
        continue;
      const node = instance.nodeAt(el);
      if (node)
        return node;
    }
    return null;
  }

  /** A click that lands on no node — the designer's own padding — leaves the selection alone, so the
   * context panel and the status bar never disagree about what is current. */
  onCanvasClick(e: MouseEvent): void {
    if (this._host.dragged())
      return;
    if (e.ctrlKey || e.metaKey) {
      const node = this.hit(e);
      if (node)
        this.toggleMulti(node);
      return;
    }
    if (this._forwardHeader(e))
      return;
    const node = this.hit(e);
    if (node)
      this.select(node);
  }

  /** Canvas Ctrl/Cmd+click (F5): toggles the node's membership, the lead following the toggle —
   * the accepted asymmetry is that the tree keeps highlighting the lead's row only. */
  toggleMulti(node: SpecNode): void {
    const multi = this.multi.filter((member) => member !== node);
    if (multi.length === this.multi.length)
      multi.push(node);
    else if (multi.length === 0)
      return;
    if (multi.length === 1) {
      this.select(multi[0]);
      return;
    }
    this.selectMulti(multi, multi.includes(node) ? node : multi[multi.length - 1]);
  }

  /** Tabs and accordion headers stay live at design time: the glass pane owns the click, so a hit
   * on a header activates the pane it names and selects that pane's node — activation plus
   * selection is what a follow-up hit-test would land on anyway. */
  private _forwardHeader(e: MouseEvent): boolean {
    const instance = this._host.instance();
    if (!instance)
      return false;
    const el = document.elementsFromPoint(e.clientX, e.clientY)
      .find((at) => this._host.surface.contains(at));
    if (!el)
      return false;
    const tab = el.closest('.u2-tabs-tab') as HTMLElement | null;
    const header = tab ?? el.closest('.u2-accordion-header') as HTMLElement | null;
    if (!header)
      return false;
    const node = instance.nodeAt(header);
    const built = node ? instance.nodes().get(node) : undefined;
    if (tab && built instanceof TabStrip && tab.dataset.id !== undefined) {
      built.activeTab.value = tab.dataset.id;
      // the pane index off the registration's own id scheme, `tab-${i}`
      const child = node!.children?.[Number(tab.dataset.id.slice(4))];
      this.select(child ?? node!);
      return true;
    }
    if (!tab && built instanceof Accordion) {
      const index = Array.prototype.indexOf.call(built.root.children, header.parentElement);
      const pane = built.paneAt(index);
      if (!pane)
        return false;
      pane.expanded.value = !pane.expanded.peek();
      if (pane.expanded.peek())
        this.select(node!.children?.[index] ?? node!);
      return true;
    }
    return false;
  }

  /** Re-selecting what is already selected is not a no-op: the context panel can miss an update
   * when clicks come faster than the platform renders it, and clicking the row again is how a user
   * asks for it back — so `shell.o` is always re-assigned. */
  select(node: SpecNode): void {
    const instance = this._host.instance();
    if (!instance)
      return;
    this.selected = node;
    this.multi = [node];
    this._reveal(node);
    this.reposition();
    this._host.tray.update(instance, node);
    this._host.report();
    this._host.touch();
    grok.shell.windows.showContextPanel = true;
    this._issue(instance, node);
    const id = this._host.model()?.ids.get(node);
    if (id !== undefined)
      void this._host.tree.expandPath(idPath(id));
  }

  private _issue(instance: SpecInstance, node: SpecNode): void {
    const ref = new SpecNodeRef(instance, node, this._host.editor.peek(), () => this._host.actions());
    this._shownBroken = ref.error() !== null;
    grok.shell.o = ref;
  }

  private _issueMulti(instance: SpecInstance, nodes: SpecNode[]): void {
    grok.shell.o = new SpecNodesRef(nodes.map((node) =>
      new SpecNodeRef(instance, node, this._host.editor.peek())));
  }

  /** The panel again for what is already selected — one assignment: after the Design/Run toggle
   * rebuilt what it renders from (a source, a viewer; or the platform's own gear put a viewer there
   * in Run mode), and after a re-render turned the lead from broken to built or back. The panel on
   * screen is brought up to date in place first: the platform renders the write one step behind. */
  reassert(): void {
    const instance = this._host.instance();
    if (!instance || !this.selected)
      return;
    if (this.multi.length > 1) {
      this._issueMulti(instance, this.multi);
      return;
    }
    refreshPanel(new SpecNodeRef(instance, this.selected));
    this._issue(instance, this.selected);
  }

  /** A panel built for a placeholder shows no status and a stale error once the node builds; one
   * built for a component still says the node is fine once it breaks. */
  followRerender(): void {
    const instance = this._host.instance();
    if (!instance || !this.selected || this.multi.length > 1)
      return;
    if ((new SpecNodeRef(instance, this.selected).error() !== null) !== this._shownBroken)
      this.reassert();
  }

  /** The tree reporting a selection it made itself — including the row `select` just expanded to,
   * which would otherwise assign `grok.shell.o` a second time and make the platform rebuild the
   * property panel twice per patch. A user re-click bypasses this: it goes straight to `select`. */
  onTreeSelection(node: SpecNode | null): void {
    if (!node)
      return;
    if (node !== this.selected)
      this.select(node);
    else if (this.multi.length > 1) {
      // a tree rebuild (setItems) collapsed its set under a live multi-selection — the designer
      // must follow, or Delete would act on members no longer shown anywhere. Deferred past the
      // event dispatch: when a plain re-click on the lead is what collapsed the set, its own
      // listener re-asserts, and reconciling here too would assign `shell.o` twice per gesture.
      queueMicrotask(() => {
        if (this.multi.length > 1 && node === this.selected &&
          this._host.tree.selectedNodes.peek().length <= 1)
          this.select(node);
      });
    }
  }

  /** A multi-selection gesture: the state, the adorner pool, the status and ONE `shell.o`
   * assignment — a `SpecNodesRef` over the members. Only the lead is revealed. */
  selectMulti(nodes: SpecNode[], lead: SpecNode | null): void {
    const instance = this._host.instance();
    if (!instance || nodes.length < 2)
      return;
    this.selected = lead && nodes.includes(lead) ? lead : nodes[nodes.length - 1];
    this.multi = nodes;
    this._reveal(this.selected);
    this.reposition();
    this._host.tray.update(instance, this.selected);
    this._host.report();
    this._host.touch();
    grok.shell.windows.showContextPanel = true;
    this._issueMulti(instance, nodes);
  }

  /** Esc over a multi-selection: back to the lead alone, the tree's own set included. */
  collapse(): void {
    const lead = this.selected!;
    this.multi = [lead];
    this._host.tree.clearSelection();
    this.select(lead);
  }

  /** A selection inside a hidden tab pane or a collapsed accordion pane would adorn a zero rect:
   * every such ancestor opens onto the path before the adorner is placed. */
  private _reveal(node: SpecNode): void {
    const instance = this._host.instance()!;
    for (let at = node, parent = instance.parentOf(at); parent; at = parent, parent = instance.parentOf(at)) {
      const built = instance.nodes().get(parent);
      const index = (parent.children ?? []).indexOf(at);
      if (built instanceof TabStrip)
        built.activeTab.value = `tab-${index}`;
      else if (built instanceof Accordion) {
        const pane = built.paneAt(index);
        if (pane)
          pane.expanded.value = true;
      }
    }
  }

  setHovered(node: SpecNode | null): void {
    if (node === this.hovered)
      return;
    this.hovered = node;
    this.reposition();
  }

  reposition(): void {
    const design = this._host.pane.isConnected;
    const nodes = design && this.selected ? this.multi : [];
    while (this._adorners.length < nodes.length) {
      const box = div([], 'u2-designer-selected');
      this._adorners.push(box);
      this._host.root.append(box);
    }
    for (let i = 0; i < this._adorners.length; i++)
      this._outline(this._adorners[i], nodes[i] ?? null);
    this._outline(this.hover,
      design && this.hovered && !nodes.includes(this.hovered) ? this.hovered : null);
    this._observe(nodes.length > 0 ? this._elementOf(nodes[0]) : null);
  }

  private _observe(el: Element | null): void {
    if (el === this._observed)
      return;
    if (this._observed)
      this._resize.unobserve(this._observed);
    this._observed = el;
    if (el)
      this._resize.observe(el);
  }

  private _elementOf(node: SpecNode | null): HTMLElement | null {
    const instance = this._host.instance();
    return (node && instance ? new SpecNodeRef(instance, node).element() : undefined) ?? null;
  }

  private _outline(box: HTMLElement, node: SpecNode | null): void {
    const el = this._elementOf(node);
    place(this._host.root, box, el && el.isConnected ? el.getBoundingClientRect() : null);
  }
}
