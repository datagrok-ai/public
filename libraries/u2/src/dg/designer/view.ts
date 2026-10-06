/* The designer as a citizen of the shell (DD1): the canvas is a live `SpecInstance`, the toolbox is
   the palette over its registry plus its structure tree, the ribbon holds the actions and the
   Design/Run toggle, and the status bar carries the selection path. Design mode is one glass pane
   over the canvas — the components render exactly as they run, and the pane hit-tests through
   `SpecInstance.nodeAt` instead of fighting each control's own handlers. Every edit is a patch
   through `SpecEditor`: this file is the gesture layer, the engine owns validation and undo. */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Control} from '../../core/component.js';
import {signal, ReadonlySignal} from '../../core/signals.js';
import {div, divV, h3} from '../../core/elements.js';
import type {Action} from '../../components/actions/actions.js';
import {actionsMenu} from '../../components/actions/actions.js';
import {ButtonGroup} from '../../components/actions/button-group.js';
import {TreeNode, VirtualTree} from '../../components/collections/tree.js';
import {SpecContext, SpecInstance, parseSpec, renderSpec} from '../../spec/spec.js';
import type {Spec, SpecNode} from '../../spec/spec.js';
import {registry as globalRegistry} from '../../spec/registry.js';
import type {Registry} from '../../spec/registry.js';
import {SpecEditor} from '../../spec/editor.js';
import type {AppliedPatch, SpecPatch} from '../../spec/editor.js';
import {appView} from '../shell/app-view.js';
import {disposePanel} from './handler.js';
import {accepts} from './dnd.js';
import {makeDesignerDroppable} from './drop.js';
import {DragLayer} from './drag.js';
import type {Drag} from './drag.js';
import {Selection} from './selection.js';
import {DesignerActions} from './actions.js';
import {dragKey, keyDown} from './keys.js';
import {Ribbon} from './ribbon.js';
import {Palette} from './palette.js';
import {Tray} from './tray.js';
import {SpecNodeRef, SpecTree, brokenCount, specTree} from './node-ref.js';

export interface DesignerViewOptions {
  name?: string;
  /** What `$.` paths and `cmd:` names resolve against while the spec runs. */
  ctx?: SpecContext;
  /** The components the canvas may build and the palette offers; the global registry by default. */
  registry?: Registry;
}

/** The canvas: a rendered spec plus the design-mode chrome over it. */
export class SpecDesigner extends Control {
  readonly status: ReadonlySignal<string>;
  readonly toolbox: HTMLElement;

  private readonly _ctx: SpecContext;
  private readonly _registry: Registry;
  private readonly _surface = div([], 'u2-designer-surface');
  private readonly _pane = div([], 'u2-designer-pane');
  private readonly _tree: VirtualTree<SpecNode>;
  private readonly _tray: Tray;
  private readonly _palette: Palette;
  private readonly _mode: ButtonGroup;
  private readonly _selection: Selection;
  private readonly _dragLayer: DragLayer;
  private readonly _verbs: DesignerActions;
  private readonly _ribbon: Ribbon;
  private readonly _status = signal('');
  private readonly _outlines = signal(true);
  private readonly _editor = signal<SpecEditor | undefined>(undefined);
  /** Bumped by anything that changes what the actions would apply to — the ribbon reads it. */
  private readonly _revision = signal(0);
  private _instance: SpecInstance | undefined;
  private _model: SpecTree | undefined;
  /** What the gesture in flight wants selected once its patch has landed — the node it inserts, or
   * the parent of the node it removes. The rebuild reads it, so one patch is one selection. */
  private _pendingSelect: SpecNode | null = null;

  constructor(spec: string | object, ctx?: SpecContext, registry?: Registry) {
    super();
    this.status = this._status;
    this._ctx = ctx ?? new SpecContext();
    this._registry = registry ?? globalRegistry;
    this.root.classList.add('u2-designer');
    this.root.dataset.u2 = 'designer';
    // the shortcuts below are the view's own: the root has to be able to hold focus for them
    this.root.tabIndex = 0;
    this._tray = this.runInScope(() => new Tray({
      onSelect: (node) => this._select(node),
      onContext: (node, x, y) => {
        this._select(node);
        actionsMenu(this._verbs.handBack(this._actions())).show({x, y});
      },
      onAdd: (base, seed) => this._dragLayer.addComponent(base, seed),
    }));
    this._tree = this.runInScope(() => new VirtualTree<SpecNode>({
      onRename: (node, label) => this._rename(node, label),
      contextActions: (node) => this._rowActions(node),
    }));
    this._mode = this.runInScope(() => new ButtonGroup({
      items: [{id: 'design', label: 'Design'}, {id: 'run', label: 'Run'}],
      toggle: 'single',
      density: 'ribbon',
    }));
    this._mode.selected.value = ['design'];

    const host = {
      root: this.root,
      surface: this._surface,
      pane: this._pane,
      tree: this._tree,
      tray: this._tray,
      mode: this._mode.root,
      outlines: this._outlines,
      editor: this._editor,
      instance: () => this._instance,
      model: () => this._model,
      selected: () => this._selected,
      multi: () => this._multi,
      revision: () => this._revision.value,
      inDrag: () => this._drag !== undefined,
      dragged: () => this._dragLayer.dragged,
      hit: (e: MouseEvent) => this._selection.hit(e),
      select: (node: SpecNode) => this._select(node),
      collapse: () => this._selection.collapse(),
      endDrag: (commit: boolean) => this._endDrag(commit),
      actions: () => this._actions(),
      run: (name: string) => this._run(name),
      pending: (node: SpecNode | null) => this._pendingSelect = node,
      apply: (patch: SpecPatch) => this._apply(patch),
      refresh: (node: SpecNode) => this._verbs.refresh(node),
      open: (spec: string | object) => this.open(spec),
      dump: () => this.dump(),
      effect: (fn: () => void) => this.effect(fn),
      own: (dispose: () => void) => this.own(dispose),
      report: () => this._report(),
      touch: () => this._touch(),
      refocus: () => this._refocus(),
    };
    this._selection = new Selection(host);
    this._dragLayer = new DragLayer(host);
    this._verbs = new DesignerActions(host);
    this._ribbon = new Ribbon(host);

    // in flow, after the surface: the tray sticks to the bottom edge and hides itself in Run mode
    this.root.append(this._surface, this._tray.root, this._selection.box, this._selection.hover,
      this._dragLayer.into, this._dragLayer.line);
    // the platform's drag channel, on the root rather than the pane: the pane is removed in Run
    // mode and the dart-side registration has no counterpart to undo, so it must be wired on the
    // element that lives as long as the designer — and every callback asks `active()` first
    makeDesignerDroppable({
      element: this.root,
      active: () => !this.scope.isDisposed && this._pane.isConnected && this._drag === undefined,
      onDragActive: (on) => {
        this._dragLayer.platformDrag = on;
        this._tray.setDropTarget(on);
        if (on)
          this._selection.setHovered(null);
      },
      onDrop: (reading) => this._dragLayer.dropSources(reading),
    });

    // a click on the row that is already selected moves no signal, so the selection effect below
    // never re-fires — the click itself has to re-assert it, past the guard that effect needs;
    // a modifier click is a multi-selection gesture, and re-asserting would collapse the set (F5)
    this._listen(this._tree.root, 'click', (e) => {
      const me = e as MouseEvent;
      if (me.ctrlKey || me.metaKey || me.shiftKey)
        return;
      const node = this._tree.selectedNode.peek()?.data;
      if (node)
        this._select(node);
    });
    this._listen(this._pane, 'click', (e) => this._selection.onCanvasClick(e as MouseEvent));
    this._listen(this._pane, 'mousedown', (e) => this._dragLayer.onPaneDown(e as MouseEvent));
    this._listen(this._pane, 'contextmenu',
      (e) => this._verbs.onContextMenu(e as MouseEvent, this._selection.hit(e as MouseEvent)));
    this._listen(this._pane, 'mousemove',
      (e) => this._selection.setHovered(this._dragging ? null : this._selection.hit(e as MouseEvent)));
    this._listen(this._pane, 'mouseleave', () => this._selection.setHovered(null));
    // the second hover source: a tree row highlights its node on the canvas (run mode draws
    // nothing — `reposition` suppresses the adorners while the pane is detached)
    this._listen(this._tree.root, 'mouseover', (e) =>
      this._selection.setHovered(this._dragging ? null : this._tree.nodeForRow(e.target)?.data ?? null));
    this._listen(this._tree.root, 'mouseleave', () => this._selection.setHovered(null));
    // the platform wraps every JS view in a widget whose ancestor-level bubble contextmenu
    // listener opens its own 'Properties...' menu — nothing may propagate past the designer
    // while it is designing; Run mode keeps the platform menu
    // preventDefault too: a KEYBOARD contextmenu (Shift+F10) targets the focused designer root
    // itself, below the pane's own handler, and would otherwise pop the browser menu
    this._listen(this.root, 'contextmenu', (e) => {
      if (this._pane.isConnected) {
        e.preventDefault();
        e.stopPropagation();
      }
    });
    this._listen(this._surface, 'scroll', () => this._selection.reposition(), true);
    this._listen(window, 'resize', () => this._selection.reposition());
    // a drag that starts on a palette row leaves the pane at once, and Escape has to reach it
    // wherever the focus is; all three do nothing while no drag is in flight
    this._listen(document, 'mousemove', (e) => this._dragLayer.onMove(e as MouseEvent));
    this._listen(document, 'mouseup', () => this._endDrag(true));
    this._listen(document, 'keydown', (e) => dragKey(e as KeyboardEvent, host));

    this.open(spec);

    this._palette = this.runInScope(() => new Palette(this._instance?.registry ?? this._registry));
    this.toolbox = divV([h3('Palette'), this._palette.root, h3('Structure'), this._tree.root],
      'u2-designer-toolbox');
    this._listen(this._palette.root, 'mousedown', (e) => this._dragLayer.onPaletteDown(e as MouseEvent));
    this._listen(this.toolbox, 'contextmenu', (e) => {
      // keyboard contextmenu with focus on the tree/palette chrome — no menu of any kind
      e.preventDefault();
      e.stopPropagation();
    });
    this._listen(this._mode.root, 'click', () => this._refocus());
    this._listen(this.root, 'keydown', (e) => keyDown(e as KeyboardEvent, host));
    this._listen(this.toolbox, 'keydown', (e) => keyDown(e as KeyboardEvent, host));

    this.effect(() => {
      const rows = this._tree.selectedNodes.value;
      const lead = this._tree.selectedNode.value?.data ?? null;
      if (rows.length > 1)
        this._selection.selectMulti(rows.map((row) => row.data!), lead);
      else
        this._selection.onTreeSelection(lead);
    });
    this.effect(() => this._setDesign(this._mode.selected.value[0] !== 'run'));
    this.effect(() => this.root.classList.toggle('u2-designer-outlines', this._outlines.value));
    this.effect(() => {
      const applied = this._editor.value?.onDidApply.value;
      if (applied)
        this._afterPatch(applied);
    });
    // the context panel the handler rendered for this designer's selections dies with it —
    // otherwise its sync effect keeps refreshing a detached form for the rest of the session
    this.own(disposePanel);
  }

  ribbon(): HTMLElement[] {
    return this.runInScope(() => this._ribbon.build());
  }

  /** Renders `spec` in place of what the canvas holds, with a fresh editor — an instance and its
   * edit history die together. A spec that is not a spec is refused and the current one stays up.
   * What the canvas holds is always this designer's own copy: the editor patches the document in
   * place (DD10), so a caller's object stays the constant it was — `dump()` is how the edits leave. */
  open(spec: string | object): void {
    let parsed: Spec;
    try {
      parsed = parseSpec(spec);
      // a string parsed into an object nobody else holds; only an object input needs copying
      if (typeof spec !== 'string')
        parsed = JSON.parse(JSON.stringify(parsed)) as Spec;
    } catch (e) {
      grok.shell.error(e instanceof Error ? e.message : String(e));
      return;
    }
    this._instance?.dispose();
    this._surface.textContent = '';
    this._selection.selected = null;
    this._selection.multi = [];
    this._hovered = null;
    this._pendingSelect = null;
    this._instance = this.runInScope(() => renderSpec(parsed, this._ctx, this._registry,
      {designTime: this._mode.selected.peek()[0] !== 'run'}));
    this._surface.append(this._instance.root);
    this._editor.value = new SpecEditor(this._instance);
    this._model = specTree(this._instance);
    // cleared BEFORE setRoots: the list preserves selection by row key on items change, and
    // tree ids are ordinals that collide across specs — the old selection's id would silently
    // re-select an arbitrary node of the new spec (and announce it) before root wins
    this._tree.clearSelection();
    this._tree.setRoots(this._model.roots);
    this._stampDesign();
    this._tray.update(this._instance, null);
    this._report();
    // the root is selected outright: the panel shows something from the first frame, and the tree
    // never carries a row selected against the spec that has just been replaced
    const root = this._model.roots[0]?.data;
    if (root)
      this._select(root);
    else
      this._selection.reposition();
  }

  dump(): object | null {
    return this._instance ? this._instance.dump() : null;
  }

  private get _drag(): Drag | undefined { return this._dragLayer.drag; }
  private set _drag(drag: Drag | undefined) { this._dragLayer.drag = drag; }
  private get _dragging(): boolean { return this._dragLayer.dragging; }
  private get _selected(): SpecNode | null { return this._selection.selected; }
  private get _multi(): SpecNode[] { return this._selection.multi; }
  private get _hovered(): SpecNode | null { return this._selection.hovered; }
  private set _hovered(node: SpecNode | null) { this._selection.hovered = node; }
  private _select(node: SpecNode): void { this._selection.select(node); }
  private _endDrag(commit: boolean): void { this._dragLayer.end(commit); }
  private _actions(): Action[] { return this._verbs.actions(); }
  private _run(name: string): void { this._verbs.run(name); }
  private _rowActions(row: TreeNode<SpecNode>): Action[] { return this._verbs.rowActions(row); }

  /** Every applied patch lands here, from a gesture, a keystroke or undo alike: the adorners and
   * the status line always follow, the tree and the selection only when the structure moved. */
  private _afterPatch(applied: AppliedPatch): void {
    if (applied.structural)
      this._rebuildTree();
    else
      this._selection.followRerender();
    this._selection.reposition();
    this._tray.update(this._instance, this._selected);
    this._report();
    this._stampDesign();
    this._touch();
  }

  private _rebuildTree(): void {
    const instance = this._instance;
    if (!instance)
      return;
    this._model = specTree(instance);
    this._tree.clearSelection();
    this._tree.setRoots(this._model.roots);
    const survivor = this._survivor(instance);
    this._pendingSelect = null;
    if (survivor)
      this._select(survivor);
  }

  /** What stays selected across a structural patch: what the gesture asked for, the same node where
   * it is still rendered otherwise, and the root as the last resort. */
  private _survivor(instance: SpecInstance): SpecNode | null {
    const nodes = instance.nodes();
    if (this._pendingSelect && nodes.has(this._pendingSelect))
      return this._pendingSelect;
    if (this._selected && nodes.has(this._selected))
      return this._selected;
    return nodes.has(instance.spec.root) ? instance.spec.root : null;
  }

  /** The design-time classes on what the spec rendered: the root's always-on bounds, the
   * per-control outline (toggleable via `u2-designer-outlines` on the designer root) and the
   * empty-container affordance that makes an invisible drop target visible. Run mode keeps the
   * classes and drops the styling (css/designer.css). */
  private _stampDesign(): void {
    const instance = this._instance;
    if (!instance)
      return;
    for (const [node, built] of instance.nodes()) {
      const el = SpecInstance.elementOf(built);
      // a tray component has no element to stamp
      if (!el)
        continue;
      el.classList.add('u2-designer-node');
      el.classList.toggle('u2-designer-root', node === instance.spec.root);
      el.classList.toggle('u2-designer-empty',
        (node.children?.length ?? 0) === 0 && accepts(instance.registry, node));
    }
  }

  /** Never `_revision.value++`: an effect that reads the signal to rebuild the ribbon would count
   * the read-modify-write as a dependency on itself. */
  private _touch(): void {
    this._revision.value = this._revision.peek() + 1;
  }

  /** A ribbon or menu click parks focus on chrome outside the designer, where the next Ctrl+Z
   * belongs to the platform ("Nothing to undo") — a finished action hands focus back. */
  private _refocus(): void {
    this.root.focus({preventScroll: true});
  }

  /** Every gesture goes through here: the engine is the one authority on what may apply, and its
   * refusal is what the user is told. */
  private _apply(patch: SpecPatch): boolean {
    const editor = this._editor.peek();
    if (!editor)
      return false;
    const refusal = editor.canApply(patch);
    if (refusal) {
      this._pendingSelect = null;
      grok.shell.warning(refusal);
      return false;
    }
    editor.apply(patch);
    return true;
  }

  private _rename(row: TreeNode<SpecNode>, label: string): void {
    if (!row.data)
      return;
    // the tree writes the new label as it commits — a refused rename is undone by rebuilding it
    if (!this._apply({op: 'rename', node: row.data, name: label}))
      this._rebuildTree();
  }

  private _setDesign(design: boolean): void {
    if (design)
      this.root.append(this._pane);
    else {
      this._pane.remove();
      this._hovered = null;
    }
    this.root.classList.toggle('u2-designer-running', !design);
    // the sources build for the mode they are in (DD9) — and everything bound to one is rendered
    // again with it, so the toggle really does swap in the live context
    const toggled = this._instance !== undefined && this._instance.designTime !== design;
    this._instance?.setDesignTime(design);
    this._stampDesign();
    this._tray.update(this._instance, this._selected);
    this._selection.reposition();
    // the panel was rendered from what the old mode built — a source whose status it still reads,
    // or whatever the platform's own gear put there in Run mode
    if (toggled)
      this._selection.reassert();
  }

  private _report(): void {
    if (!this._instance) {
      this._status.value = 'Nothing rendered';
      return;
    }
    const count = this._instance.nodes().size;
    const parts = [this._multi.length > 1 ? `${this._multi.length} selected` :
      this._selected ? new SpecNodeRef(this._instance, this._selected).path() : 'No selection',
    `${count} node${count === 1 ? '' : 's'}`];
    const broken = brokenCount(this._instance);
    if (broken > 0)
      parts.push(`${broken} broken`);
    if (this._editor.peek()?.dirty.peek())
      parts.push('modified');
    this._status.value = parts.join(' · ');
  }

  private _listen(target: EventTarget, type: string, handler: (e: Event) => void, capture = false): void {
    target.addEventListener(type, handler, capture);
    this.own(() => target.removeEventListener(type, handler, capture));
  }
}

/** The designer as a platform view: canvas as content, palette and structure tree in the toolbox,
 * actions and the Design/Run toggle in the ribbon, selection path in the status bar. */
export function designerView(spec: string | object, options: DesignerViewOptions = {}): DG.ViewBase {
  const designer = new SpecDesigner(spec, options.ctx, options.registry);
  return appView({
    name: options.name ?? 'Designer',
    content: designer,
    toolbox: designer.toolbox,
    ribbon: [designer.ribbon()],
    status: designer.status,
  });
}
