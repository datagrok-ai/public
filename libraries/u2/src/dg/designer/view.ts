/* The designer as a citizen of the shell (DD1): the canvas is a live `SpecInstance`, the toolbox is
   its structure tree, the ribbon holds the actions and the Design/Run toggle, and the status bar
   carries the selection path. Design mode is one glass pane over the canvas — the components render
   exactly as they run, and the pane hit-tests through `SpecInstance.nodeAt` instead of fighting each
   control's own handlers. Read-only at P2: selecting is the only gesture. */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Control} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {signal, ReadonlySignal} from '../../core/signals.js';
import {button, div, divV, h3} from '../../core/elements.js';
import {ButtonGroup} from '../../components/button-group.js';
import {Dialog} from '../../components/dialog.js';
import {TextArea} from '../../components/text-input.js';
import {VirtualTree} from '../../components/tree.js';
import {SpecContext, SpecInstance, parseSpec, renderSpec} from '../../spec/spec.js';
import type {Spec, SpecNode} from '../../spec/spec.js';
import {appView} from '../app-view.js';
import {SpecNodeRef, SpecTree, brokenCount, idPath, specTree} from './node-ref.js';

export interface DesignerViewOptions {
  name?: string;
  /** What `$.` paths and `cmd:` names resolve against while the spec runs. */
  ctx?: SpecContext;
}

/** The canvas: a rendered spec plus the design-mode chrome over it. */
export class SpecDesigner extends Control {
  readonly status: ReadonlySignal<string>;
  readonly toolbox: HTMLElement;

  private readonly _ctx: SpecContext;
  private readonly _surface = div([], 'u2-designer-surface');
  private readonly _pane = div([], 'u2-designer-pane');
  private readonly _selection = div([], 'u2-designer-selected');
  private readonly _hover = div([], 'u2-designer-hover');
  private readonly _tree: VirtualTree<SpecNode>;
  private readonly _mode: ButtonGroup;
  private readonly _status = signal('');
  private _instance: SpecInstance | undefined;
  private _model: SpecTree | undefined;
  private _selected: SpecNode | null = null;
  private _hovered: SpecNode | null = null;

  constructor(spec: string | object, ctx?: SpecContext) {
    super();
    this.status = this._status;
    this._ctx = ctx ?? new SpecContext();
    this.root.classList.add('u2-designer');
    this.root.dataset.u2 = 'designer';
    this.root.append(this._surface, this._selection, this._hover);

    this._tree = this.run(() => new VirtualTree<SpecNode>({}));
    this.toolbox = divV([h3('Structure'), this._tree.root], 'u2-designer-toolbox');
    this._mode = this.run(() => new ButtonGroup({
      items: [{id: 'design', label: 'Design'}, {id: 'run', label: 'Run'}],
      toggle: 'single',
      density: 'ribbon',
    }));
    this._mode.selected.value = ['design'];

    // a click on the row that is already selected moves no signal, so the selection effect below
    // never re-fires — the click itself has to re-assert it
    this._listen(this._tree.root, 'click',
      () => this._onTreeSelection(this._tree.selectedNode.peek()?.data ?? null));
    this._listen(this._pane, 'click', (e) => this._onCanvasClick(e as MouseEvent));
    this._listen(this._pane, 'mousemove', (e) => this._setHovered(this._hit(e as MouseEvent)));
    this._listen(this._pane, 'mouseleave', () => this._setHovered(null));
    this._listen(this._surface, 'scroll', () => this._reposition(), true);
    this._listen(window, 'resize', () => this._reposition());

    this.open(spec);
    this.effect(() => this._onTreeSelection(this._tree.selectedNode.value?.data ?? null));
    this.effect(() => this._setDesign(this._mode.selected.value[0] !== 'run'));
  }

  ribbon(): HTMLElement[] {
    return this.run(() => [
      button('Open…', () => this._openDialog()),
      button('Copy spec', () => this._copy()),
      this._mode.root,
    ]);
  }

  /** Renders `spec` in place of what the canvas holds. A spec that is not a spec is refused and the
   * current one stays up — the designer never leaves an empty canvas behind. */
  open(spec: string | object): void {
    let parsed: Spec;
    try {
      parsed = parseSpec(spec);
    } catch (e) {
      grok.shell.error(e instanceof Error ? e.message : String(e));
      return;
    }
    this._instance?.dispose();
    this._surface.textContent = '';
    this._selected = null;
    this._hovered = null;
    this._instance = this.run(() => renderSpec(parsed, this._ctx));
    this._surface.append(this._instance.root);
    this._model = specTree(this._instance);
    // cleared BEFORE setRoots: the list preserves selection by row key on items change, and
    // tree ids are ordinals that collide across specs — the old selection's id would silently
    // re-select an arbitrary node of the new spec (and announce it) before root wins
    this._tree.clearSelection();
    this._tree.setRoots(this._model.roots);
    this._report();
    // the root is selected outright: the panel shows something from the first frame, and the tree
    // never carries a row selected against the spec that has just been replaced
    const root = this._model.roots[0]?.data;
    if (root)
      this._select(root);
    else
      this._reposition();
  }

  dump(): object | null {
    return this._instance ? this._instance.dump() : null;
  }

  /** A click that lands on no node — the designer's own padding — leaves the selection alone, so the
   * context panel and the status bar never disagree about what is current. */
  private _onCanvasClick(e: MouseEvent): void {
    const node = this._hit(e);
    if (node)
      this._select(node);
  }

  /** Re-selecting what is already selected is not a no-op: the context panel can miss an update
   * when clicks come faster than the platform renders it, and clicking the row again is how a user
   * asks for it back — so `shell.o` is always re-assigned. */
  private _select(node: SpecNode): void {
    if (!this._instance)
      return;
    this._selected = node;
    this._reposition();
    this._report();
    grok.shell.windows.showContextPanel = true;
    grok.shell.o = new SpecNodeRef(this._instance, node);
    const id = this._model?.ids.get(node);
    if (id !== undefined)
      void this._tree.expandPath(idPath(id));
  }

  private _onTreeSelection(node: SpecNode | null): void {
    if (node)
      this._select(node);
  }

  private _setHovered(node: SpecNode | null): void {
    if (node === this._hovered)
      return;
    this._hovered = node;
    this._reposition();
  }

  /** The topmost element under the pointer that the spec owns — the pane itself sits above the
   * canvas, so it is skipped, and the adorners never take pointer events at all. */
  private _hit(e: MouseEvent): SpecNode | null {
    if (!this._instance)
      return null;
    for (const el of document.elementsFromPoint(e.clientX, e.clientY)) {
      if (!this._surface.contains(el))
        continue;
      const node = this._instance.nodeAt(el);
      if (node)
        return node;
    }
    return null;
  }

  private _setDesign(design: boolean): void {
    if (design)
      this.root.append(this._pane);
    else {
      this._pane.remove();
      this._hovered = null;
    }
    this.root.classList.toggle('u2-designer-running', !design);
    this._reposition();
  }

  private _reposition(): void {
    const design = this._pane.isConnected;
    this._outline(this._selection, design ? this._selected : null);
    this._outline(this._hover, design && this._hovered !== this._selected ? this._hovered : null);
  }

  private _outline(box: HTMLElement, node: SpecNode | null): void {
    const el = node && this._instance ? new SpecNodeRef(this._instance, node).element() : undefined;
    if (!el || !el.isConnected) {
      box.style.display = 'none';
      return;
    }
    const rect = el.getBoundingClientRect();
    const base = this.root.getBoundingClientRect();
    box.style.display = 'block';
    box.style.left = `${rect.left - base.left + this.root.scrollLeft}px`;
    box.style.top = `${rect.top - base.top + this.root.scrollTop}px`;
    box.style.width = `${rect.width}px`;
    box.style.height = `${rect.height}px`;
  }

  private _report(): void {
    if (!this._instance) {
      this._status.value = 'Nothing rendered';
      return;
    }
    const count = this._instance.nodes().size;
    const parts = [this._selected ? new SpecNodeRef(this._instance, this._selected).path() : 'No selection',
      `${count} node${count === 1 ? '' : 's'}`];
    const broken = brokenCount(this._instance);
    if (broken > 0)
      parts.push(`${broken} broken`);
    this._status.value = parts.join(' · ');
  }

  private _openDialog(): void {
    const scope = new Scope();
    const spec = this.dump();
    const area = Scope.runWith(scope, () => new TextArea({
      label: 'Spec',
      value: spec ? JSON.stringify(spec, null, 2) : '',
      placeholder: 'Paste a dg-ui/1 spec',
    }));
    area.root.classList.add('u2-designer-spec-editor');
    Scope.runWith(scope, () => new Dialog('Open spec'))
      .add(area)
      .onOK(() => {
        this.open(area.value.peek());
        scope.dispose();
      })
      .onCancel(() => scope.dispose())
      .show({modal: true, width: 560});
  }

  private _copy(): void {
    const spec = this.dump();
    if (!spec) {
      grok.shell.info('Nothing to copy — the canvas is empty.');
      return;
    }
    navigator.clipboard.writeText(JSON.stringify(spec, null, 2))
      .then(() => grok.shell.info('Spec copied to the clipboard.'),
        (e: unknown) => grok.shell.error(`Could not copy: ${e instanceof Error ? e.message : String(e)}`));
  }

  private _listen(target: EventTarget, type: string, handler: (e: Event) => void, capture = false): void {
    target.addEventListener(type, handler, capture);
    this.own(() => target.removeEventListener(type, handler, capture));
  }
}

/** The designer as a platform view: canvas as content, structure tree in the toolbox, actions and
 * the Design/Run toggle in the ribbon, selection path in the status bar. */
export function designerView(spec: string | object, options: DesignerViewOptions = {}): DG.ViewBase {
  const designer = new SpecDesigner(spec, options.ctx);
  return appView({
    name: options.name ?? 'Designer',
    content: designer,
    toolbox: designer.toolbox,
    ribbon: [designer.ribbon()],
    status: designer.status,
  });
}
