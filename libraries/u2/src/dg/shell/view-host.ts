import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Component, Control, component} from '../../core/component.js';
import type {BuiltRoot} from '../../core/component.js';
import {Scope, Usage} from '../../core/scope.js';
import {Ownership} from '../../core/owner.js';
import {signal, computed, Signal, ReadonlySignal} from '../../core/signals.js';
import {dfBindings} from '../../sources/df-bindings.js';
import type {DataFrameLike} from '../../sources/df-bindings.js';
import type {Dialog, DialogShowOptions} from '../../components/containers/dialog.js';
import {TableViewChrome} from './table-view-chrome.js';
import type {RibbonOptions} from './table-view-chrome.js';
import {host, isView} from './host.js';

/** A view's state as signals: {@link ViewHost.signals}. */
export interface ViewSignals {
  /** The shell shows this view now. */
  readonly isCurrent: ReadonlySignal<boolean>;
  /** A TableView's grid: a new one after a layout is applied, which an effect over it re-subscribes for; null
   * for any other view. */
  readonly grid: ReadonlySignal<DG.Grid | null>;
  /** A TableView's table; null for any other view. */
  readonly dataFrame: ReadonlySignal<DG.DataFrame | null>;
  /** The table's current row (-1 for none); a write moves it. */
  readonly currentRowIdx: Signal<number>;
  /** The selection and the filter: the table's live bitsets, firing on every change of them. */
  readonly selection: ReadonlySignal<DG.BitSet | undefined>;
  readonly filter: ReadonlySignal<DG.BitSet | undefined>;
  /** The column names, in order: fires when columns are added, removed or renamed. */
  readonly columns: ReadonlySignal<string[]>;
  readonly rowCount: ReadonlySignal<number>;
}

/** What the shell's chrome around a view holds: a control, disposed with the view, or an element. */
export type ChromeItem = Control | HTMLElement;

/** What {@link ViewHost.create} makes a view of. */
export interface ViewHostOptions {
  name: string;
  /** Builds the view under the host, so what it builds goes with the view: the content, or the content with the
   * shell's chrome around it ({@link ViewContent}). A control it answers is the view's root; an element, a text or
   * several go into a `div`. The view is made from the content: use the host's view, signals, `show` and `dock` once
   * `create` has answered (in a handler, a `host.task`). */
  content: (host: ViewHost) => BuiltRoot | ViewContent;
}

/** The view's content and the shell's chrome around it, as a {@link ViewHostOptions.content} builder answers them:
 * the chrome is built with the content, under the host, and goes with it. */
export interface ViewContent {
  root: BuiltRoot;
  /** Ribbon panel groups → `view.setRibbonPanels`. Main view controls (filter, refresh, mode switches, menu bar)
   * belong here, not inside the content area. */
  ribbon?: ChromeItem[][];
  /** → the view's toolbox pane: navigation and structure, not view controls. */
  toolbox?: ChromeItem;
  /** → the shell's per-view status bar. A signal renders as a live text panel; elements and controls are placed as
   * given. */
  status?: ReadonlySignal<string> | ChromeItem | ChromeItem[];
  /** The view's URL path, written verbatim to `view.path`. For an app, include the app root (`/apps/<Package>/<App>`):
   * the platform prefixes nothing for a JS app. A signal is mirrored on every change, so in-app navigation round-trips
   * to the address bar (the shell replaces rather than pushes, so it costs no history entries). */
  path?: string | ReadonlySignal<string>;
}

/** How {@link ViewHost.dock} docks a panel. */
export interface ViewDockOptions {
  /** The panel's caption. */
  title?: string;
  /** The share of the area it takes (0.5 by default). */
  ratio?: number;
}

/** A platform view as the owner of what code living in it builds: its content, its ribbon group, its signals, its
 * dialogs and its docked panels all go with the view. `ViewHost.of(view)` answers one host per view, whichever
 * wrapper the platform handed out (it makes a new one each time); `ViewHost.create` makes a new view. The owner verbs
 * are the view's: `host.build(fn)`, `host.task(async (t) => …)` (after an `await`, `t.build(…)`),
 * `host.subscribe(stream, handler)`, `host.effect(fn)`, `host.own(cleanup)`. The host is disposed when the shell
 * removes the view (and a created view's, when the platform kills its content); `dispose()` ends it early.
 * @example
 * const host = ViewHost.of(tv);
 * host.ribbon({actions: [{icon: 'sync', tooltip: 'Refresh Campaign', run: () => refresh()}]});
 * host.effect(() => watch(host.signals.grid.value));
 * host.show(new Dialog('Progress tracker'), {modal: true}); */
export class ViewHost extends Component {
  static isRoot = true;
  private static readonly _hosts = new WeakMap<object, ViewHost>();
  private _view: DG.ViewBase | undefined;
  private _signals: ViewSignals | undefined;

  private constructor() {
    super();
  }

  /** The host of `view`, made on first ask: the same one for every wrapper of the view. */
  static of(view: DG.ViewBase): ViewHost {
    return ViewHost._hosts.get(ViewHost._key(view)) ?? new ViewHost()._attach(view);
  }

  /** A new platform view: what `options.content` builds — the content, and the shell's chrome around it (ribbon,
   * toolbox, status bar) with the view's path. Answers its host; `host.view` is the view to answer from an app
   * function or to add to the shell. The content is released through the platform's kill channel when the view
   * closes, and the chrome with it (the shell keeps the chrome outside the view root, where the kill never reaches).
   * @example
   * const host = ViewHost.create({name: 'Reports', content: () => {
   *   const page = reportsPage();
   *   return {root: page, ribbon: [[page.search]], status: page.summary, path: '/apps/Reports'};
   * }});
   * grok.shell.addView(host.view); */
  static create(options: ViewHostOptions): ViewHost {
    const self = new ViewHost();
    let built: BuiltRoot | ViewContent;
    try {
      built = self.build(() => options.content(self));
    }
    catch (e) {
      self.dispose();
      throw e;
    }
    const parts = ViewHost._isContent(built) ? built : {root: built};
    const root = parts.root;
    const content = Control.is(root) ? root : self.build(() => component(() => root));
    self.own(() => content.dispose());
    content.own(() => self.dispose());
    // the app's automation scope root: page-level selectors compose from the view down; a named content wins
    if (content.name === undefined)
      content.name = options.name;
    const view = DG.View.fromRoot(host(content));
    view.name = options.name;
    self._attach(view);
    self._chrome(view, parts);
    return self;
  }

  /** The platform view. */
  get view(): DG.ViewBase {
    return this._view ?? Usage.fail('view-not-ready', 'ViewHost.create() makes the view from its content once the ' +
      'content is built: use the host\'s view, ribbon, signals, show and dock after create() answers — in a handler ' +
      'or a `host.task`');
  }

  /** The view's state as signals, owned by the host and made on first ask: whether it is current, a TableView's grid
   * (an effect over it re-subscribes per grid), table, current row (a write moves it), selection, filter, columns and
   * row count. */
  get signals(): ViewSignals {
    return this._signals ??= viewSignalsOn(this.view, this.scope.child('view-signals'));
  }

  /** The package's ribbon group over the view — icons, buttons, the crumb back, the Save policy and a status line
   * ({@link TableViewChrome}) — disposed with the view, or before it by its own `dispose()`. */
  ribbon(options: RibbonOptions = {}): TableViewChrome {
    return new TableViewChrome(this, options);
  }

  /** Shows `dialog` pinned to the view: hidden while another view is current, closed when the view goes. What it
   * follows ends with the showing; an open dialog is left as it is. Answers `dialog`. */
  show(dialog: Dialog, options: Omit<DialogShowOptions, 'visible' | 'owner'> = {}): Dialog {
    if (dialog.isOpen.peek())
      return dialog;
    const view = this.view;
    const pin = this.scope.child('view-pin');
    const visible = signal(isView(grok.shell.v, view));
    pin.subscribe(grok.events.onCurrentViewChanged, () => visible.value = isView(grok.shell.v, view));
    // the dialog's own listener: nothing a builder disposes may cut the pin's end loose
    const off = dialog.build(() => dialog.on('close', () => {
      off();
      pin.dispose();
    }));
    dialog.show({...options, visible, owner: pin});
    if (!dialog.isOpen.peek()) {
      off();
      pin.dispose();
    }
    return dialog;
  }

  /** Docks `control` in the view at `place` (`'left'`, `'right'`, `'up'`, `'down'`, `'fill'`: `DG.DOCK_TYPE`) — in a
   * TableView's own dock manager, else in the shell's. It is disposed when its panel is closed (its ✕, or the dock
   * manager closing it from code) and when the view goes; disposing it closes its panel. */
  dock(control: Control, place: DG.DockType, options: ViewDockOptions = {}): void {
    const view = this.view;
    const manager = 'dockManager' in view ? (view as DG.TableView).dockManager : grok.shell.dockManager;
    const root = host(control);
    const holds = (el: HTMLElement) => el === root || el.contains(root);
    control.own(Ownership.adopt(this, () => control.dispose()));
    control.subscribe(manager.onClosed, (el) => {
      if (holds(el))
        control.dispose();
    });
    // closed from code, a panel is only hidden: one the layout no longer has is gone
    control.subscribe(manager.onPanelHidden, (el) => {
      if (holds(el) && manager.findNode(root) == null)
        control.dispose();
    });
    control.own(() => {
      if (manager.findNode(root) != null)
        manager.close(root);
    });
    manager.dock(root, place, null, options.title, options.ratio ?? 0.5);
  }

  protected get ambientAdoption(): boolean {
    return false;
  }

  private _attach(view: DG.ViewBase): this {
    const key = ViewHost._key(view);
    this._view = view;
    ViewHost._hosts.set(key, this);
    this.own(() => {
      if (ViewHost._hosts.get(key) === this)
        ViewHost._hosts.delete(key);
    });
    this.subscribe(grok.events.onViewRemoved, (removed) => {
      if (isView(removed, view))
        this.dispose();
    });
    return this;
  }

  private _chrome(view: DG.ViewBase, options: ViewContent): void {
    const place = (item: ChromeItem): HTMLElement => {
      if (!Control.is(item))
        return item;
      this.own(() => item.dispose());
      return item.root;
    };
    const panel = (item: ChromeItem): HTMLDivElement => {
      const el = place(item);
      if (el instanceof HTMLDivElement)
        return el;
      const wrap = document.createElement('div');
      wrap.append(el);
      return wrap;
    };
    if (options.ribbon) {
      const groups = options.ribbon.map((group) => group.map(place));
      view.setRibbonPanels(groups);
      // a panel whose items are all hidden goes, with its divider
      const fold = () => {
        for (const group of groups) {
          const panel = group[0]?.closest<HTMLElement>('.d4-ribbon-panel');
          if (panel == null)
            continue;
          panel.style.display = group.every((el) => el.hidden) ? 'none' : '';
          const divider = panel.previousElementSibling;
          if (divider instanceof HTMLElement && divider.classList.contains('d4-ribbon-separator'))
            divider.style.display = panel.style.display;
        }
      };
      const observer = new MutationObserver(fold);
      for (const el of groups.flat())
        observer.observe(el, {attributes: true, attributeFilter: ['hidden']});
      this.own(() => observer.disconnect());
      fold();
    }
    if (options.toolbox)
      view.toolbox = place(options.toolbox);
    const status = options.status;
    if (status instanceof Signal) {
      const el = document.createElement('div');
      this.effect(() => {
        el.textContent = status.value;
      });
      view.statusBarPanels = [el];
    }
    else if (status != null) {
      const items = (Array.isArray(status) ? status : [status]) as ChromeItem[];
      view.statusBarPanels = items.map(panel);
    }
    const path = options.path;
    if (typeof path === 'string')
      (view as DG.View).path = path;
    else if (path != null) {
      this.effect(() => {
        (view as DG.View).path = path.value;
      });
    }
  }

  private static _key(view: DG.ViewBase): object {
    return view.dart ?? view;
  }

  private static _isContent(x: BuiltRoot | ViewContent): x is ViewContent {
    return typeof x === 'object' && x !== null && !Array.isArray(x) && !(x instanceof HTMLElement) && !Control.is(x) &&
      'root' in x;
  }
}

/** `view`'s state as signals, its subscriptions on `scope`: what {@link ViewHost.signals} is made of. */
export function viewSignalsOn(view: DG.ViewBase, scope: Scope): ViewSignals {
  const table = 'dataFrame' in view ? view as DG.TableView : null;
  const isCurrent = signal(isView(grok.shell.v, view));
  const grid = signal<DG.Grid | null>(table?.grid ?? null);
  const dataFrame = signal<DG.DataFrame | null>(table?.dataFrame ?? null);
  scope.subscribe(grok.events.onCurrentViewChanged, () => isCurrent.value = isView(grok.shell.v, view));
  // a layout replaces the grid, and may hand the view another table
  scope.subscribe(grok.events.onViewLayoutApplied, (info) => {
    if (table === null || !isView(info?.view, view))
      return;
    if (table.grid?.dart !== grid.peek()?.dart)
      grid.value = table.grid ?? null;
    if (table.dataFrame?.dart !== dataFrame.peek()?.dart)
      dataFrame.value = table.dataFrame ?? null;
  });
  const frame = computed(() => dataFrame.value ?? undefined) as unknown as ReadonlySignal<DataFrameLike | undefined>;
  const steps = dfBindings(frame, scope);
  const step = <T>(name: string) => steps.bindStep(name) as unknown as T;
  return {isCurrent, grid, dataFrame, currentRowIdx: step('currentRowIdx'), selection: step('selection'),
    filter: step('filter'), columns: step('columns'), rowCount: step('rowCount')};
}
