/* ViewHost.ribbon: a package's chrome over a platform TableView, owned by the view. Its ribbon group sits after
   a rule at the ribbon's end, pinned there while the ribbon overflows; the small icons fold into one "…" menu
   when the ribbon cannot hold them, the primary buttons always show. Plus the crumb back to the app (the shell
   renders it: `TableView.showParentCrumb`), the project Save left off (`TableView.showSaveButton`), and a status
   line. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import type * as DG from 'datagrok-api/dg';
import {Control} from '../../core/component.js';
import {Ownership} from '../../core/owner.js';
import {signal, computed, Signal, ReadonlySignal} from '../../core/signals.js';
import {button, div, span} from '../../core/elements.js';
import type {Reason} from '../../core/reason.js';
import type {LiveOption} from '../../core/input-base.js';
import {Input} from '../../core/input-base.js';
import type {IItemsState} from '../../core/widget-like.js';
import {iconButton} from '../../components/actions/buttons.js';
import {Menu} from '../../components/navigation/menu.js';
import {Toolbar} from '../../components/navigation/toolbar.js';
import {isView} from './host.js';
import type {ViewHost} from './view-host.js';

/** One of the package's small ribbon icons. */
export interface RibbonAction {
  /** A Font Awesome name. */
  icon: string;
  /** The icon's tooltip and accessible name ("Refresh Campaign"). */
  tooltip: string;
  /** Its item in the "…" menu: the tooltip by default; end it in "..." when a dialog follows. */
  text?: string;
  /** Returns a promise or nothing, as an {@link Action}'s `run`. */
  run: () => void | Promise<void>;
  /** Offered only while true (a signal is followed): no edit right, no layout to apply. */
  visible?: LiveOption<boolean>;
  /** A note while non-null ("New changes available. Click to refresh campaign."): the icon pulses and says it
   * as its tooltip; folded, the "…" button carries a dot and the menu item `hintDetail`. */
  hint?: ReadonlySignal<string | null>;
  /** The menu item's detail while hinted ("New changes"). */
  hintDetail?: string;
  /** The automation id: its `data-u2-name` and its item key in the group's state (the tooltip by default). */
  name?: string;
}

/** One of the group's buttons — shown whatever the width, after the icons. */
export interface RibbonButton {
  text: string;
  /** Returns a promise or nothing, as an {@link Action}'s `run`. */
  run: () => void | Promise<void>;
  /** The ribbon's raised button (default true). */
  primary?: boolean;
  /** Why it is off (its tooltip), `null` while on; a signal is followed. */
  disabled?: Reason;
  visible?: LiveOption<boolean>;
  /** A class for the button element (a selector older code and tests know it by). */
  cls?: string;
  name?: string;
}

/** The crumb before the view's own name ("Hit Design › KHL-3"), which the shell renders for the view's parent. */
export interface ChromeCrumbs {
  /** The view it goes back to. */
  parent: DG.ViewBase;
  /** What a press on it does besides making the parent current (the shell's part): close the view, say. */
  onClick?: () => void;
}

/** What {@link ViewHost.ribbon} puts over the view. */
export interface RibbonOptions {
  /** The small icons, in order. */
  actions?: RibbonAction[];
  /** The buttons after them. */
  buttons?: RibbonButton[];
  /** Fold the icons into "…": when the ribbon cannot hold them (`'auto'`, the default), always or never. */
  fold?: 'auto' | 'always' | 'never';
  /** Keep the group at the ribbon's visible end while the ribbon overflows (default true). */
  pin?: boolean;
  /** Leave the project's Save off the ribbon: for an app that saves on its own. */
  hideSave?: boolean;
  breadcrumbs?: ChromeCrumbs;
  /** A line in the view's status bar. */
  status?: ReadonlySignal<string>;
  /** The group's automation id (`data-u2-name`). */
  name?: string;
  /** The group's accessible name: 'Actions' by default. */
  label?: string;
}

interface Placed<T> {
  spec: T;
  el: HTMLButtonElement;
  shown: Signal<boolean>;
}

/** A package's ribbon group, crumb, Save policy and status line over a platform TableView, owned by the view's
 * {@link ViewHost} (`host.ribbon(options)` builds it): disposed with the view, or before it by `dispose()`. Its state
 * block lists the icons and buttons as items (`rendered: false` once folded) and says `extra.folded`. */
export class TableViewChrome extends Control {
  /** The icons are folded into "…". */
  readonly folded: ReadonlySignal<boolean>;
  readonly view: DG.TableView;

  private readonly _options: RibbonOptions;
  private readonly _folded = signal(false);
  private readonly _icons = span('', 'u2-ribbon-icons');
  private readonly _actions: Placed<RibbonAction>[] = [];
  private readonly _buttons: Placed<RibbonButton>[] = [];
  private _more!: HTMLButtonElement;
  private _ribbon: HTMLElement | undefined;
  private _unlistenCrumb: (() => void) | undefined;
  private _panel: HTMLElement | undefined;
  private _queued = false;
  private _status: HTMLDivElement | undefined;
  /** The view's Save and crumb as the chrome found them: what it puts back when it goes before the view. */
  private _found: {showSaveButton: boolean, parentView: DG.ViewBase, showParentCrumb: boolean} | undefined;
  private readonly _ribbonResize = new ResizeObserver(() => this._queueRefit());
  /** The platform rebuilds a view's ribbon panels (twice in the first 2 s on the stand): the pin goes back on
   * before the next paint, the fold is measured again. */
  private readonly _ribbonTree = new MutationObserver(() => {
    this._attach();
    this._queueRefit();
  });

  constructor(host: ViewHost, options: RibbonOptions = {}) {
    super(document.createElement('div'));
    const view = this.view = host.view as DG.TableView;
    this._options = options;
    this.folded = this._folded;
    const root = this.root;
    root.className = 'u2-ribbon-group d4-ribbon-wrapper';
    root.dataset.u2 = 'ribbon-group';
    root.setAttribute('role', 'toolbar');
    root.setAttribute('aria-label', options.label ?? 'Actions');
    if (options.name !== undefined)
      this.name = options.name;
    this.build(() => this._build());

    const own = new ResizeObserver(() => this._queueRefit());
    own.observe(root);
    this.own(() => {
      own.disconnect();
      this._ribbonResize.disconnect();
      this._ribbonTree.disconnect();
    });
    this.own(Ownership.adopt(host, () => this.dispose()));
    this.subscribe(grok.events.onCurrentViewChanged, () => {
      if (isView(grok.shell.v, view))
        this._queueRefit();
    });
    this.own(() => this._remove());
    this._place();
  }

  /** Measures the ribbon now and folds or unfolds the icons; it runs by itself on every resize. */
  refit(): void {
    this._attach();
    const ribbon = this._ribbon;
    const fold = this._options.fold ?? 'auto';
    const any = this._actions.some((a) => a.shown.peek());
    if (fold !== 'auto' || ribbon === undefined || ribbon.clientWidth === 0) {
      this._setFolded(fold === 'always' && any);
      return;
    }
    this._setFolded(false);
    this._setFolded(any && ribbon.scrollWidth > ribbon.clientWidth + 1);
    // the platform's icons run under the pinned group: an edge shows they go on
    if (this._panel !== undefined && this._options.pin !== false)
      this._panel.classList.toggle('u2-ribbon-overflow', ribbon.scrollWidth > ribbon.clientWidth + 1);
  }

  /** The icons, then the buttons; a folded icon is not `rendered`, a disabled button `disabled`. */
  getItems(opts?: {all?: boolean}): IItemsState {
    const folded = this._folded.peek();
    const list = [
      ...this._actions.filter((a) => a.shown.peek()).map((a) =>
        ({key: TableViewChrome._key(a.spec.name, a.spec.tooltip), label: a.spec.tooltip, rendered: !folded})),
      ...this._buttons.filter((b) => b.shown.peek()).map((b) => ({key: TableViewChrome._key(b.spec.name, b.spec.text),
        label: b.spec.text, rendered: true, ...b.el.disabled ? {disabled: true} : {}})),
    ].map((item, index) => ({...item, index}));
    return Input.itemsState(list, opts?.all);
  }

  getState(): ReturnType<Control['getState']> {
    return {...super.getState(), extra: {folded: this._folded.peek()}};
  }

  private _build(): void {
    const options = this._options;
    const rule = span('', 'u2-ribbon-rule');
    rule.setAttribute('aria-hidden', 'true');
    for (const action of options.actions ?? []) {
      const hint = action.hint;
      const el = iconButton(action.icon, () => action.run(),
        {tooltip: hint === undefined ? action.tooltip : computed(() => hint.value ?? action.tooltip)});
      el.dataset.u2Item = TableViewChrome._key(action.name, action.tooltip);
      if (action.name !== undefined)
        el.dataset.u2Name = action.name;
      const placed = {spec: action, el, shown: signal(true)};
      this._live(action.visible, (on) => {
        el.hidden = !on;
        placed.shown.value = on;
      });
      if (hint !== undefined) {
        this.effect(() => {
          const note = hint.value;
          el.classList.toggle('u2-ribbon-hinted', note !== null);
          if (note === null || this._folded.peek() || !el.isConnected)
            return;
          const target = ui.hints.addHintIndicator(el, true, 10000);
          return () => ui.hints.remove(target);
        });
      }
      this._icons.append(el);
      this._actions.push(placed);
    }
    const hinted = () => this._actions.find((a) => a.shown.value && a.spec.hint?.value != null);
    this._more = iconButton(Toolbar.OVERFLOW_ICON, () => {}, {
      tooltip: computed(() => {
        const action = hinted();
        return action === undefined ? Toolbar.OVERFLOW_LABEL :
          `${Toolbar.OVERFLOW_LABEL} — ${action.spec.hintDetail ?? action.spec.hint!.peek()}`;
      }),
      menu: () => {
        const menu = Menu.popup();
        for (const a of this._actions.filter((x) => x.shown.peek())) {
          const note = a.spec.hint?.peek() ?? null;
          menu.item(a.spec.text ?? a.spec.tooltip, () => a.spec.run(),
            {icon: a.spec.icon, ...note === null ? {} : {detail: a.spec.hintDetail ?? note}});
        }
        return menu;
      },
    });
    this._more.classList.add('u2-ribbon-more');
    this._more.dataset.u2Part = 'more';
    this._more.hidden = true;
    const dot = span('', 'u2-ribbon-hint-dot');
    this._more.append(dot);
    this.effect(() => {
      dot.hidden = hinted() === undefined;
    });
    this.root.append(rule, this._icons, this._more);
    for (const spec of options.buttons ?? []) {
      const el = button(spec.text, () => spec.run(), {primary: spec.primary ?? true, disabled: spec.disabled});
      el.dataset.u2Item = TableViewChrome._key(spec.name, spec.text);
      if (spec.name !== undefined)
        el.dataset.u2Name = spec.name;
      if (spec.cls !== undefined)
        el.classList.add(spec.cls);
      const placed = {spec, el, shown: signal(true)};
      this._live(spec.visible, (on) => {
        el.hidden = !on;
        placed.shown.value = on;
      });
      this.root.append(el);
      this._buttons.push(placed);
    }
    const status = options.status;
    if (status !== undefined)
      this._status = div([status], 'u2-ribbon-status');
  }

  /** Into the ribbon, after the platform's panels: one panel holding the group. */
  private _place(): void {
    const view = this.view;
    this._found = {showSaveButton: view.showSaveButton, parentView: view.parentView,
      showParentCrumb: view.showParentCrumb};
    if (this._options.hideSave === true)
      view.showSaveButton = false;
    view.setRibbonPanels([[this.root]], false);
    const crumbs = this._options.breadcrumbs;
    if (crumbs !== undefined) {
      view.parentView = crumbs.parent;
      view.showParentCrumb = true;
    }
    if (this._status !== undefined)
      view.statusBarPanels = [...view.statusBarPanels ?? [], this._status];
    this._queueRefit();
  }

  /** Off the ribbon again, the platform's Save and crumb as they were, when the chrome goes before its view. The view's
   * removal brings it here through the host with the view already out of the shell: nothing is written to a closed
   * view. */
  private _remove(): void {
    const view = this.view;
    if (![...grok.shell.views].some((v) => isView(v, view)))
      return;
    const root = this.root;
    const panels = view.getRibbonPanels();
    if (panels.some((group) => group.some((el) => el === root || el.contains(root))))
      view.setRibbonPanels(panels.filter((group) => !group.some((el) => el === root || el.contains(root))));
    const found = this._found;
    if (found !== undefined && this._options.hideSave === true)
      view.showSaveButton = found.showSaveButton;
    if (found !== undefined && this._options.breadcrumbs !== undefined) {
      view.parentView = found.parentView;
      view.showParentCrumb = found.showParentCrumb;
    }
    if (this._status !== undefined)
      view.statusBarPanels = (view.statusBarPanels ?? []).filter((el) => el !== this._status);
  }

  private _queueRefit(): void {
    if (this._queued || this.scope.isDisposed)
      return;
    this._queued = true;
    // a frame later: folding inside the observer's own callback would resize what it observes
    this.scope.requestFrame(() => {
      this._queued = false;
      this.refit();
    });
  }

  /** Finds the view's ribbon once the shell has put the group in it: pins the group's panel, follows the
   * ribbon's width and hears the press of the parent crumb the shell puts in it. */
  private _attach(): void {
    let ribbon = this.root.closest<HTMLElement>('.d4-ribbon') ?? undefined;
    const outer = (el: HTMLElement | undefined) => el?.parentElement?.closest<HTMLElement>('.d4-ribbon') ?? undefined;
    for (let up = outer(ribbon); up !== undefined; up = outer(up))
      ribbon = up;
    if (ribbon !== this._ribbon) {
      this._ribbonResize.disconnect();
      this._ribbonTree.disconnect();
      this._unlistenCrumb?.();
      this._ribbon = ribbon;
      if (ribbon !== undefined) {
        this._ribbonResize.observe(ribbon);
        this._ribbonTree.observe(ribbon, {childList: true, subtree: true});
        const onClick = this._options.breadcrumbs?.onClick;
        if (onClick !== undefined) {
          // the shell's own listener on the crumb makes the parent current and renders that view's ribbon, which
          // takes the crumb out of this one: the press is told on its way down, acted on once the shell is done —
          // a microtask later, since the shell's handler is async: a real click runs it between the listeners, a
          // scripted one only after the dispatch
          let pressed = false;
          const offDown = this.scope.listen(ribbon, 'click', (e) => {
            pressed = (e.target as Element | null)?.closest?.(TableViewChrome.PARENT_CRUMB) != null;
          }, {capture: true});
          const offUp = this.scope.listen(ribbon, 'click', () => {
            if (pressed) {
              pressed = false;
              queueMicrotask(() => {
                try {
                  onClick();
                }
                catch (e) {
                  this.scope.reportError(e);
                }
              });
            }
          });
          this._unlistenCrumb = () => {
            offDown();
            offUp();
          };
        }
      }
    }
    const panel = this.root.closest<HTMLElement>('.d4-ribbon-panel') ?? undefined;
    if (panel !== undefined && panel !== this._panel && this._options.pin !== false) {
      panel.classList.add('u2-ribbon-pinned');
      this._panel = panel;
    }
  }

  private _setFolded(on: boolean): void {
    this._icons.hidden = on;
    this._more.hidden = !on;
    if (this._folded.peek() !== on)
      this._folded.value = on;
  }

  private _live(option: LiveOption<boolean> | undefined, apply: (on: boolean) => void): void {
    if (option instanceof Signal) {
      this.effect(() => {
        apply(option.value);
        this._queueRefit();
      });
    }
    else if (option !== undefined)
      apply(option as boolean);
  }

  /** The crumb the shell renders for the current view's parent: before the view's name, in its ribbon. */
  private static readonly PARENT_CRUMB = '.d4-view-ribbon > .d4-ribbon-name:not(.grok-ns-entity-name)';

  private static _key(name: string | undefined, text: string): string {
    return name ?? text;
  }
}
