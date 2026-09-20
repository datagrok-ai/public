import {DockView, View, ViewBase} from "./view";
import {TabControl, TabPane} from "../widgets";
import {IDartApi} from "../api/grok_api.g";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

class EmptyView extends View {
  constructor() {
    super(api.grok_View());
  }
}

type ViewDescription = {factory: () => ViewBase, allowClose: boolean};
type ViewFactory = () => ViewBase;

export interface MultiViewOptions {
  viewFactories?: {[name: string]: ViewFactory | ViewDescription};
}

export class MultiView extends ViewBase {
  _views: Map<String, ViewBase> = new Map();
  _options: MultiViewOptions;
  _currentView: ViewBase = new EmptyView();
  tabs: TabControl = TabControl.create();
  private _fixedName: string | undefined;
  private _ribbonItems: WeakMap<HTMLElement, HTMLElement> = new WeakMap();
  private _activated: WeakSet<ViewBase> = new WeakSet();

  constructor(options?: MultiViewOptions) {
    super({});
    this.box = true;

    this._options = options ?? {viewFactories: {}};
    this._options.viewFactories ??= {};
    this.root.appendChild(this.tabs.root);
    this.tabs.onTabChanged.subscribe((_) => this.currentView = this.getView(this.tabs.currentPane.name));

    this.tabs.onTabRemoved.subscribe((tab) => {
      if (this._views.get(tab.name) == this._currentView)
        this.innerView = null;
      this._views.delete(tab.name);
      delete this._options.viewFactories![tab.name];
    });
    if (options?.viewFactories) {
      for (let [name] of Object.entries(options.viewFactories)) {
        this._addNewViewTab(name, false);
      }
    }
  }

  private _addNewViewTab(name: string, activate: boolean): TabPane {
    let allowClose = false;
    if ((<ViewDescription>(this._options?.viewFactories![name]!)).allowClose)
      allowClose = true;
    let tab = this.tabs.addPane(name, () => this._paneContent(name), null, {allowClose: allowClose});
    if (activate)
      this.tabs.currentPane = tab;
    return tab;
  }

  addView(name: string, desc: ViewFactory | ViewDescription, activate: boolean) {
    (<any>this._options?.viewFactories)![name] = desc;
    this._addNewViewTab(name, activate);
  }

  _getFactory(factory: ViewFactory | ViewDescription) {
    let _factory: ViewFactory = (<ViewDescription>factory).factory;
    if (_factory != undefined)
      factory = _factory;
    return factory;
  }

  getView(name: string): ViewBase {
    if (!this._views.has(name)) {
      let factory = this._getFactory(this._options?.viewFactories![name]!);
      this._views.set(name, (<ViewFactory>factory)()!);
    }
    return this._views.get(name)!;
  }

  get name(): string {
    return this._fixedName ?? super.name;
  }

  set name(s: string) {
    this._fixedName = s;
    super.name = s;
  }

  /** A pane child stretches only if it declares a ui layout class. A Dart view root
   * gets one from the shell when it is opened there. */
  private _paneContent(name: string): HTMLElement {
    let root = this.getView(name).root!;
    if (!root.classList.contains('ui-box') && !root.classList.contains('ui-panel') && !root.classList.contains('ui-div'))
      root.classList.add('ui-box');
    return root;
  }

  private _ribbonPanelsOf(x: ViewBase): HTMLElement[][] {
    return x.getRibbonPanels().map((panel) => panel.map((item) => {
      if (!item.classList.contains('d4-ribbon-item'))
        return item;
      let content = this._ribbonItems.get(item);
      if (content == null && item.firstElementChild instanceof HTMLElement) {
        content = item.firstElementChild as HTMLElement;
        this._ribbonItems.set(item, content);
      }
      return content ?? item;
    }));
  }

  get currentView(): ViewBase { return this._currentView; }
  set currentView(x) {
    this._currentView = x;
    this.innerView = x;
    this.toolbox = x.toolbox;
    this.setRibbonPanels(this._ribbonPanelsOf(x));
    this.ribbonMenu = x.ribbonMenu;
    this.statusBarPanels = x.statusBarPanels;
    this.name = x.name;

    if (x instanceof View && !this._activated.has(x)) {
      this._activated.add(x);
      x._onAdded();
    }

    if (x instanceof DockView)
      x._handleResize();
  }
}
