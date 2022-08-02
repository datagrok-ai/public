import {DockView, View, ViewBase} from "./view";
import {TabControl, TabPane} from "../widgets";
let api = <any>window;

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
  _options?: MultiViewOptions;
  _currentView: ViewBase = new EmptyView();
  tabs: TabControl = TabControl.create();
  private _fixedName: string | undefined;

  constructor(options?: MultiViewOptions) {
    super({});
    this.box = true;

    this._options = options;
    this.root.appendChild(this.tabs.root);
    this.tabs.onTabChanged.subscribe((_) => this.currentView = this.getView(this.tabs.currentPane.name));

    this.tabs.onTabRemoved.subscribe((tab) => {
      this._views.delete(tab.name);
      (<any>this._options!.viewFactories!)[tab.name] = undefined;
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
    let tab = this.tabs.addPane(name, () => this.getView(name).root!, null, {allowClose: allowClose});
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
    let view = this._views.get(name)!;
    return view;
  }

  get name(): string {
    return this._fixedName ?? super.name;
  }

  set name(s: string) {
    this._fixedName = s;
  }

  get currentView(): ViewBase { return this._currentView; }
  set currentView(x) {
    this._currentView = x;
    this.toolbox = x.toolbox;
    this.setRibbonPanels(x.getRibbonPanels());
    this.name = x.name;
    if (x instanceof DockView) {
      console.log('bingo');
      x.initDock();
      x._onAdded();
      x._handleResize();
    }
  }
}