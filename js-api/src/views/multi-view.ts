import {View} from "./view";
import {ObjectHandler} from "../../ui";
import {toJs} from "../wrappers";
import {TabControl} from "../widgets";
let api = <any>window;

class EmptyView extends View {
  constructor() {
    super(api.grok_View());
  }
}

export interface MultiViewOptions {
  viewFactories?: {[name: string]: () => View};
}

export class MultiView extends View {
  _views: Map<String, View> = new Map();
  _options?: MultiViewOptions;
  _currentView: View = new EmptyView();
  tabs: TabControl = TabControl.create();

  constructor(options?: MultiViewOptions) {
    super(api.grok_View());

    this._options = options;
    this.root.appendChild(this.tabs.root);
    this.tabs.onTabChanged.subscribe((_) => this.currentView = this.getView(this.tabs.currentPane.name));

    if (options?.viewFactories) {
      for (let [name, factory] of Object.entries(options.viewFactories))
        this.tabs.addPane(name, () => factory().root!);
    }
  }

  getView(name: string): View {
    if (!this._views.has(name))
      this._views.set(name, this._options?.viewFactories![name]()!);

    return this._views.get(name)!;
  }

  get currentView(): View { return this._currentView; }
  set currentView(x) {
    this._currentView = x;
    this.toolbox = x.toolbox;
    this.setRibbonPanels(x.getRibbonPanels());
  }
}