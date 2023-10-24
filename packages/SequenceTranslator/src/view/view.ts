/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TAB, APP} from './const/view';
import {MainTabUI} from './tabs/main';
import {SdfTabUI} from './tabs/sdf';
import {AxolabsTabUI} from './tabs/axolabs';
import {MonomerLibViewer} from './monomer-lib-viewer/viewer';
import {_package} from '../package';

type ViewFactories = {[name: string]: () => DG.View};

export class AppMultiView {
  constructor(externalViewFactories?: ViewFactories) {
    this.externalViewFactories = externalViewFactories;
    const factories = this.getViewFactories();
    console.log(`factories:`, factories);
    this.multiView = new DG.MultiView({viewFactories: factories}); 
  }

  private multiView: DG.MultiView;
  private externalViewFactories?: ViewFactories;

  private getViewFactories(): ViewFactories {
    function viewFactory(uiConstructor: new (view: DG.View) => AppUIBase): () => DG.View {
      const view = DG.View.create();
      const translateUI = new uiConstructor(view);
      // intentonally don't await for the promise
      translateUI.initView();
      return () => translateUI.appView;
    }

    let result: {[key: string]: () => DG.View } = {
      [TAB.TRANSLATOR]: viewFactory(OligoTranslatorUI),
      [TAB.PATTERN]: viewFactory(OligoPatternUI),
      [TAB.STRUCTRE]: viewFactory(OligoStructureUI),
    }

    if (this.externalViewFactories)
      result = Object.assign({}, result, this.externalViewFactories);

    return result
  }

  private getPath(): string {
    let name = this.multiView.tabs.currentPane.name;
    name = name.charAt(0).toUpperCase() + name.substring(1).toLowerCase();
    const path = `/apps/${_package.name}/OligoToolkit/${name}`;
    return path;
  }

  private setUrl(): void {
    this.multiView.path = this.getPath();
  }

  async createLayout(): Promise<void> {
    this.multiView.tabs.onTabChanged.subscribe(() => this.setUrl());
    this.setUrl();
    grok.shell.addView(this.multiView);
  }
}

export abstract class AppUIBase {
  constructor(private appName: string, private parentAppName?: string) {
  }

  protected view: DG.View;
  protected abstract getHtml(): Promise<HTMLDivElement>;

  async initView(): Promise<void> {
    const html = await this.getHtml();
    this.view.append(html);
  }

  /** Create master layout of the app  */
  async createLayout(): Promise<void> {
    await this.initView();
    const name = this.parentAppName ? this.parentAppName + '/' + this.appName : this.appName;
    this.view.path = `/apps/${_package.name}/${name.replace(/\s/g, '')}/`;
    grok.shell.addView(this.view);
  }
  
  protected setupView(): void {
    this.view.box = true;
    this.view.name = this.appName;

    const windows = grok.shell.windows;
    windows.showProperties = false;
    windows.showToolbox = false;
    windows.showHelp = false;
  }

  get appView(): DG.View {
    return this.view;
  }
}

abstract class SimpleAppUIBase extends AppUIBase {
  constructor(appName: string) {
    super(appName);
    this.view = DG.View.create();
    this.setupView();
  }
}

/** For plugins from external packages */
export class ExternalPluginUI extends SimpleAppUIBase {
  constructor(viewName: string, layout: HTMLDivElement) {
    super(viewName);
    this.layout = layout;
  }
  private layout: HTMLDivElement;

  protected getHtml(): Promise<HTMLDivElement> {
    return Promise.resolve(this.layout);
  }
}

export class OligoTranslatorUI extends SimpleAppUIBase {
  constructor() {
    super(APP.TRANSLATOR);

    const viewMonomerLibIcon = ui.iconFA('book', MonomerLibViewer.view, 'View monomer library');
    this.topPanel = [
      viewMonomerLibIcon,
    ];
    this.view.setRibbonPanels([this.topPanel]);
    this.ui = new MainTabUI();
  }

  private readonly topPanel: HTMLElement[];
  private readonly ui: MainTabUI;

  protected getHtml(): Promise<HTMLDivElement> {
    return this.ui.getHtmlElement();
  };
}

export class OligoPatternUI extends SimpleAppUIBase {
  constructor() {
    super(APP.PATTERN);
    this.ui = new AxolabsTabUI();
  }
  private readonly ui: AxolabsTabUI;
  protected getHtml(): Promise<HTMLDivElement> {
    return Promise.resolve(this.ui.htmlDivElement);
  }
}

export class OligoStructureUI extends SimpleAppUIBase {
  constructor() {
    super(APP.STRUCTRE)
    this.ui = new SdfTabUI();
  }
  private readonly ui: SdfTabUI;

  protected getHtml(): Promise<HTMLDivElement> {
    return this.ui.getHtmlDivElement();
  }
}
