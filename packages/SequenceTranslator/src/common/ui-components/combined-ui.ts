/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TAB, APP} from './const';
import {TranslatorLayoutHandler} from '../../translator-app/view/app-ui';
import {StructureLayoutHandler} from '../../structure-app/view/app-ui';
import {PatternLayoutController} from '../../pattern-app/view/layout';
import {MonomerLibViewer} from '../monomer-lib/viewer';
import {_package} from '../../package';
import {tryCatch} from '../model/helpers';

type ViewFactories = {[name: string]: () => DG.View};

export abstract class AppUIBase {
  constructor(protected appName: string, protected parentAppName?: string) { }
  abstract addView(): Promise<void>;

  async createAppLayout(): Promise<void> {
    const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create(`Loading ${this.appName}...`);

    let currentView = grok.shell.v?.root;
    if (currentView)
      ui.setUpdateIndicator(currentView, true);

    await tryCatch(async () => {
      await this.addView();
    }, () => pi.close());

    if (currentView)
      ui.setUpdateIndicator(currentView, false);
  }
}

abstract class SimpleAppUIBase extends AppUIBase {
  constructor(appName: string) {
    super(appName);
    this.view = DG.View.create();
    this.setupView();
  }

  protected view: DG.View;
  async addView(): Promise<void> {
    await this.initView();
    const name = this.parentAppName ? this.parentAppName + '/' + this.appName : this.appName;
    this.view.path = `/apps/${_package.name}/${name.replace(/\s/g, '')}/`;
    grok.shell.addView(this.view);
  }

  protected abstract getHtml(): Promise<HTMLDivElement>;
  async initView(): Promise<void> {
    const html = await this.getHtml();
    this.view.append(html);
  }

  protected setupView(): void {
    this.view.box = true;
    this.view.name = this.appName;

    const windows = grok.shell.windows;
    windows.showProperties = false;
    windows.showToolbox = false;
    windows.showHelp = false;
  }

  getView(): DG.View {
    return this.view;
  }
}

export class CombinedAppUI extends AppUIBase {
  constructor(externalViewFactories: ViewFactories) {
    super(APP.COMBINED)
    this.externalViewFactories = externalViewFactories;
    const factories = this.getViewFactories();
    this.multiView = new DG.MultiView({viewFactories: factories}); 
  }

  private multiView: DG.MultiView;
  private externalViewFactories?: ViewFactories;

  private getViewFactories(): ViewFactories {
    function viewFactory(uiConstructor: new (view: DG.View) => SimpleAppUIBase): () => DG.View {
      const view = DG.View.create();
      const translateUI = new uiConstructor(view);
      // intentonally don't await for the promise
      translateUI.initView();
      return () => translateUI.getView();
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

  async addView(): Promise<void> {
    this.multiView.tabs.onTabChanged.subscribe(() => this.setUrl());
    this.setUrl();
    grok.shell.addView(this.multiView);
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

export class AppUIFactory {
  private constructor() {}

  static getUI(appName: string): SimpleAppUIBase {
    switch (appName) {
      case APP.TRANSLATOR:
        return new OligoTranslatorUI();
      case APP.PATTERN:
        return new OligoPatternUI();
      case APP.STRUCTRE:
        return new OligoStructureUI();
      default:
        throw new Error(`Unknown app name: ${appName}`);
    }
  }
}

class OligoTranslatorUI extends SimpleAppUIBase {
  constructor() {
    super(APP.TRANSLATOR);

    const viewMonomerLibIcon = ui.iconFA('book', MonomerLibViewer.view, 'View monomer library');
    this.topPanel = [
      viewMonomerLibIcon,
    ];
    this.view.setRibbonPanels([this.topPanel]);
    this.ui = new TranslatorLayoutHandler();
  }

  private readonly topPanel: HTMLElement[];
  private readonly ui: TranslatorLayoutHandler;

  protected getHtml(): Promise<HTMLDivElement> {
    return this.ui.getHtmlElement();
  };
}

class OligoPatternUI extends SimpleAppUIBase {
  constructor() {
    super(APP.PATTERN);
    this.ui = new PatternLayoutController();
  }
  private readonly ui: PatternLayoutController;
  protected getHtml(): Promise<HTMLDivElement> {
    return Promise.resolve(this.ui.htmlDivElement);
  }
}

class OligoStructureUI extends SimpleAppUIBase {
  constructor() {
    super(APP.STRUCTRE)
    this.ui = new StructureLayoutHandler();
  }
  private readonly ui: StructureLayoutHandler;

  protected getHtml(): Promise<HTMLDivElement> {
    return this.ui.getHtmlDivElement();
  }
}
