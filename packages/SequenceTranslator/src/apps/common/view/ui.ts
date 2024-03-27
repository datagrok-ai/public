/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {_package} from '../../../package';
import {PatternAppLayout} from '../../pattern/view/ui';
import {StructureAppLayout} from '../../structure/view/app-ui';
import {TranslatorAppLayout} from '../../translator/view/app-ui';
import {tryCatch} from '../model/helpers';
import {MonomerLibViewer} from './monomer-lib-viewer';
import {APP_NAME, TAB_NAME} from './const';

type ViewFactories = {[name: string]: () => DG.View};

export abstract class AppUIBase {
  constructor(protected appName: string, protected parentAppName?: string) { }
  abstract addView(): Promise<void>;

  async initializeAppLayout(): Promise<void> {
    const progressIndicator: DG.TaskBarProgressIndicator =
      DG.TaskBarProgressIndicator.create(`Loading ${this.appName}...`);

    const currentView = grok.shell.v?.root;
    if (currentView)
      ui.setUpdateIndicator(currentView, true);

    await tryCatch(async () => {
      await this.addView();
    }, () => progressIndicator.close());

    if (currentView)
      ui.setUpdateIndicator(currentView, false);
  }
}

abstract class IsolatedAppUIBase extends AppUIBase {
  constructor(appName: string) {
    super(appName);
    this.view = DG.View.create();
    this.configureView();
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

  protected configureView(): void {
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
    super(APP_NAME.COMBINED);
    this.externalViewFactories = externalViewFactories;
    const factories = this.getViewFactories();
    this.multiView = new DG.MultiView({viewFactories: factories});
  }

  private multiView: DG.MultiView;
  private externalViewFactories?: ViewFactories;

  private getViewFactories(): ViewFactories {
    function viewFactory(UiConstructor: new (view: DG.View) => IsolatedAppUIBase): () => DG.View {
      const view = DG.View.create();
      const translateUI = new UiConstructor(view);
      translateUI.initView().catch(
        (e) => console.error(`Failed to initialize ${UiConstructor.name}: ${e}`)
      );
      return () => translateUI.getView();
    }

    let result: {[key: string]: () => DG.View } = {
      [TAB_NAME.TRANSLATOR]: viewFactory(OligoTranslatorUI),
      [TAB_NAME.PATTERN]: viewFactory(OligoPatternUI),
      [TAB_NAME.STRUCTURE]: viewFactory(OligoStructureUI),
    };

    if (this.externalViewFactories)
      result = Object.assign({}, result, this.externalViewFactories);

    return result;
  }

  private getCurrentPanePath(): string {
    let name = this.multiView.tabs.currentPane.name;
    name = name.charAt(0).toUpperCase() + name.substring(1).toLowerCase();
    const path = `/apps/${_package.name}/OligoToolkit/${name}`;
    return path;
  }

  private setUrl(): void {
    this.multiView.path = this.getCurrentPanePath();
  }

  async addView(): Promise<void> {
    this.multiView.tabs.onTabChanged.subscribe(() => this.setUrl());
    this.setUrl();
    grok.shell.addView(this.multiView);
  }
}

/** For plugins from external packages */
export class ExternalPluginUI extends IsolatedAppUIBase {
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

  static createAppUIInstance(appName: string): IsolatedAppUIBase {
    switch (appName) {
    case APP_NAME.TRANSLATOR:
      return new OligoTranslatorUI();
    case APP_NAME.PATTERN:
      return new OligoPatternUI();
    case APP_NAME.STRUCTURE:
      return new OligoStructureUI();
    default:
      throw new Error(`Unknown app name: ${appName}`);
    }
  }
}

class OligoTranslatorUI extends IsolatedAppUIBase {
  private readonly topPanel: HTMLElement[];
  private readonly layout = new TranslatorAppLayout();

  constructor() {
    super(APP_NAME.TRANSLATOR);

    const viewMonomerLibIcon = ui.iconFA('book', MonomerLibViewer.view, 'View monomer library');
    this.topPanel = [
      viewMonomerLibIcon,
    ];
    this.view.setRibbonPanels([this.topPanel]);
  }

  protected getHtml(): Promise<HTMLDivElement> {
    return this.layout.getHtmlElement();
  };
}

class OligoPatternUI extends IsolatedAppUIBase {
  constructor() {
    super(APP_NAME.PATTERN);
  }
  protected getHtml(): Promise<HTMLDivElement> {
    return PatternAppLayout.generateHTML();
  }
}

class OligoStructureUI extends IsolatedAppUIBase {
  constructor() {
    super(APP_NAME.STRUCTURE);
    this.layout = new StructureAppLayout();
  }
  private readonly layout: StructureAppLayout;

  protected getHtml(): Promise<HTMLDivElement> {
    return this.layout.getHtmlDivElement();
  }
}
