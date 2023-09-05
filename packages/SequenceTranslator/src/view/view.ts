/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TRANSLATOR_TAB, PATTERN_TAB, STRUCTRE_TAB, STRUCTRE_APP_NAME, PATTERN_APP_NAME, TRANSLATOR_APP_NAME} from './const/view';
import {MainTabUI} from './tabs/main';
import {SdfTabUI} from './tabs/sdf';
import {AxolabsTabUI} from './tabs/axolabs';
import {MonomerLibViewer} from './monomer-lib-viewer/viewer';

type ViewFactories = {[name: string]: () => DG.View};

export class AppMultiView {
  constructor(externalViewFactories?: ViewFactories) {
    this.externalViewFactories = externalViewFactories;
    this.multiView = new DG.MultiView({viewFactories: this.viewFactories});
  }

  private multiView: DG.MultiView;
  private externalViewFactories?: ViewFactories;

  get viewFactories() {
    function viewFactory(uiConstructor: new (view: DG.View) => AppUI): () => DG.View {
      const view = DG.View.create();
      const translateUI = new uiConstructor(view);
      // intentonally don't await for the promise
      translateUI.initView();
      return () => view;
    }

    let result: {[key: string]: () => DG.View } = {
      [TRANSLATOR_TAB]: viewFactory(OligoTranslatorUI),
      [PATTERN_TAB]: viewFactory(OligoPatternUI),
      [STRUCTRE_TAB]: viewFactory(OligoStructureUI),
    }

    if (this.externalViewFactories)
      result = Object.assign({}, result, this.externalViewFactories);

    return result
  }

  async createLayout(): Promise<void> {
    grok.shell.addView(this.multiView);
  }
}

export abstract class AppUI {
  constructor(view: DG.View, viewName: string) {
    this.view = view;
    this.view.box = true;
    this.view.name = viewName;

    const windows = grok.shell.windows;
    windows.showProperties = false;
    windows.showToolbox = false;
    windows.showHelp = false;
  }

  protected readonly view: DG.View;
  protected abstract getHtml(): Promise<HTMLDivElement>;

  async initView(): Promise<void> {
    const html = await this.getHtml();
    this.view.append(html);
  }

  /** Create master layout of the app  */
  async createLayout(): Promise<void> {
    await this.initView();
    grok.shell.addView(this.view);
  }
}

/** For plugins from external packages */
export class ExternalPluginUI extends AppUI {
  constructor(view: DG.View, viewName: string, layout: HTMLDivElement) {
    super(view, viewName);
    this.layout = layout;
  }
  private layout: HTMLDivElement;

  protected getHtml(): Promise<HTMLDivElement> {
    return Promise.resolve(this.layout);
  }
}

export class OligoTranslatorUI extends AppUI {
  constructor(view: DG.View) {
    super(view, TRANSLATOR_APP_NAME);

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

export class OligoPatternUI extends AppUI {
  constructor(view: DG.View) {
    super(view, PATTERN_APP_NAME);
    this.ui = new AxolabsTabUI();
  }
  private readonly ui: AxolabsTabUI;
  protected getHtml(): Promise<HTMLDivElement> {
    return Promise.resolve(this.ui.htmlDivElement);
  }
}

export class OligoStructureUI extends AppUI {
  constructor(view: DG.View) {
    super(view,  STRUCTRE_APP_NAME)
    this.ui = new SdfTabUI();
  }
  private readonly ui: SdfTabUI;

  protected getHtml(): Promise<HTMLDivElement> {
    return this.ui.getHtmlDivElement();
  }
}

/** For proprietary plugins  */
export class ExternalPlugin extends AppUI {
  constructor(view: DG.View, app_name: string, getHtml: () => Promise<HTMLElement>) {
    super(view,  app_name)
  }

  protected getHtml(): Promise<HTMLDivElement> {
    return this.getHtml();
  }
}

// class URLRouter {
//   constructor(private readonly view: DG.View) {
//     this.pathParts = window.location.pathname.split('/');
//     this.searchParams = new URLSearchParams(window.location.search);
//   }

//   searchParams: URLSearchParams;
//   private pathParts: string[];

//   get urlParamsString(): string {
//     return Object.entries(this.searchParams)
//       .map(([key, value]) => `${key}=${encodeURIComponent(value)}`).join('&');
//   }

//   updatePath(control: DG.TabControl) {
//     const urlParamsTxt: string = Object.entries(this.searchParams)
//       .map(([key, value]) => `${key}=${encodeURIComponent(value)}`).join('&');
//     this.view.path = '/apps/SequenceTranslator' + `/${control.currentPane.name}`;
//     if (urlParamsTxt)
//       this.view.path += `/?${urlParamsTxt}`;
//   }

//   get tabName(): string {
//     const idx = this.pathParts.findIndex((el) => el === 'SequenceTranslator');
//     if (idx === -1) // todo: remove after verification of validity condition
//       return '';
//     return this.pathParts[idx + 1];
//   }
// }
