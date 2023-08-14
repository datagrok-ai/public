/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TRANSLATION_TAB, AXOLABS_TAB, DUPLEX_TAB} from './const/view';
import {MainTabUI} from './tabs/main';
import {SdfTabUI} from './tabs/sdf';
import {AxolabsTabUI} from './tabs/axolabs';
import {MonomerLibViewer} from './monomer-lib-viewer/viewer';

export class AppMultiView {
  constructor() {
    this.multiView = new DG.MultiView({viewFactories: this.viewFactories});
  }

  private multiView: DG.MultiView;

  get viewFactories() {
    function viewFactory(uiConstructor: new (view: DG.View) => AppUI): () => DG.View {
      const view = DG.View.create();
      const translateUI = new uiConstructor(view);
      // intentonally don't await for the promise
      translateUI.initView();
      return () => view;
    }

    return {
      [TRANSLATION_TAB]: viewFactory(TranslateSequenceUI),
      [AXOLABS_TAB]: viewFactory(AxolabsUI),
      [DUPLEX_TAB]: viewFactory(DuplexUI),
    }
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

export class TranslateSequenceUI extends AppUI {
  constructor(view: DG.View) {
    super(view,  'Translate Sequence');

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

export class AxolabsUI extends AppUI {
  constructor(view: DG.View) {
    super(view, 'Sequence Design');
    this.ui = new AxolabsTabUI();
  }
  private readonly ui: AxolabsTabUI;
  protected getHtml(): Promise<HTMLDivElement> {
    return Promise.resolve(this.ui.htmlDivElement);
  }
}

export class DuplexUI extends AppUI {
  constructor(view: DG.View) {
    super(view,  'Visualize Structure')
    this.ui = new SdfTabUI();
  }
  private readonly ui: SdfTabUI;

  protected getHtml(): Promise<HTMLDivElement> {
    return this.ui.getHtmlDivElement();
  }
}

class TabLayout {
  constructor(
    readonly mainTab: MainTabUI,
    readonly axolabsTab: AxolabsTabUI,
    readonly sdfTab: SdfTabUI,
    private readonly urlRouter: URLRouter
  ) {}

  private control: DG.TabControl | null = null;

  async getControl(): Promise<DG.TabControl> {
    if (this.control)
      return this.control;
    const control = ui.tabControl({
      [TRANSLATION_TAB]: await this.mainTab.getHtmlElement(),
      [AXOLABS_TAB]: this.axolabsTab.htmlDivElement,
      [DUPLEX_TAB]: await this.sdfTab.getHtmlDivElement(),
    });

    const sdfPane = control.getPane(DUPLEX_TAB);
    ui.tooltip.bind(sdfPane.header, 'Get atomic-level structure for SS + AS/AS2 and save SDF');

    const mainPane = control.getPane(TRANSLATION_TAB);
    ui.tooltip.bind(mainPane.header, 'Translate across formats');

    const axolabsPane = control.getPane(AXOLABS_TAB);
    ui.tooltip.bind(axolabsPane.header, 'Create modification pattern for SS and AS');

    control.onTabChanged.subscribe(() => {
      if (control.currentPane.name !== TRANSLATION_TAB)
        this.urlRouter.searchParams.delete('seq');
      else {
        console.log('sequence:', this.mainTab.sequence);
        this.urlRouter.searchParams.set('seq', this.mainTab.sequence);
        console.log('searchParams:', Object.entries(this.urlRouter.searchParams));
      }
      this.urlRouter.updatePath(control);
    });

    this.control = control;

    return this.control;
  }
}

class URLRouter {
  constructor(private readonly view: DG.View) {
    this.pathParts = window.location.pathname.split('/');
    this.searchParams = new URLSearchParams(window.location.search);
  }

  searchParams: URLSearchParams;
  private pathParts: string[];

  get urlParamsString(): string {
    return Object.entries(this.searchParams)
      .map(([key, value]) => `${key}=${encodeURIComponent(value)}`).join('&');
  }

  updatePath(control: DG.TabControl) {
    const urlParamsTxt: string = Object.entries(this.searchParams)
      .map(([key, value]) => `${key}=${encodeURIComponent(value)}`).join('&');
    this.view.path = '/apps/SequenceTranslator' + `/${control.currentPane.name}`;
    if (urlParamsTxt)
      this.view.path += `/?${urlParamsTxt}`;
  }

  get tabName(): string {
    const idx = this.pathParts.findIndex((el) => el === 'SequenceTranslator');
    if (idx === -1) // todo: remove after verification of validity condition
      return '';
    return this.pathParts[idx + 1];
  }
}
