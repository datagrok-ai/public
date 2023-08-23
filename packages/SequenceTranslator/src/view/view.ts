/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MAIN_TAB, AXOLABS_TAB, SDF_TAB} from './const/view';
import {MainTabUI} from './tabs/main';
import {SdfTabUI} from './tabs/sdf';
import {AxolabsTabUI} from './tabs/axolabs';
import {MonomerLibViewer} from './monomer-lib-viewer/viewer';

export class SequenceTranslatorUI {
  constructor() {
    this.view = DG.View.create();
    this.urlRouter = new URLRouter(this.view);
    this.view.box = true;
    this.view.name = 'Sequence Translator';

    const windows = grok.shell.windows;
    windows.showProperties = false;
    windows.showToolbox = false;
    windows.showHelp = false;

    const viewMonomerLibIcon = ui.iconFA('book', MonomerLibViewer.view, 'View monomer library');
    this.topPanel = [
      viewMonomerLibIcon,
    ];
    this.view.setRibbonPanels([this.topPanel]);

    this.tabs = new TabLayout(
      new MainTabUI(),
      new AxolabsTabUI(),
      new SdfTabUI(),
      this.urlRouter
    );
  }

  private readonly view: DG.View;
  public readonly tabs: TabLayout;
  private readonly topPanel: HTMLElement[];
  private readonly urlRouter: URLRouter;

  /** Create master layout of the app  */
  public async createLayout(): Promise<void> {
    const tabControl = await this.tabs.getControl();

    const tabName = this.urlRouter.tabName;
    if (tabName)
      tabControl.currentPane = tabControl.getPane(tabName);
    this.view.append(tabControl);

    grok.shell.addView(this.view);
  }
};

class TabLayout {
  constructor(
    public readonly mainTab: MainTabUI,
    public readonly axolabsTab: AxolabsTabUI,
    public readonly sdfTab: SdfTabUI,
    private readonly urlRouter: URLRouter
  ) {}

  private control: DG.TabControl | null = null;

  async getControl(): Promise<DG.TabControl> {
    if (this.control)
      return this.control;
    const control = ui.tabControl({
      [MAIN_TAB]: await this.mainTab.getHtmlElement(),
      [AXOLABS_TAB]: this.axolabsTab.htmlDivElement,
      [SDF_TAB]: await this.sdfTab.getHtmlDivElement(),
    });

    const sdfPane = control.getPane(SDF_TAB);
    ui.tooltip.bind(sdfPane.header, 'Get atomic-level structure for SS + AS/AS2 and save SDF');

    const mainPane = control.getPane(MAIN_TAB);
    ui.tooltip.bind(mainPane.header, 'Translate across formats');

    const axolabsPane = control.getPane(AXOLABS_TAB);
    ui.tooltip.bind(axolabsPane.header, 'Create modification pattern for SS and AS');

    control.onTabChanged.subscribe(() => {
      if (control.currentPane.name !== MAIN_TAB)
        this.urlRouter.searchParams.delete('seq');
      else {
        this.urlRouter.searchParams.set('seq', this.mainTab.sequence);
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

  public searchParams: URLSearchParams;
  private pathParts: string[];

  get urlParamsString(): string {
    return Object.entries(this.searchParams)
      .map(([key, value]) => `${key}=${encodeURIComponent(value)}`).join('&');
  }

  public updatePath(control: DG.TabControl) {
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
