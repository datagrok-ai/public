/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MAIN_TAB, AXOLABS_TAB, SDF_TAB} from './const/view';
import {MainTabUI} from './tabs/main';
import {SdfTabUI} from './tabs/sdf';
import {AxolabsTabUI} from './tabs/axolabs';
import {viewMonomerLib} from '../model/monomer-lib-viewer';

export class SequenceTranslatorUI {
  constructor() {
    this.router = new URLRouter();

    this.view = DG.View.create();
    this.view.box = true;
    this.view.name = 'Sequence Translator';

    // todo: should this code be here?
    const windows = grok.shell.windows;
    windows.showProperties = false;
    windows.showToolbox = false;
    windows.showHelp = false;

    // top panel icons
    const viewMonomerLibIcon = ui.iconFA('book', viewMonomerLib, 'View monomer library');
    const viewHint = ui.iconFA('lightbulb', () => {}, 'About the app');
    this.topPanel = [
      viewMonomerLibIcon,
      // viewHint
    ];
    this.view.setRibbonPanels([this.topPanel]);

    this.tabs = new TabLayout(
      new MainTabUI(),
      new AxolabsTabUI(),
      new SdfTabUI()
    );
  }

  /** Master view containing app's main interface elements */
  private readonly view: DG.View;
  /** Control for 3 master tabs: Main, Axolabs, SDF */
  private readonly tabs: TabLayout;
  private readonly topPanel: HTMLElement[];
  private readonly router: URLRouter;

  /** Create master layout of the app  */
  public async createLayout(): Promise<void> {
    const tabControl = await this.tabs.getControl();
    // at this point we should perform the manipulations with the tab control

    this.view.append(tabControl);
    grok.shell.addView(this.view);
  }
};

class TabLayout {
  constructor(
    private readonly mainTab: MainTabUI,
    private readonly axolabsTab: AxolabsTabUI,
    private readonly sdfTab: SdfTabUI
  ) {}

  async getControl(): Promise<DG.TabControl> {
    const control = ui.tabControl({
      [MAIN_TAB]: await this.mainTab.getHtmlElement(),
      [AXOLABS_TAB]: this.axolabsTab.htmlDivElement,
      [SDF_TAB]: this.sdfTab.htmlDivElement,
    });

    // bind tooltips to each tab
    // todo: bind tooltips to Main and Axolabs
    const sdfPane = control.getPane(SDF_TAB);
    ui.tooltip.bind(sdfPane.header, 'Get atomic-level structure for SS + AS/AS2 and save SDF');

    // control.onTabChanged.subscribe(() => {
    //   if (control.currentPane.name !== MAIN_TAB)
    //     urlParams.delete('seq');
    //   else
    //     urlParams.set('seq', mainSeq);
    //   updatePath();
    // });

    return control;
  }
}

class URLRouter {
  constructor() {
    this.pathParts = window.location.pathname.split('/');
    this.searchParams = new URLSearchParams(window.location.search);
  }

  private searchParams: URLSearchParams;
  private pathParts: string[];

  get urlParamsString(): string {
    return Object.entries(this.searchParams)
      .map(
        ([key, value]) => `${key}=${encodeURIComponent(value)}`
      ).join('&');
  }

  // public updatePath(control: DG.TabControl) {
  //   const urlParamsTxt: string = Object.entries(this._searchParams)
  //     .map(([key, value]) => `${key}=${encodeURIComponent(value)}`).join('&');
  //   view.path = '/apps/SequenceTranslator' + `/${control.currentPane.name}/?${urlParamsTxt}`;
  // }

  // if (this._pathParts.length >= 4) {
  //   const tabName: string = pathParts[3];
  //   tabControl.currentPane = tabControl.getPane(tabName);
  // }
}
