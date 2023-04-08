/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MAIN_TAB, AXOLABS_TAB, SDF_TAB} from './const/view-const';
import {MainTabUI} from './tabs-ui/main';
import {SdfTabUI} from './tabs-ui/sdf';
import {AxolabsTabUI} from './tabs-ui/axolabs';
import {viewMonomerLib} from '../utils/monomer-lib-viewer';

export class SequenceTranslatorUI {
  constructor() {
    this._router = new URLRouter();

    this._view = DG.View.create();
    this._view.box = true;
    this._view.name = 'Sequence Translator';

    // todo: should this code be here?
    const windows = grok.shell.windows;
    windows.showProperties = false;
    windows.showToolbox = false;
    windows.showHelp = false;

    // top panel icons
    const viewMonomerLibIcon = ui.iconFA('book', viewMonomerLib, 'View monomer library');
    const viewHint = ui.iconFA('lightbulb', () => {}, 'About the app');
    this._topPanel = [
      viewMonomerLibIcon,
      // viewHint
    ];
    this._view.setRibbonPanels([this._topPanel]);

    this._tabs = new TabLayout(
      new MainTabUI((seq: string) => {}),
      new AxolabsTabUI(),
      SdfTabUI.init()
    );
  }

  /** Master view containing app's main interface elements */
  private readonly _view: DG.View;
  /** Control for 3 master tabs: Main, Axolabs, SDF */
  private readonly _tabs: TabLayout;
  private readonly _topPanel: HTMLElement[];
  private readonly _router: URLRouter;

  /** Create master layout of the app  */
  public async createLayout(): Promise<void> {
    const tabControl = await this._tabs.getControl();
    // at this point we should perform the manipulations with the tab control

    this._view.append(tabControl);
    grok.shell.addView(this._view);
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
    this._pathParts = window.location.pathname.split('/');
    this._searchParams = new URLSearchParams(window.location.search);
  }

  private _searchParams: URLSearchParams;
  private _pathParts: string[];

  get urlParamsString(): string {
    return Object.entries(this._searchParams)
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
