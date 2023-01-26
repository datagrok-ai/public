/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MAIN_TAB, AXOLABS_TAB, SDF_TAB} from './view-const';

import {MainTabUI} from './main-tab-ui';
import {SdfTabUI} from './sdf-tab-ui';
import {AxolabsTabUI} from './axolabs-tab-ui';
import {viewMonomerLib} from '../utils/monomer-lib-viewer';

/** Class responsible for the UI of the application */
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
    // const viewHint = ui.iconFA('lightbulb', () => {}, 'hhh');
    // const appHint = ui.iconFA('circle-question', () => {}, 'About the app');
    // const info = ui.iconFA('circle-info', () => {}, 'About the app');
    this._topPanel = [viewMonomerLibIcon];
    this._view.setRibbonPanels([this._topPanel]);

    this._tabs = new TabLayout(
      new MainTabUI((seq: string) => {}), new AxolabsTabUI(), new SdfTabUI()
    );

    this._view.append(this._tabs.control);
  }

  /** Master view containing app's main interface elements */
  private readonly _view: DG.View;
  /** Control for 3 master tabs: Main, Axolabs, SDF */
  private readonly _tabs: TabLayout;
  private readonly _topPanel: HTMLElement[];
  private readonly _router: URLRouter;

  get view() { return this._view; };

  private updatePath() {
    this._view.path = this._router.urlParamsString;
  }

  public addView(): void { grok.shell.addView(this._view); }
};

class TabLayout {
  constructor(mainTab: MainTabUI, axolabsTab: AxolabsTabUI, sdfTab: SdfTabUI) {
    this._control = ui.tabControl({
      [MAIN_TAB]: mainTab.htmlDivElement,
      [AXOLABS_TAB]: axolabsTab.htmlDivElement,
      [SDF_TAB]: sdfTab.htmlDivElement,
    });

    // todo: all tab tooltips to be defined here
    const sdfPane = this._control.getPane(SDF_TAB);
    ui.tooltip.bind(sdfPane.header, 'Get atomic-level structure for SS + AS/AS2 and save SDF');
  }

  private readonly _control: DG.TabControl;
  // todo: port to main tab object
  // private inputSequence: string;

  get control(): DG.TabControl { return this._control; }

  get currentTab(): string { return this._control.currentPane.name; }

  get onTabChanged() { return this._control.onTabChanged; }
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
}
