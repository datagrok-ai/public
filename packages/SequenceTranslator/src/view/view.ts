/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MAIN_TAB, AXOLABS_TAB, SDF_TAB} from './view-const';

import {MainTabUI} from './main-tab-ui';
import {SdfTabUI} from './sdf-tab-ui';
import {AxolabsTabUI} from './axolabs-tab-ui';

/** Class responsible for the UI of the application */
export class SequenceTranslatorUI {
  constructor() {
    this._router = new URLRouter();

    this._view = grok.shell.newView('Sequence Translator', []);
    this._view.box = true;

    // todo: should this code be here?
    const windows = grok.shell.windows;
    windows.showProperties = false;
    windows.showToolbox = false;
    windows.showHelp = false;

    // const mainTab = new MainTabUI(
    //   (seq) => {
    //     mainTab.inputSequence = seq;
    //     urlParams = new URLSearchParams();
    //     urlParams.set('seq', mainSeq);
    //     this.updatePath();
    //   });

    this._tabs = new TabLayout(
      new MainTabUI((seq) => {}), new AxolabsTabUI, new SdfTabUI
    );

    this._tabs.onTabChanged.subscribe(() => {
      // if (this._tabs.currentTab !== MAIN_TAB)
      //   urlParams.delete('seq');
      // else
      //   urlParams.set('seq', this.inputSequence);
      // this.updatePath();
    });

    this._view.append(this._tabs.control);
  }

  /** The master view containing app's main interface elements */
  private readonly _view: DG.View;
  /** Tabs control of the master view */
  private readonly _tabs: TabLayout;
  private readonly _router: URLRouter;

  get view() { return this._view; };

  private updatePath() {
    this._view.path = this._router.urlParamsString;
  }
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
