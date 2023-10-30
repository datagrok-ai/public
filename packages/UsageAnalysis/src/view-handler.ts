import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {UaView} from './tabs/ua';
import {UaToolbox} from './ua-toolbox';
import {EventsView} from './tabs/events';
import {PackagesView} from './tabs/packages';
import {FunctionsView} from './tabs/functions';
import {OverviewView} from './tabs/overview';
import {LogView} from './tabs/log';
import {TestsView} from './tabs/tests';

const APP_PREFIX: string = `/apps/UsageAnalysis/`;


export class ViewHandler {
  private static instance: ViewHandler;
  private urlParams: Map<string, string> = new Map<string, string>();
  public static UAname = 'Usage Analysis';
  static UA: DG.MultiView;
  dockFilters: DG.DockNode | null = null;

  public static getInstance(): ViewHandler {
    if (!ViewHandler.instance)
      ViewHandler.instance = new ViewHandler();
    return ViewHandler.instance;
  }

  async init() {
    ViewHandler.UA = new DG.MultiView({viewFactories: {}});
    ViewHandler.UA.parentCall = grok.functions.getCurrentCall();
    const toolbox = await UaToolbox.construct();
    const params = this.getSearchParameters();
    // [ErrorsView, FunctionsView, UsersView, DataView];
    const viewClasses: (typeof UaView)[] = [OverviewView, PackagesView, FunctionsView, EventsView, LogView, TestsView];
    for (let i = 0; i < viewClasses.length; i++) {
      const currentView = new viewClasses[i](toolbox);
      currentView.tryToinitViewers();
      ViewHandler.UA.addView(currentView.name, () => currentView, false);
    }
    const paramsHaveDate = params.has('date');
    const paramsHaveUsers = params.has('groups');
    const paramsHavePackages = params.has('packages');
    if (paramsHaveDate || paramsHaveUsers || paramsHavePackages) {
      if (paramsHaveDate)
        toolbox.setDate(params.get('date')!);
      if (paramsHaveUsers)
        toolbox.setGroups(params.get('groups')!);
      if (paramsHavePackages)
        toolbox.setPackages(params.get('packages')!);
      toolbox.applyFilter();
    }
    let helpShown = false;
    ViewHandler.UA.tabs.onTabChanged.subscribe((tab) => {
      const view = ViewHandler.UA.currentView;
      ViewHandler.UA.path = ViewHandler.UA.path.replace(/(UsageAnalysis\/)([a-zA-Z]+)/, '$1' + view.name);
      if (view instanceof UaView) {
        for (const viewer of view.viewers) {
          if (!viewer.activated) {
            viewer.activated = true;
            viewer.reloadViewer();
          }
        }
      }
      if (!helpShown) {
        if (ViewHandler.UA.currentView instanceof PackagesView || ViewHandler.UA.currentView instanceof FunctionsView) {
          grok.shell.windows.showToolbox = true;
          grok.shell.windows.showContextPanel = true;
          const info = ui.divText(`To learn more about an event, click the corresponding point.\
          To select multiple points, use CTRL + Click or SHIFT + Mouse Drag. Once you've made your selection,\
          see the detailed information on the Context Panel`);
          info.classList.add('ua-hint');
          grok.shell.o = info;
        }
        helpShown = true;
      }
      grok.shell.windows.showToolbox = view.name !== 'Tests';
    });
    ViewHandler.UA.name = ViewHandler.UAname;
    ViewHandler.UA.box = true;
    const urlTab = window.location.pathname.match(/UsageAnalysis\/([a-zA-Z]+)/)?.[1];
    ViewHandler.UA.path = APP_PREFIX + (urlTab ?? 'Overview');
    if (urlTab) ViewHandler.changeTab(urlTab);
    grok.shell.addView(ViewHandler.UA);
  }

  public static getView(name: string) {
    return ViewHandler.UA.getView(name) as UaView;
  }

  public static getCurrentView(): UaView {
    return ViewHandler.UA.currentView as UaView;
  }

  public static changeTab(name: string) {
    ViewHandler.UA.tabs.currentPane = ViewHandler.UA.tabs.getPane(name);
  }

  getSearchParameters() : Map<string, string> {
    const prmstr = window.location.search.substring(1);
    return new Map<string, string>(Object.entries(prmstr ? this.transformToAssocArray(prmstr) : {}));
  }

  transformToAssocArray(prmstr: string) {
    const params: {[key: string]: string} = {};
    const prmarr = prmstr.split('&');
    for (let i = 0; i < prmarr.length; i++) {
      const tmparr = prmarr[i].split('=');
      params[decodeURI(tmparr[0])] = decodeURI(tmparr[1]);
    }
    return params;
  }

  setUrlParam(key: string, value: string, saveDuringChangingView: boolean = false) {
    const urlParams = new URLSearchParams(window.location.search);
    urlParams.set(key, value);

    const params: string[] = [];

    for (const keyAndValue of urlParams.entries())
      params.push(encodeURI(keyAndValue.join('=')));

    if (saveDuringChangingView)
      this.urlParams.set(key, value);

    ViewHandler.UA.path = `${APP_PREFIX}${ViewHandler.getCurrentView().name}?${params.join('&')}`;
  }
}
