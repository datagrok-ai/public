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
import {TestsView, filters} from './tabs/tests';
import {ReportsView} from "./tabs/reports";
import {ErrorsView} from "./tabs/errors"

export class ViewHandler {
  private static instance: ViewHandler;
  private urlParams: Map<string, string> = new Map<string, string>();
  public static UAname = 'Usage Analysis';
  static UA: DG.MultiView;

  public static getInstance(): ViewHandler {
    if (!ViewHandler.instance)
      ViewHandler.instance = new ViewHandler();
    return ViewHandler.instance;
  }

  async init(params: {date?: string, groups?: string, packages?: string}, path?: string): Promise<void> {
    ViewHandler.UA = new DG.MultiView({viewFactories: {}});
    ViewHandler.UA.parentCall = grok.functions.getCurrentCall();
    const toolbox = await UaToolbox.construct();
    toolbox.filters.root.after(filters);
    // [ErrorsView, FunctionsView, UsersView, DataView];
    const viewClasses: (typeof UaView)[] = [OverviewView, PackagesView, FunctionsView, EventsView, LogView, TestsView, ErrorsView, ReportsView];
    const inits: Promise<void>[] = [];
    for (let i = 0; i < viewClasses.length; i++) {
      const currentView = new viewClasses[i](toolbox);
      inits.push(currentView.tryToinitViewers(path));
      ViewHandler.UA.addView(currentView.name, () => currentView, false);
    }
    await Promise.all(inits);
    const paramsHaveDate = params.date != undefined;
    const paramsHaveUsers = params.groups != undefined;
    const paramsHavePackages = params.packages != undefined;
    if (paramsHaveDate || paramsHaveUsers || paramsHavePackages) {
      if (paramsHaveDate)
        toolbox.setDate(params.date!);
      if (paramsHaveUsers)
        toolbox.setGroups(params.groups!);
      if (paramsHavePackages)
        toolbox.setPackages(params.packages!);
      toolbox.applyFilter();
    }
    let helpShown = false;
    const puButton = ui.bigButton('Usage', () => {
      const v = ViewHandler.getCurrentView();
      v.switchRout();
      ViewHandler.updatePath();
      v.viewers[1].root.style.display = 'none';
      v.viewers[0].root.style.display = 'flex';
      puButton.disabled = true;
      piButton.disabled = false;
    });
    const piButton = ui.bigButton('Installation time', () => {
      const v = ViewHandler.getCurrentView();
      v.switchRout();
      ViewHandler.updatePath();
      v.viewers[0].root.style.display = 'none';
      v.viewers[1].root.style.display = 'flex';
      puButton.disabled = false;
      piButton.disabled = true;
    });
    puButton.disabled = true;
    const pButtons = ui.divH([puButton, piButton], 'ua-packages-buttons');
    pButtons.style.display = 'none';
    toolbox.filters.root.before(pButtons);

    const fuButton = ui.bigButton('Usage', () => {
      const v = ViewHandler.getCurrentView() as FunctionsView;
      v.switchRout();
      ViewHandler.updatePath();
      v.functionsExecTime.style.display = 'none';
      v.viewers[0].root.style.display = 'flex';
      fuButton.disabled = true;
      feButton.disabled = false;
    });
    const feButton = ui.bigButton('Execution time', () => {
      const v = ViewHandler.getCurrentView() as FunctionsView;
      v.switchRout();
      ViewHandler.updatePath();
      v.viewers[0].root.style.display = 'none';
      v.functionsExecTime.style.display = 'flex';
      fuButton.disabled = false;
      feButton.disabled = true;
    });
    fuButton.disabled = true;
    const fButtons = ui.divH([fuButton, feButton], 'ua-packages-buttons');
    fButtons.style.display = 'none';
    toolbox.filters.root.before(fButtons);

    ViewHandler.UA.tabs.onTabChanged.subscribe((tab) => {
      const view = ViewHandler.UA.currentView;
      // ViewHandler.UA.path = ViewHandler.UA.path.replace(/(UsageAnalysis\/)([a-zA-Z/]+)/, '$1' + view.name);
      ViewHandler.updatePath();
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
      if (view.name === 'Tests') {
        toolbox.filters.expanded = false;
        filters.style.display = 'flex';
      } else {
        toolbox.filters.expanded = true;
        filters.style.display = 'none';
      }
      if (view.name === 'Packages')
        pButtons.style.display = 'flex';
      else
        pButtons.style.display = 'none';
      if (view.name === 'Functions')
        fButtons.style.display = 'flex';
      else
        fButtons.style.display = 'none';
    });
    ViewHandler.UA.name = ViewHandler.UAname;
    ViewHandler.UA.box = true;
    let urlTab = 'Overview';

    if (path != undefined && path.length > 1) {
      const segments = path.split('/').filter((s) => s != '');
      if (segments.length > 0) {
        urlTab = segments[0];
        urlTab = urlTab[0].toUpperCase() + urlTab.slice(1);
      }
    }
    if (viewClasses.some((v) => v.name === `${urlTab}View`))
      ViewHandler.changeTab(urlTab);
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

    ViewHandler.UA.path = `/${ViewHandler.getCurrentView().name}?${params.join('&')}`.toLowerCase();
  }

  static updatePath(): void {
    const v = ViewHandler.getCurrentView();
    const s = ViewHandler.UA.path.split('?');
    const params = s.length === 2 ? s[1] : null;
    ViewHandler.UA.path = `/${v.name}${v.rout ?? ''}${params ? '?' + params : ''}`.toLowerCase();
  }
}
