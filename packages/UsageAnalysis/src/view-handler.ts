import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {UaView} from './views/ua-view';
import {UaToolbox} from './ua-toolbox';
import {OverviewView} from './views/overview-view';
import {EventsView} from './views/events-view';
import {ErrorsView} from './views/errors-view';
import {FunctionsView} from './views/function-errors-view';
import {UsersView} from './views/users-view';
import {DataView} from './views/data-view';
import {PackagesView} from './views/packages-view';

const APP_PREFIX: string = `/apps/UsageAnalysis/`;
export class ViewHandler {
  private static instance: ViewHandler;
  private urlParams: Map<string, string> = new Map<string, string>();

  public static getInstance(): ViewHandler {
    if (!ViewHandler.instance)
      ViewHandler.instance = new ViewHandler();


    return ViewHandler.instance;
  }

  private constructor() { }

  async init() {
    const pathSplits = decodeURI(window.location.pathname).split('/');
    const params = this.getSearchParameters();

    const toolbox = await UaToolbox.construct();

    const tabs = ui.tabControl();
    tabs.root.style.width = 'inherit';
    tabs.root.style.height = 'inherit';

    const viewClasses: (typeof UaView)[] = [OverviewView, EventsView,
      ErrorsView, FunctionsView, UsersView, DataView, PackagesView];

    for (let i = 0; i < viewClasses.length; ++i) {
      tabs.addPane(viewClasses[i].viewName, () => {
        const currentView = new viewClasses[i](toolbox);
        currentView.tryToinitViewers();
        return currentView.root;
      });
    }

    if (pathSplits.length > 3 && pathSplits[3] != '') {
      const viewName = pathSplits[3];

      if (tabs.panes.map((p) => p.name).includes(viewName))
        tabs.currentPane = tabs.getPane(viewName);
      else
        tabs.currentPane = tabs.getPane(viewClasses[0].viewName);
    } else
      tabs.currentPane = tabs.getPane(viewClasses[0].viewName);

    const paramsHaveDate = params.has('date');
    const paramsHaveUsers = params.has('users');
    if (paramsHaveDate || paramsHaveUsers) {
      if (paramsHaveDate)
        toolbox.setDate(params.get('date')!);


      if (paramsHaveUsers)
        await toolbox.setUsers(params.get('users')!);


      await toolbox.applyFilter();
    }

    grok.shell.newView('Usage Analysis', [tabs]);
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

    grok.shell.v.path = `${APP_PREFIX}${grok.shell.v.name}?${params.join('&')}`;
  }
}
