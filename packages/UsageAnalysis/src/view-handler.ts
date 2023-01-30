import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { UaView } from './views/ua-view';
import { UaToolbox } from './ua-toolbox';
import { OverviewView } from './views/overview-view';
import { EventsView } from './views/events-view';
import { ErrorsView } from './views/errors-view';
import { FunctionErrorsView } from './views/function-errors-view';
import { UsersView } from './views/users-view';
import { DataView } from './views/data-view';

const APP_PREFIX: string = `/apps/UsageAnalysis/`;
export class ViewHandler {
  private static instance: ViewHandler;
  private urlParams: Map<string, string>  = new Map<string, string>();

  public static getInstance(): ViewHandler {
      if (!ViewHandler.instance) {
        ViewHandler.instance = new ViewHandler();
      }

      return ViewHandler.instance;
  }

  private constructor() { }

  async init() {
    let pathSplits = decodeURI(window.location.pathname).split('/');
    let params = this.getSearchParameters();

    let toolbox = await UaToolbox.construct();

    const tabs = ui.tabControl();
    tabs.root.style.width = 'inherit';
    tabs.root.style.height = 'inherit';

    let viewConstructors = [OverviewView, EventsView, ErrorsView, FunctionErrorsView, UsersView, DataView];

    const views = ['Overview', 'Events', 'Errors', 'Functions', 'Users', 'Data'];
  
    for (let i = 0; i < viewConstructors.length; ++i) {
      tabs.addPane(views[i], () => {
        const currentView = new viewConstructors[i](toolbox);
        currentView.tryToinitViewers();
        return currentView.root;
      });
    }

    if (pathSplits.length > 3 &&  pathSplits[3] != '') {
      let viewName = pathSplits[3];

      if (tabs.panes.map(p => p.name).includes(viewName))
        tabs.currentPane = tabs.getPane(viewName);
      else
        tabs.currentPane = tabs.getPane(views[0]);
    } else
      tabs.currentPane = tabs.getPane(views[0]);

    let paramsHaveDate = params.has('date');
    let paramsHaveUsers = params.has('users');
    if (paramsHaveDate || paramsHaveUsers) {
      if (paramsHaveDate) {
        //@ts-ignore
        toolbox.setDate(params.get('date'));
      }

      if (paramsHaveUsers) {
        //@ts-ignore
        await toolbox.setUsers(params.get('users'));
      }

      await toolbox.applyFilter();
    }
  
    grok.shell.newView('Usage Analysis', [tabs]);
  }

  getSearchParameters() : Map<string, string> {
    var prmstr = window.location.search.substr(1);
    return new Map<string, string>(Object.entries(prmstr != null && prmstr != "" ? this.transformToAssocArray(prmstr) : {}));
  }

  transformToAssocArray( prmstr: string ) {
    var params: {[key: string]: string} = {};
    var prmarr = prmstr.split("&");
    for ( var i = 0; i < prmarr.length; i++) {
        var tmparr = prmarr[i].split("=");
        params[decodeURI(tmparr[0])] = decodeURI(tmparr[1]);
    }
    return params;
  }

  setUrlParam(key: string, value: string, saveDuringChangingView: boolean = false) {
    let urlParams = new URLSearchParams(window.location.search);
    urlParams.set(key, value);

    let params = [];

    for (let keyAndValue of urlParams.entries())
      params.push(encodeURI(keyAndValue.join('=')));

    if (saveDuringChangingView)
      this.urlParams.set(key, value);

    grok.shell.v.path = `${APP_PREFIX}${grok.shell.v.name}?${params.join('&')}`;
  }


}