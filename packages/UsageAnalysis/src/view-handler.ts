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
import {ViewBase} from "datagrok-api/dg";

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

  private constructor() {
  }

  async init() {
    let pathSplits = decodeURI(window.location.pathname).split('/');
    let params = this.getSearchParameters();

    let toolbox = await UaToolbox.construct();

    let overviewView = new OverviewView(toolbox);
  
    let views = [
      overviewView,
      new EventsView(toolbox),
      new ErrorsView(toolbox),
      new FunctionErrorsView(toolbox),
      new UsersView(toolbox),
      new DataView(toolbox)
    ];
  
    let neededViews = [...views];
  
    for (let v of grok.shell.views) {
      if (!(v instanceof UaView))
        continue;
      
  
      for(let i = 0; i < neededViews.length; i++) {
        if (v.name === neededViews[i].name){
          neededViews.splice(i, 1);
  
          if (v instanceof OverviewView)
            overviewView = v;
        }
      }
      
    }
  
    for (let v of neededViews) {
      grok.shell.addView(v);
    }

    let onCurrentViewChanged = (view: ViewBase, urlParams: Map<string, string>) => {
      if (view instanceof UaView) {
        view.tryToinitViewers();
        view.path = `${APP_PREFIX}${grok.shell.v.name}`;
        if (urlParams.size > 0)
          view.path += `?${Array.from(urlParams, ([k, v]) => `${k}=${v}`).join('&')}`
      }
    };

    if (pathSplits.length > 3 &&  pathSplits[3] != '') {
      let viewName = pathSplits[3];

      for (let v of views) {
        if (v.name === viewName) {
          grok.shell.v = v;
          break;
        }
      }
    }
    else
      grok.shell.v = overviewView;

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

    onCurrentViewChanged(grok.shell.v, params);

    if (grok.shell.v instanceof UaView)
      grok.shell.v.handleUrlParams(params);

    grok.events.onEvent('d4-current-view-changed').subscribe(() => onCurrentViewChanged(grok.shell.v, this.urlParams));
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