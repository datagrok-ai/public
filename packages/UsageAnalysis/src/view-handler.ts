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

  public static getInstance(): ViewHandler {
      if (!ViewHandler.instance) {
        ViewHandler.instance = new ViewHandler();
      }

      return ViewHandler.instance;
  }

  private constructor() {
    let pathSplits = decodeURI(window.location.pathname).split('/');
    let params = this.getSearchParameters();

    let toolbox = new UaToolbox();

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

    let onCurrentViewChanged = (view: ViewBase) => {
      if (view instanceof UaView) {
        view.tryToinitViewers();
        view.path = `${APP_PREFIX}${grok.shell.v.name}`;
      }
    };

    if (pathSplits.length > 3 &&  pathSplits[3] != '') {
      let viewName = pathSplits[3];

      for (let v of views) {
        if (v.name === viewName) {
          grok.shell.v = v;
          v.handleUrlParams(params);
          break;
        }
      }
    }
    else
      grok.shell.v = overviewView;

    onCurrentViewChanged(grok.shell.v);

    grok.events.onEvent('d4-current-view-changed').subscribe(() => onCurrentViewChanged(grok.shell.v));
  }

  getSearchParameters() {
    var prmstr = window.location.search.substr(1);
    return prmstr != null && prmstr != "" ? this.transformToAssocArray(prmstr) : {};
  }

  transformToAssocArray( prmstr: string ) {
    var params: {[key: string]: string} = {};
    var prmarr = prmstr.split("&");
    for ( var i = 0; i < prmarr.length; i++) {
        var tmparr = prmarr[i].split("=");
        params[tmparr[0]] = tmparr[1];
    }
    return params;
  }

  setUrlParam(key: string, value: string) {
    let urlParams = new URLSearchParams(window.location.search);
    urlParams.set(key, value);

    let params = [];

    for (let keyAndValue of urlParams.entries())
      params.push(keyAndValue.join('='));

    grok.shell.v.path = `${APP_PREFIX}${grok.shell.v.name}?${encodeURI(params.join('&'))}`;
  }


}