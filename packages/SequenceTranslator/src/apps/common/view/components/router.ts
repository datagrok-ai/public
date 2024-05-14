/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

class MultiviewRouter {
  constructor(private view: DG.View) {
    this.pathParts = window.location.pathname.split('/');
  }

  private appName: string;
  private tabName: string;
  private searchParams: URLSearchParams;
  private pathParts: string[];

  navigate() {
    this.parseUrl();
    this.updateHistory();
    this.updateView();
  }

  parseUrl() {
    this.appName = this.pathParts[1];
    this.tabName = this.pathParts[2];
    this.searchParams = new URLSearchParams(window.location.search);
  }

  private getUrl() {
    return '/' + this.appName + '/' + this.tabName + '?' + this.searchParams.toString();
  }

  updateHistory() {
    window.history.pushState({}, '', this.getUrl());
  }

  updateView() {
    this.view.path = this.getUrl();
  }
}

// class Router {
//   constructor(private readonly view: DG.View) {
//     this.pathParts = window.location.pathname.split('/');
//     console.log(`pathParts:`, this.pathParts);
//     this.searchParams = new URLSearchParams(window.location.search);
//     console.log(`searchParams:`, this.searchParams);
//   }

//   searchParams: URLSearchParams;
//   private pathParts: string[];

//   getUrlParamsString(): string {
//     return Object.entries(this.searchParams)
//       .map(([key, value]) => `${key}=${encodeURIComponent(value)}`).join('&');
//   }

//   updatePath(control: DG.TabControl) {
//     const urlParamsTxt: string = Object.entries(this.searchParams)
//       .map(([key, value]) => `${key}=${encodeURIComponent(value)}`).join('&');
//     this.view.path = '/apps/SequenceTranslator' + `/${control.currentPane.name}`;
//     if (urlParamsTxt)
//       this.view.path += `/?${urlParamsTxt}`;
//   }

//   get tabName(): string {
//     const idx = this.pathParts.findIndex((el) => el === 'SequenceTranslator');
//     if (idx === -1) // todo: remove after verification of validity condition
//       return '';
//     return this.pathParts[idx + 1];
//   }
// }
