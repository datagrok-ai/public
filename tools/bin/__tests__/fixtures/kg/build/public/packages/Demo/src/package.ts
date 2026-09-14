import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() : void {
  grok.shell.info(_package.webRoot);
}

//name: Demo App
export function demoAppLegacy() : void {
  PackageFunctions.demoApp();
}

export class PackageFunctions {
  @grok.decorators.func({name: 'Demo App'})
  static demoApp() : any {
    return null;
  }
}
