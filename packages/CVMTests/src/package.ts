/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export * from './package.g';
export const _package = new DG.Package();

export class PackageFunctions {
  @grok.decorators.func()
  static info() : void {
    grok.shell.info(_package.webRoot);
  }
}
